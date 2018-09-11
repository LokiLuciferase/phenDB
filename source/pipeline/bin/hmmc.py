#!/usr/bin/env python3

"""hmmc.py: Client for hmmpgmd servers (part of HMMER3)."""

__author__      = "Thomas Rattei"
__copyright__   = "Copyright 2017, University of Vienna"

import sys, os.path, argparse, socket, gzip, threading, time
from Bio import SeqIO
from struct import unpack
from math import exp
from operator import itemgetter

## Argument handling
parser = argparse.ArgumentParser(description='Client for hmmpgmd servers')

parser.add_argument('-v', '--verbose', action="store_true",
   help="verbose output")

parser.add_argument('-i', '--infile', type=argparse.FileType('r'), default=sys.stdin,
   help="optional input multiple protein FASTA file (default is stdin)")

parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout,
   help="optional output file (tab separated values: protein name, match name, evalue; default is stdout)")

parser.add_argument('-d', '--dbfile', type=argparse.FileType('r'), required=True,
   help="database file (sequence or HMM file, see hmmpgmd docs)")

parser.add_argument('-t', '--dbtype', default="h", choices=['h', 's'],
   help="database type (h for HMM [default], s for sequence, see hmmpgmd docs)")

parser.add_argument('-n', '--numthreads', default=1, type=int,
   help="number of parallel queries")

parser.add_argument('-q', '--queriesperthread', default=10, type=int,
   help="queries per thread (default=10)")

parser.add_argument('-s', '--servers', nargs = '*', default="127.0.0.1",
   help="TCPIP address of one or multiple servers (blank delimited, default=127.0.0.1)")

parser.add_argument('-p', '--port', default=51371, type=int,
   help="TCP port of the server (default=51371)")

parser.add_argument('-f', '--finishtimeout', default=60, type=int,
   help="maximal seconds to finish (default=60)")

parser.add_argument('-e', '--maxevalue', default=0.05, type=float,
   help="Maximal evalue to report (default=0.05)")

parser.add_argument('-m', '--maxhits', default=0, type=int,
   help="Maximal number of hits to report (default=0, which means no limit)")

parser.add_argument('-y', '--modeY', action='store_true',
   help="Downstream processing allowing multiple non-overlapping hits per sequence")

args = parser.parse_args()


## Functions and classes
def init_dbnames(dbfilename, verbose):
  namesfilename = "%s.names.gz" % dbfilename
  names=[]
  if not os.path.isfile(namesfilename):
    sys.stderr.write("Names file %s does not exist. Please extract all HMM or sequence names, e.g. by running:\n" % namesfilename)
    sys.stderr.write("grep '^ACC   ' %s | cut -c7- | gzip >%s" % (dbfilename, namesfilename))
    sys.exit(1)
  if verbose:
    sys.stderr.write("Reading names filename %s..." % namesfilename)
  with gzip.open(namesfilename, "rt") as infile:
    for line in infile:
      names.append(line.strip())
  if verbose:
    sys.stderr.write("done. Retrieved %i names.\n" % len(names))
  return names


def send_hmmpgmd(s, dbtype, queryname, querysequence):
  # Prepare query
  if dbtype == "s":
    requestheader = "@--seqdb 1"
  else:
    requestheader = "@--hmmdb 1"
  request = "%s\n>%s\n%s\n//" % (requestheader, queryname, querysequence)
  # Perform query
  s.send(bytes(request, encoding="UTF-8"))


def receive_hmmpgmd(s):
  (status, msglen) = unpack('I4xQ', s.recv(16, socket.MSG_WAITALL))
  if status:
    errmsg = unpack('%is' % msglen, s.recv(msglen, socket.MSG_WAITALL))
    sys.stderr.write("WARNING (trying another server with the same data): There was an error reported by HMMER: %s" % (errmsg[0].decode("utf-8", "ignore")))
    return None
  else:
    return s.recv(msglen, socket.MSG_WAITALL)

class serverPool:
   def __init__(self, servers):
      self.connections_by_server={}
      for server in servers:
        self.connections_by_server[server]=0

   def report_down(self, server):
      if server in self.connections_by_server:
        sys.stderr.write("WARNING: Server %s not responding. Trying other servers.\n" % (server))
        self.connections_by_server.pop(server)
      if not self.connections_by_server:
        sys.stderr.write("FATAL: No more available servers. Aborting.\n")
        sys.exit(1)

   def get_server(self):
      if not self.connections_by_server:
        return None
      (server, connections) = sorted(self.connections_by_server.items(), key=itemgetter(1))[0]
      self.connections_by_server[server] += 1
      return server

   def release_server(self, server):
      if server in self.connections_by_server:
        self.connections_by_server[server] -= 1

class hmmpgmdThread (threading.Thread):
   def __init__(self, serverpool, port, verbose, queries, dbtype, maxhits, maxevalue, dbnames, outfile, errors, modeY):
      threading.Thread.__init__(self)
      self.verbose = verbose
      self.modeY = modeY
      self.serverpool = serverpool
      self.port = port
      self.queries = queries
      self.dbtype = dbtype
      self.maxhits = maxhits
      self.maxevalue = maxevalue
      self.dbnames = dbnames
      self.outfile = outfile
      self.errors = errors

   def run(self):
      success=False
      while not success:
        threadLock.acquire()
        server = self.serverpool.get_server()
        threadLock.release()
        if not server:
          self.errors.append("Did not get server from serverpool. Aborting.\n")
          break

        try:
          # Open connection
          if self.verbose:
            sys.stderr.write("Initializing connection to host %s at port %i.\n" % (server, self.port))
          s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
          s.connect((server, self.port))
        except:
          threadLock.acquire()
          self.serverpool.report_down(server)
          threadLock.release()
          continue

        for (queryname, querysequence) in self.queries:
          send_hmmpgmd(s, self.dbtype, queryname, querysequence)
          data = receive_hmmpgmd(s)
          if not data:
            threadLock.acquire()
            self.serverpool.report_down(server)
            threadLock.release()
            success=False
            break
          else:
            if self.verbose:
              sys.stderr.write("Received: %i bytes.\n" % len(data))
            # Unpack and extract result data
            statusvalues=unpack('5d2I9q', data[:120])
            nmodels = statusvalues[7]
            nhits = statusvalues[13]
            nreported = statusvalues[14]
            dbsize = statusvalues[3]
            dom_dbsize = statusvalues[4]
            if self.verbose:
              sys.stderr.write("Models: %i  Hits: %i  Dbsize: %i\n" % (nmodels, nhits, dbsize))
            if nmodels != len(self.dbnames):
              sys.stderr.write("FATAL: Number of models in database %i and number of database names from names file %i are different. Aborting.\n" % (nmodels, len(self.dbnames)))
              sys.exit(1)
            if not len(statusvalues) == 16:
              sys.stderr.write("FATAL: Missmatch between the number of stats data elements recieved [%i] and expected [%i]" % (len(statusvalues),16))
              sys.exit(1)
            
            selectedhits = nhits
            if self.maxhits:
              selectedhits=min(self.maxhits, nhits)
            hits=[]
            for i in range(nhits):
              offset=120+i*152
              hitvalues=unpack('3QI4xd3f4x3df9I4Q', data[offset:offset+152])
              dbid=hitvalues[0]
              dbname=self.dbnames[dbid-1]
              evalue=exp(hitvalues[8]) * dbsize
              score=round(hitvalues[5],1)
              ndom=hitvalues[16]
              nrep=hitvalues[18]
              hits.append({'queryname':queryname,'dbid':dbid,'dbname':dbname,'evalue':evalue,'ndom':ndom,'nrep':nrep,'score':score})
            if not self.modeY:
              for i in range(selectedhits):
                dbid=hits[i]['dbid']
                dbname=self.dbnames[dbid-1]
                evalue=hits[i]['evalue']
                if self.verbose:
                  sys.stderr.write("Query: %s Model: %i %s  Evalue: %1.1e\n" % (queryname, dbid, dbname, evalue))
                if evalue <= self.maxevalue:
                  threadLock.acquire()
                  self.outfile.write("%s\t%s\t%1.1e\n" % (queryname, dbname, evalue))
                  threadLock.release()
            else:
              hits_list=[]
              all_domains=[]
              outfmt=['queryname','dbname','evalue','score','ievalue','cevalue','tlen','hmmfrom','hmmto','qlen','sqfrom','sqto']
              offset=120+nhits*152
              for hit in hits:
                hit_domains=[]
                print_next=False
                for i in range(hit['ndom']):
                  domainvalues=unpack('4i5f4xd2iQ8x', data[offset:offset+72])
                  if domainvalues[9] > sys.float_info.max_10_exp:
                    local_evalue=exp(min(domainvalues[9],sys.float_info.max_10_exp))
                    sys.stderr.write("WARNING: evalue-exponent (%i) exeeds sys.max, unpacking binary data erroneous?\n" % domainvalues[9])
                  else:
                    local_evalue=exp(domainvalues[9])
                  domain_dict=hit.copy()
                  add_dict={'ievalue':local_evalue*dbsize,'cevalue':local_evalue*dom_dbsize}
                  domain_dict.update(add_dict)
                  hit_domains.append(domain_dict)
                  offset=offset+72
                all_domains.append(hit_domains)
                for domain in hit_domains:
                  alivalues=unpack("7QI4x3Q3I4x6QI4xQ", data[offset:offset+168])
                  memsize=alivalues[20]
                  offset=offset+168+memsize
                  dict_ali={'tlen':alivalues[13],'hmmfrom':alivalues[11],'hmmto':alivalues[12],'qlen':alivalues[19],'sqfrom':alivalues[17],'sqto':alivalues[18]}
                  domain.update(dict_ali)
                  hits_list.append(domain)
                  #self.outfile.write("%s\n" % "\t".join([str(domain[key]) for key in outfmt]))

              #downstream processing for mode Y here:
              #Yanbin Yin script from 07/21/2015 adapted by
              #Patrick Hyden on 03/27/2018 to fit our needs

              # filtering hits below (adaptive) cutoff (1e-3 or 1e-5) and target coverage (0.3):
              del_ids=[]
              for i in range(0,len(hits_list)):
                hit=hits_list[i]
                if hit['ievalue']>1e-3:
                  del_ids.append(i)
                elif hit['sqto']-hit['sqfrom']>80 and hit['ievalue']>1e-5:
                  del_ids.append(i)
                #delete hits with low target coverage
                #elif float(hit['hmmto']-hit['hmmfrom'])/hit['tlen']<0.3:
                  #del sorted_hits[i]
                  #i=i-1

              #delete them
              for i in reversed(del_ids):
                del hits_list[i]

              if len(hits_list) > 1:
                sorted_hits=sorted(hits_list, key=lambda k: k['sqfrom'])
              else:
                sorted_hits=hits_list
              i = 1

              # filtering overlapping hits
              while i < len(sorted_hits):
                hit1=sorted_hits[i-1]
                hit2=sorted_hits[i]
                
                len_h1=hit1['sqto']-hit1['sqfrom']
                len_h2=hit2['sqto']-hit2['sqfrom']
                overlap=hit1['sqto']-hit2['sqfrom']

                if len_h1 == 0:
                  del sorted_hits[i-1]
                  i = i - 1
                elif len_h2 == 0:
                  del sorted_hits[i]
                  i = i - 1
                elif float(overlap)/len_h1 > 0.5 or float(overlap)/len_h2 > 0.5:
                  if hit1['ievalue']>hit2['ievalue']:
                    del sorted_hits[i-1]
                  elif hit1['ievalue']<hit2['ievalue']:
                    del sorted_hits[i]
                  elif hit1['score']<hit2['score']:
                    del sorted_hits[i-1]
                  else:
                    del sorted_hits[i]
                  i=i-1
                i=i+1
              hits_list=sorted_hits

              checked_ids={}
              del_ids=[]
              for i in range(0,len(hits_list)):
                if not i in checked_ids:
                  hit=hits_list[i]
                  hit_len=hit['hmmto']-hit['hmmfrom']
                  if hit_len/hit['tlen']<0.3:
                    #not directly delete, look for other hits with same 'dbid'
                    ids=[i]
                    for k in range(i,len(hits_list)):
                      hit_comp=hits_list[k]
                      if hit_comp['dbid']==hit['dbid']:
                        if hit_comp['hmmto']<=hit['hmmfrom'] or hit_comp['hmmfrom']>=hit['hmmto']:
                          hit_len=hit_len+hit_comp['hmmto']-hit_comp['hmmfrom']
                          ids.append(k)
                    for k in ids:
                      checked_ids[k]=1
                      if hit_len/hit['tlen']<0.3:
                        del_ids.append(k)
                      elif not k==i:
                        del_ids.append(k)

              del_ids.sort(reverse=True)
              #print(del_ids)
              for i in del_ids:
                try:
                  del_hit=hits_list[i]
                  del hits_list[i]

                except IndexError:
                  print(i, del_ids, len(hits_list))
                  print("%s" % "\t".join([str(del_hit[key]) for key in outfmt]))


              for hit in hits_list:
                self.outfile.write("%s\n" % "\t".join([str(hit[key]) for key in outfmt]))

            success=True
        threadLock.acquire()
        self.serverpool.release_server(server)
        threadLock.release()

        # Close connection
        if self.verbose:
          sys.stderr.write("Closing connection.\n")
        s.close()

def startnewthread(args, queries, dbnames, errors, serverpool):
  thread = hmmpgmdThread(serverpool, args.port, args.verbose, queries, args.dbtype, args.maxhits, args.maxevalue, dbnames, args.outfile, errors, args.modeY)
  thread.start()


### Main program

## Init database accessions
dbnames = init_dbnames(args.dbfile.name, args.verbose)

## Parse input file and submit in chunks of args.queriesperthread proteins to server
errors=[]
threadLock = threading.Lock()
serverpool = serverPool(args.servers)
queries=[]
if args.infile.name.endswith(".gz"):
  infile=gzip.open(args.infile.name, "rt")
else:
  infile=args.infile
for entry in SeqIO.parse(infile, "fasta"):
  queries.append((entry.id, entry.seq))
  if len(queries) == args.queriesperthread:
    while threading.activeCount() > args.numthreads:
      time.sleep(0.001)
    startnewthread(args, queries, dbnames, errors, serverpool)
    queries=[]
if queries:
  startnewthread(args, queries, dbnames, errors, serverpool)

## Wait for threads to complete
c = 0
while threading.activeCount() >= 2 and c < args.finishtimeout:
  time.sleep(1)
  c += 1

if errors:
  for error in errors:
    sys.stderr.write(error)
  os._exit(1)
elif c >= args.finishtimeout:
  print("Exited with timeout. {n} threads were still active.".format(n=threading.activeCount()))
  os._exit(0)
