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
   def __init__(self, serverpool, port, verbose, queries, dbtype, maxhits, maxevalue, dbnames, outfile, errors):
      threading.Thread.__init__(self)
      self.verbose = verbose
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
            nhits = statusvalues[14]
            dbsize = statusvalues[3]
            if self.verbose:
              sys.stderr.write("Models: %i  Hits: %i  Dbsize: %i\n" % (nmodels, nhits, dbsize))
            if nmodels != len(self.dbnames):
              sys.stderr.write("FATAL: Number of models in database %i and number of database names from names file %i are different. Aborting.\n" % (nmodels, len(self.dbnames)))
              sys.exit(1)
            selectedhits = nhits
            if self.maxhits:
              selectedhits=min(self.maxhits, nhits)
            for i in range(selectedhits):
              offset=120+i*152
              hitvalues=unpack('3QI4xd3f4x3df9I4Q', data[offset:offset+152])
              dbid=hitvalues[0]
              dbname=self.dbnames[dbid-1]
              evalue=exp(hitvalues[8]) * dbsize
              if self.verbose:
                sys.stderr.write("Query: %s Model: %i %s  Evalue: %1.1e\n" % (queryname, dbid, dbname, evalue))
              if evalue <= self.maxevalue:
                threadLock.acquire()
                self.outfile.write("%s\t%s\t%1.1e\n" % (queryname, dbname, evalue))
                threadLock.release()
            success=True
        threadLock.acquire()
        self.serverpool.release_server(server)
        threadLock.release()

        # Close connection
        if self.verbose:
          sys.stderr.write("Closing connection.\n")
        s.close()

def startnewthread(args, queries, dbnames, errors, serverpool):
  thread = hmmpgmdThread(serverpool, args.port, args.verbose, queries, args.dbtype, args.maxhits, args.maxevalue, dbnames, args.outfile, errors)
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
c=0
while threading.activeCount() > 1 and c < args.finishtimeout:
  time.sleep(1)
  c += 1

if errors:
  for error in errors:
    sys.stderr.write(error)
  sys.exit(1)
