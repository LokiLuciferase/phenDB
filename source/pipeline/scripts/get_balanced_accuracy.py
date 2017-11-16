#/user/bin/python
import sys
#
#contimantion to the right, completeness to the bottom

ma=[]
with open(sys.argv[1],'r') as myfile:
    for line in myfile:
        if "[completeness]" in line:
            break
        if "mean" not in line:
            ma.append(line.split())

with open(sys.argv[2],'r') as myfile:
    for line in myfile:
        linelist=line.split()
        comp=linelist[0]
        cont=linelist[1]



cont_steps=int(round(float(cont)/0.05))
completeness_steps=int(round(float(comp)/0.05))


print(ma[completeness_steps][cont_steps])
