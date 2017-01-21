#!/usr/bin/python
import sys
tani='../list/ZINC_tanimoto.dat'

with open(tani,"r") as inf:
    next(inf)
    for line in inf:
        line=line.strip().split("\t")
        qIdx=int(line[0])
        rest=line[1:]
        tIdx=int(1)
        for dis in rest:
            if tIdx > qIdx:
                sim=float(1.0-float(dis))
                if sim >= 0.5:
                    print "%s, %s, %s"%(str(qIdx),str(tIdx),str(sim))
            tIdx+=1
