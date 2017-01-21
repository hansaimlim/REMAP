#!/usr/bin/python

prots='../list/ZINC_protein_index.tsv' #list of unique proteins with protein index
blast='../list/ZINC_blast_result.dat' #list of BLAST results for ZINC proteins
idx2prot={}
prot2idx={}
with open(prots,"r") as protline:
    next(protline)
    for line in protline:
        line=line.strip().split("\t")
        idx=str(line[0])
        prot=str(line[1])
        idx2prot[idx]=prot
        prot2idx[prot]=idx

maxbitscore={}
with open(blast,"r") as bline:
    for line in bline:
        line=line.strip().split("\t")
        query=str(line[0]).strip().split("|")[1] #Query gene ID
        target=str(line[1]).strip().split("|")[1] #Target gene ID
        bitscore=float(line[-1]) #bit score for query-target pair
        if query == target:
            maxbitscore[query]=bitscore #define self-query score

with open(blast,"r") as bline:
    for line in bline:
        line=line.strip().split("\t")
        query=str(line[0]).strip().split("|")[1] #Query gene ID
        target=str(line[1]).strip().split("|")[1] #Target gene ID
        bitscore=float(line[-1]) #bit score for query-target pair
        if query == target:
            #skip self-queries
            continue
        else:
            qIdx=prot2idx[query]
            tIdx=prot2idx[target]
            simscore=float(bitscore/maxbitscore[query])
            print "%s, %s, %s"%(qIdx,tIdx,str(simscore))
