#!/usr/bin/python

dti='../ZINC_DTI.tsv' #drug-target interaction in text identifiers
chems='../list/ZINC_chemicals.tsv' #list of unique chemicals with chemical index
prots='../list/ZINC_protein_index.tsv' #list of unique proteins with protein index

idx2zinc={}
zinc2idx={}
idx2prot={}
prot2idx={}

with open(chems,"r") as chemline:
    next(chemline) #skip header line
    for line in chemline:
        line=line.strip().split("\t") #tab-separated
        idx=str(line[0])
        zinc=str(line[1])
        idx2zinc[idx]=zinc
        zinc2idx[zinc]=idx

with open(prots,"r") as protline:
    next(protline)
    for line in protline:
        line=line.strip().split("\t")
        idx=str(line[0])
        prot=str(line[1])
        idx2prot[idx]=prot
        prot2idx[prot]=idx

with open(dti,"r") as dtiline:
    next(dtiline) #skip header line
    for line in dtiline:
        line=line.strip().split("\t")
        zinc=str(line[0]) #ZINC ID
        prot=str(line[1]) #Gene ID
        chemIdx=zinc2idx[zinc] #chemical index
        protIdx=prot2idx[prot] #protein index

        print "%s,%s"%(chemIdx,protIdx)
