#!/usr/bin/python

Idx='../results/REMAP_prediction.txt' #Predicted drug-target association by indexes
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

with open(Idx,"r") as idxline:
    for line in idxline:
        line=line.strip().split("\t")
        try:
            chemIdx=str(line[0]) #ZINC ID
            protIdx=str(line[1]) #Gene ID
            score=str(line[2])
            chem=idx2zinc[chemIdx] #ZINC chemical ID
            prot=idx2prot[protIdx] #Protein ID

            print "%s, %s, %s"%(chem,prot,score)
        except:
            None
