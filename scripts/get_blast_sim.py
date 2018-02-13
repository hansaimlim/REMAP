#!/usr/bin/python

import os
import sys
import math
import subprocess
global blast_path
blast_path='/home/hlim/blast/ncbi-blast-2.7.0+/bin/'

def get_blast_sim(protinfo_file, fasta_file, threshold, outfile, num_cores=1):
    fasta_file=str(fasta_file)
    threshold=float(threshold)
    outfile=str(outfile)

    print("Process starts. protinfo file: %s"%protinfo_file)
    print("Similarity threshold: %s (minimum similarity)."%str(threshold))
    blastdbname=fasta_file+"_blastdb"
    blast_result_file="blastresult.dat"
    blastdbcommand=blast_path+"makeblastdb -in %s -dbtype prot -out %s"%(fasta_file,blastdbname)
    output=subprocess.check_output(['bash','-c',blastdbcommand])
    blastpcommand=blast_path+"blastp -query %s -db %s -evalue 1e-5 -outfmt 6 -num_threads %s > %s"%(fasta_file, blastdbname, str(num_cores), blast_result_file)
    output=subprocess.check_output(['bash','-c',blastpcommand])
    calc_sim(protinfo_file,blast_result_file, threshold, outfile)

    delete_file(blastdbname+".phr")
    delete_file(blastdbname+".pin")
    delete_file(blastdbname+".psq")
    delete_file(blast_result_file)
    print("Process Done. BLAST similarity file= %s"%outfile)
    print("Similarity threshold: %s (minimum similarity)."%str(threshold))
    
def delete_file(filename):
    delcommand="rm %s"%str(filename)
    output=subprocess.check_output(['bash','-c',delcommand])

def get_protinfo(filename):
    idx2acc={}
    acc2idx={}
    with open(filename) as inf:
        next(inf)
        for line in inf:
            line=line.strip().split("\t")
            idx=str(line[0]).strip()
            tid=str(line[1])
            acc=str(line[2]).strip()
            idx2acc[idx]=acc
            acc2idx[acc]=idx
    return idx2acc, acc2idx

def calc_sim(protinfo_file, datfile, threshold, outfile):
    #calculate tanimoto similarity
    #query index start from 1+query_idx_pad (for splitted file)
    #score lower than threshold discarded
    threshold=float(threshold)
    outfile=str(outfile)
    fout=open(outfile,"a+")
    idx2selfscore={}
    try:
        idx2acc, acc2idx=get_protinfo(protinfo_file)    
    except BaseException as e:
        print("Error occurred:\n"+str(e))
        print("Please provide valid protein info file.")
    with open(datfile,"r") as inf:
        for line in inf:
            line=line.strip().split('\t')
            qids=line[0].strip().split("|")
            tids=line[1].strip().split("|")
            score=float(line[-1])
            qaccession=str(qids[1])
            taccession=str(tids[1])
            qidx=acc2idx[qaccession]
            tidx=acc2idx[taccession]
            if qidx==tidx:
                idx2selfscore[qidx]=score
    with open(datfile,"r") as inf:
        for line in inf:
            line=line.strip().split('\t')
            qids=line[0].strip().split("|")
            tids=line[1].strip().split("|")
            score=float(line[-1])
            qaccession=str(qids[1])
            taccession=str(tids[1])
            qidx=acc2idx[qaccession]
            tidx=acc2idx[taccession]

            simscore=float(score / idx2selfscore[qidx])
            if simscore>=threshold:
                fout.write("%s, %s, %s\n"%(str(qidx),str(tidx),str(simscore)))
    fout.close()
if __name__ == '__main__':
    Args=sys.argv[1:]
    if len(Args)<4:
        print("Insufficient arguments.\nProtein_info_file, FASTA file, similarity threshold, outfile are required arguments.")
        print("(e.g. 1-core) python get_blast_sim.py protInfo.tsv proteins.fas 0.5 blast_sim.csv  (OR)")
        print("(e.g. 4-core) python get_blast_sim.py protInfo.tsv proteins.fas 0.5 blast_sim.csv 4")
        sys.exit()
    if len(Args)==5:
        get_blast_sim(Args[0],Args[1],float(Args[2]),Args[3],Args[4])
    elif len(Args)==4:
        get_blast_sim(Args[0],Args[1],float(Args[2]),Args[3])
