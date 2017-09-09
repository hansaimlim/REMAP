#!/usr/bin/python

#Use this script if you have large number of chemicals (e.g. >20000)
#It splits large smi file into pieces and get Tanimoto similarity

import os
import sys
import math
import subprocess

def get_tani_sim_large(large_smi_file, lines, threshold, outfile, screenmd_path):
    #split smi file into pieces of equal lines
    # bash command split large_smi_file -l lines
    large_smi_file=str(large_smi_file)
    lines=int(lines)
    threshold=float(threshold)
    outfile=str(outfile)
    split_smi_files=split_file(large_smi_file,lines)
    tanimoto_filenames=run_screenmd(screenmd_path,split_smi_files,large_smi_file)
    print tanimoto_filenames
    calc_sim(tanimoto_filenames, lines, threshold, outfile)
    delete_files(split_smi_files) #delete splitted files
    print "Splitted files deleted."
    delete_files(tanimoto_filenames)
    print "Screenmd files deleted."
    print "Process Done. Tanimoto similarity file= %s"%outfile
    print "Similarity threshold: %s (minimum similarity)."%str(threshold)
    
def split_file(infile, lines):
    files=[]
    wccommand="wc -l %s"%str(infile)
    output=subprocess.check_output(['bash','-c',wccommand])
    total_lines=output.strip().split(' ')[0]
    filecount=int(math.ceil(float(total_lines)/float(lines)))
    splitcommand="split %s -ed -a 3 -l %d" %(str(infile),int(lines)) #file name will be x000, x001, x002...
    output = subprocess.check_output(['bash','-c',splitcommand])
    print "File splitted into %s pieces."%str(filecount)
    suffix=0
    for i in range(0,filecount):
        str_suffix=str(format(suffix,'03d'))
        fname="x"+str_suffix
        files.append(fname)
        suffix+=1
    return files
def delete_files(filenames):
    for fname in filenames:
        delcommand="rm %s"%str(fname)
        output=subprocess.check_output(['bash','-c',delcommand])

def run_screenmd(screenmd_path,query_files,target_file):
    tanimoto_filenames=[]
    for qfile in query_files:
        tanifile=qfile+".dat"
        screenmd_command=os.path.join(screenmd_path,'screenmd')+" %s %s -k ECFP -g -c -M Tanimoto > %s"%(qfile,target_file,tanifile)
        output=subprocess.check_output(['bash','-c',screenmd_command])
        tanimoto_filenames.append(tanifile)
    return tanimoto_filenames

def calc_sim(datfiles, lines, threshold, outfile):
    #calculate tanimoto similarity
    #query index start from 1+query_idx_pad (for splitted file)
    #score lower than threshold discarded
    query_idx_pad=0
    threshold=float(threshold)
    outfile=str(outfile)
    fout=open(outfile,"w")
    for datfile in datfiles:
        with open(datfile,"r") as inf:
            next(inf)
            for line in inf:
                line=line.strip().split('\t')
                qidx=query_idx_pad+int(line[0])
                for tidx in range(1,len(line)):
                    score=float(1.0)-float(line[tidx])
                    if score>=threshold:
                        fout.write("%s, %s, %s\n"%(str(qidx),str(tidx),str(score)))
        query_idx_pad+=int(lines)
    fout.close()
if __name__ == '__main__':
    screenmd_path='/home/hlim/ChemAxon/JChem/bin/' #directory where screenmd is installed
    
    Args=sys.argv[1:]
    if len(Args)<4:
        print "Insufficient arguments.\nSMILES file, split line num, similarity threshold, outfile are required arguments."
        print "(e.g.) python get_tani_sim.py chemicals.smi 5000 0.5 tanimoto_sim.csv"
        sys.exit()
        
    get_tani_sim_large(Args[0],int(Args[1]),float(Args[2]),Args[3])
