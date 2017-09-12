#!/usr/bin/python

import os
import sys
import math
import subprocess
from multiprocessing import Pool
from multiprocessing import freeze_support

def get_tani_sim_large(large_smi_file, lines, threshold, outfile):
    #split smi file into pieces of equal lines
    # bash command split large_smi_file -l lines
    large_smi_file=str(large_smi_file)
    lines=int(lines)
    threshold=float(threshold)
    outfile=str(outfile)
    split_smi_files=split_file(large_smi_file,lines)
    filecount=len(split_smi_files)
#    query_index_padding=int(0)
#    file_done=0
    inputs=[]
    for finfo in split_smi_files:
        inp=(finfo[0],large_smi_file,finfo[1],outfile)
        inputs.append(inp)
    with Pool() as pool: #use multicore
        R=pool.starmap(smile_to_sim,inputs)
#        tani=run_screenmd('/home/hlim/ChemAxon/JChem/bin/',smi,large_smi_file)
#        calc_sim(tani, query_index_padding, threshold, outfile)
#        delete_file(smi)
#        delete_file(tani)
#        query_index_padding+=int(lines)
#        file_done+=1
#        print("%s out of %s files complete..." % (str(file_done),str(filecount)))
    print("Process Done. Tanimoto similarity file= %s"%outfile)
    print("Similarity threshold: %s (minimum similarity)."%str(threshold))
def smile_to_sim(query,target,index_padding,outfile,threshold=float(0.5)):
    print("Screening %s to %s..."%(str(query),str(target)))
    tani=run_screenmd(query,target)
    print("Screening %s to %s is Done. Calculating similarity... Threshold %s."%(str(query),str(target),str(threshold)))
    calc_sim(tani,int(index_padding),float(threshold),outfile)
    print("%s to %s Complete. output to %s."%(str(query),str(target),str(outfile)))
    delete_file(query)
    delete_file(tani)

def split_file(infile, lines):
    #returns file names with line padding
    files=[] #contains tuples (file_name, linepadding)
    wccommand="wc -l %s"%str(infile)
    output=subprocess.check_output(['bash','-c',wccommand])
    output=output.decode('utf8')
    print(output)
    total_lines=output.strip().split(' ')[0]
    print(total_lines)
    filecount=int(math.ceil(float(total_lines)/float(lines)))
    splitcommand="split %s -ed -a 3 -l %d" %(str(infile),int(lines)) #file name will be x000, x001, x002...
    output = subprocess.check_output(['bash','-c',splitcommand])
    print("File splitted into %s pieces."%str(filecount))
    suffix=int(0)
    linepadding=int(0)
    for i in range(0,filecount):
        str_suffix=str(format(suffix,'03d'))
        fname="x"+str_suffix
        tup=(fname,linepadding)
        files.append(tup)
        suffix+=1
        linepadding+=int(lines)
    return files
def delete_file(filename):
    delcommand="rm %s"%str(filename)
    output=subprocess.check_output(['bash','-c',delcommand])
def run_screenmd(query_file,target_file):
    screenmd_path='/home/hlim/ChemAxon/JChem/bin/'    #path to screenmd binary file
    tanimoto_filename=query_file+".dat"
    screenmd_command=os.path.join(screenmd_path,'screenmd')+" %s %s -k ECFP -g -c -M Tanimoto > %s"%(query_file,target_file,tanimoto_filename)
    output=subprocess.check_output(['bash','-c',screenmd_command])
    return tanimoto_filename

def calc_sim(datfile, query_idx_pad, threshold, outfile):
    #calculate tanimoto similarity
    #query index start from 1+query_idx_pad (for splitted file)
    #score lower than threshold discarded
    threshold=float(threshold)
    outfile=str(outfile)
    fout=open(outfile,"a+")
    with open(datfile,"r") as inf:
        next(inf)
        for line in inf:
            line=line.strip().split('\t')
            qidx=int(query_idx_pad)+int(line[0])
            for tidx in range(1,len(line)):
                score=float(1.0)-float(line[tidx])
                if score>=threshold:
                    fout.write("%s, %s, %8.5f\n"%(str(qidx),str(tidx),score))
    fout.close()
if __name__ == '__main__':
    Args=sys.argv[1:]
    if len(Args)<4:
        print("Insufficient arguments.\nSMILES file, split line num, similarity threshold, outfile are required arguments.")
        print("(e.g.) python3 get_tani_sim.py chemicals.smi 5000 0.5 tanimoto_sim.csv")
        sys.exit()
    freeze_support()
    get_tani_sim_large(Args[0],int(Args[1]),float(Args[2]),Args[3])

