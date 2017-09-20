#!/usr/bin/python

import os
import sys
import math
import subprocess
from multiprocessing import Pool
from multiprocessing import freeze_support, cpu_count
global screenmd_path_global
screenmd_path_global='/home/hlim/ChemAxon/JChem/bin/'  #Change this to your JChem bin directory (where screenmd binary file is found)

def get_tani_sim_large(large_smi_file, lines, threshold, outfile):
    #split smi file into pieces of equal lines
    # bash command split large_smi_file -l lines
    large_smi_file=str(large_smi_file)
    lines=int(lines)
    threshold=float(threshold)
    outfile=str(outfile)
    split_smi_files=split_file(large_smi_file,lines)
    filecount=len(split_smi_files)
    inputs=[]
    for finfo in split_smi_files:
        inp=(finfo[0],large_smi_file,finfo[1],outfile,threshold,screenmd_path_global)
        inputs.append(inp)
    with Pool(cpu_count()-1) as pool: #use multicore
        try:
            R=pool.starmap(smile_to_sim,inputs)
            print("Process Done. Tanimoto similarity file= %s"%outfile)
            print("Similarity threshold: %s (minimum similarity)."%str(threshold))
        except BaseException as e:
            print("Error occurred:\n"+str(e)+"\n")
            for finfo in split_smi_files:
                delete_file(finfo[0])
            print("Splitted smi files deleted.\n")
    
def smile_to_sim(query,target,index_padding,outfile,threshold,screenmd_path):
    print("Screening %s to %s..."%(str(query),str(target)))
    tani=run_screenmd(screenmd_path,query,target)
    calc_sim(tani,int(index_padding),float(threshold),outfile)
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
    
def run_screenmd(screenmd_path,query_file,target_file):
    tanimoto_filename=query_file+".dat"
    try:
        screenmd_command=os.path.join(screenmd_path,'screenmd')+" %s %s -k ECFP -g -c -M Tanimoto > %s"%(query_file,target_file,tanimoto_filename)
        output=subprocess.check_output(['bash','-c',screenmd_command])
        return tanimoto_filename
    except BaseException as e:
        print('Error occurred during screenmd: '+str(e)+"\n")
        print("Suggestion: Where is your JChem software installed?")
        print("Suggestion: Maybe you need to change screenmd_path_global variable at the bottom of the code.\n")

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
        print("(e.g.) python3 get_tani_sim_multicore.py chemicals.smi 5000 0.5 tanimoto_sim.csv")
        sys.exit()
    freeze_support()
    get_tani_sim_large(Args[0],int(Args[1]),float(Args[2]),Args[3])
