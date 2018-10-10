#!/usr/bin/python

import os
import sys
import math
import pickle
import random
import string
import subprocess
from multiprocessing import Pool
from multiprocessing import freeze_support, cpu_count
import psutil
def load_pickle(filename):
  with open(filename,'rb') as handle:
    data=pickle.load(handle)
  return data
def save_pickle(data,filename):
  with open(filename,'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

global madfast_path_global
madfast_path_global='/home/hansaimlim/ChemAxon/madfast-cli-0.2.3/bin/'  #Change this to your madfast bin directory

def get_idconv():
  idx2ikey=load_pickle('../idx2ikey.pickle')
  ikey2idx=load_pickle('../ikey2idx.pickle')
  idx2smi=load_pickle('../idx2smi.pickle')
  return (idx2ikey,ikey2idx,idx2smi)

def get_tani_sim_wholematrix(whole_smiles, threshold, outfile):
    #split smi file into pieces of equal lines
    # bash command split large_smi_file -l lines
    lines=2000
    npc=5 #num processes to run at a time
    threshold=float(threshold)
    outfile=str(outfile)
    split_smi_files=split_file(whole_smiles,lines)
    filecount=len(split_smi_files)
   
    inputs=[]
    for qfile in split_smi_files: #query-to-large file
      inp=(qfile,whole_smiles,threshold,outfile)
      inputs.append(inp)
    with Pool(npc) as pool: #use multicore
      try:
        print("{} total processes. Run {} processes at a time...".format(len(inputs)+1,npc))
        R=pool.starmap(get_madfast_result,inputs)
      except BaseException as e:
        print("Error occurred:\n"+str(e)+"\n")
    return split_smi_files


def get_madfast_result(qfile,large_smi_file,threshold,outfile):
  idx2ikey,ikey2idx,idx2smi=get_idconv()
  print("MadFast search for {} started...".format(qfile))
  OF=run_madfast(qfile,large_smi_file,threshold)
  print("Dissimilarity matrix {} created".format(OF))
#  calc_sim(OF, ikey2idx, threshold, outfile)
  parse_sim(OF, ikey2idx, outfile)
  print("Similarity information appended to {}".format(outfile))
  delete_file(OF)
  print("Dissimilarity matrix {} deleted".format(OF))
 
def split_file(infile, lines):
    #split file by given lines and return file names
    files=[] #contains filenames
    wccommand="wc -l %s"%str(infile)
    output=subprocess.check_output(['bash','-c',wccommand])
    output=output.decode('utf8')
    total_lines=output.strip().split(' ')[0]
    print(total_lines)
    filecount=int(math.ceil(float(total_lines)/float(lines)))
    splitcommand="split %s -d -a 4 -l %d" %(str(infile),int(lines)) #file name will be x0000, x0001, x0002...
    output = subprocess.check_output(['bash','-c',splitcommand])
    print("File splitted into %s pieces."%str(filecount))
    suffix=int(0)
    for i in range(0,filecount):
        str_suffix=str(format(suffix,'04d'))
        fname="x"+str_suffix
        files.append(fname)
        suffix+=1
    return files

def delete_file(filename):
    delcommand="rm %s"%str(filename)
    output=subprocess.check_output(['bash','-c',delcommand])
    
def run_madfast(query_file,target_file,cutoff):
    cutoff=float(cutoff)
    #random named dissimilarity matrix to prevent collision
    tanimoto_filename=query_file+''.join(random.choice(string.ascii_uppercase + string.digits+ string.ascii_lowercase) for _ in range(20))
    try:
        madfast_command=os.path.join(madfast_path_global,'searchStorage.sh')+" -context createSimpleEcfp4Context -qmf {} -qidname\
 -tmf {} -tidname -mode MOSTSIMILARS -maxdissim {} -out {}".format(query_file,target_file,(1.0-cutoff),tanimoto_filename)
        output=subprocess.check_output(['bash','-c',madfast_command])
        return tanimoto_filename
    except BaseException as e:
        print('Error occurred during searchStorage: '+str(e)+"\n")
        with open('./madfast_cutoff_error_command.log','a') as errorlog:
          errorlog.write("{}\n".format(madfast_command))

def parse_sim(datfile, ikey2idx, outfile):
    #parse dissimilarity list
    #assume all listed pairs are filtered by cutoff, no cutoff needed here
    outfile=str(outfile)
    fout=open(outfile,"a+")
    linenum=0
    with open(datfile,"r") as inf:
      next(inf)
      for line in inf:
        if line=='': #skip empty lines
          continue
        line=line.strip().split('\t')
        ikey1=str(line[0])
        ikey2=str(line[1])
        dissim=float(line[2])
        sim=float(1.0-dissim)
        try:
          idx1=ikey2idx[ikey1]
          idx2=ikey2idx[ikey2]
        except:
          continue
        fout.write("{0:},{1:},{2:.10f}\n".format(idx1,idx2,sim))
    fout.close()
def calc_sim(datfile, ikey2idx, threshold, outfile):
    #calculate tanimoto similarity
    #query index start from 1+query_idx_pad (for splitted file)
    #score lower than threshold discarded
    threshold=float(threshold)
    outfile=str(outfile)
    fout=open(outfile,"a+")
    linenum=0
    with open(datfile,"r") as inf:
        header=inf.readline()
        header=header.strip().split('\t')
        queries=header[1:]
        colidx=0
        colidx2qidx={}
        for query in queries:
          qikey=query.strip().split(' ')
          qikey=qikey[1].strip()
          qidx=ikey2idx[qikey]
          colidx2qidx[colidx]=qidx
          colidx+=1
        next(inf) #contains empty line
        for line in inf:
            linenum+=1
            if line=='': #skip empty lines
              continue
            line=line.strip().split('\t')
            targetikey=str(line[0])
            if targetikey=='':
              continue
            try:
              tidx=ikey2idx[targetikey]
            except:
              print("Error occurred on line {}, file={}, InChIKey={}".format(linenum,datfile,targetikey))
              continue
            scores=line[1:]
            for col in range(len(scores)):
                score=float(1.0)-float(scores[col])
                qidx=colidx2qidx[col]
                if score>=threshold:
                    fout.write("%s,%s,%8.5f\n"%(str(qidx),str(tidx),score))
    fout.close()
    
if __name__ == '__main__':
  psutil.cpu_count()
  p = psutil.Process()
  p.cpu_affinity()  # get
  p.cpu_affinity([14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39])  # set; from now on, process will run on CPU #15~ only
  threshold=0.3
  outfile='./chembl24plus_chem_sim_threshold03.csv'
#  outfile='./testout_threshold03.csv'
  whole_smiles='./chembl-24_plus.smi'#assume you already have all-to-all similarity for this large file
#  whole_smiles='./testsmi.smi'#assume you already have all-to-all similarity for this large file
  split_smi_files=get_tani_sim_wholematrix(whole_smiles, threshold, outfile)
  for qfile in split_smi_files:
    delete_file(qfile)
    print("{} deleted".format(qfile))
