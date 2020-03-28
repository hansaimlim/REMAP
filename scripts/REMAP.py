import numpy as np
from scipy import sparse
from numba import jit

import sys
import argparse

def parse_args():
    parser=argparse.ArgumentParser(description="REMAP")
    parser.add_argument('--path',nargs='?', default='/datadir/',help='Input data path.')
    parser.add_argument('--R',nargs='?', default='chem_prot.csv')
    parser.add_argument('--chemsim',nargs='?', default='chem_chem.csv')
    parser.add_argument('--protsim',nargs='?', default='prot_prot.csv')
    parser.add_argument('--low_rank', type=int, default=100,
                        help='Low rank parameter.')
    parser.add_argument('--max_iter', type=int, default=100,
                        help='Maximum iteration.')
    parser.add_argument('--weight', type=float, default=0.1,
                        help='Global weight parameter.')
    parser.add_argument('--imp', type=float, default=0.1,
                        help='Global Imputation for unobserved associations.')
    parser.add_argument('--reg', type=float, default=0.1,
                        help='Regularization parameter.')
    parser.add_argument('--weight_chem', type=float, default=0.75,
                        help='Importance weight for chem-chem.')
    parser.add_argument('--weight_prot', type=float, default=0.25,
                        help='Importance weight for prot-prot.')
    parser.add_argument('--seed', type=int, default=1987,
                        help='Random seed.')
    return parser.parse_args()

def REMAP(R,chem_chem,prot_prot,args):
    
    m,n=R.shape
    Dm=sparse.diags(np.sum(chem_chem,1),0,shape=(m,m))
    Dn=sparse.diags(np.sum(prot_prot,1),0,shape=(n,n))
    Lu=Dm-chem_chem
    Lv=Dn-prot_prot

    U,V=updateUV(R,Lu,Lv,args)

    return (U,V)

def updateUV(R, Lu, Lv, args):
    np.random.seed=args.seed
    lowrank=args.low_rank
    maxite=args.max_iter
    weight=args.weight
    imp=args.imp
    reg=args.reg
    pchem=args.weight_chem
    pprot=args.weight_prot
    m,n=R.shape

    U0=np.asmatrix(np.random.rand(m,lowrank),dtype=np.float64)
    V0=np.asmatrix(np.random.rand(n,lowrank),dtype=np.float64)
    Lu_plus=0.5*(np.abs(Lu)+Lu)
    Lu_minus=0.5*(np.abs(Lu)-Lu)

    Lv_plus=0.5*(np.abs(Lv)+Lv)
    Lv_minus=0.5*(np.abs(Lv)-Lv)

    logging.debug("Preparing indicator matrix...")
    Nzs=np.zeros((m,n))
    row,col=np.nonzero(inputMatrix)
    for i in range(len(row)):
        Nzs[row[i],col[i]]=1
    logging.debug("updateUV started...")

    if R.dtype != np.float64:
        R=R.astype(np.float64)
    if R.__class__ in [sparse.coo_matrix,
                       sparse.csr_matrix,
                       sparse.csc_matrix,
                       sparse.bsr_matrix
                      ]:
        R=R.todense()
    since=time.time()
    countEvery=10
    for ite in range(0,maxite):
        UVT=get_UVT(Nzs,U0,V0)
        U0=fill_na(updateU(R,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,pchem),val=0)
        V0=fill_na(updateU(R.transpose(),UVT.transpose(),weight,imp,Lv_plus,Lv_minus,V0,U0,reg,pprot),val=0)
        if (ite+1)%countEvery==0:
            elapsed=time.time()-since
            logging.debug("{:.4f} minutes elapsed for {} iterations: iter {}".format(elapsed/60.0,countEvery,ite+1))
            since=time.time()
    return (U0,V0)

def get_UVT(Nzs, U, V):
    return jit_mult(Nzs, jit_dot(U,V.T) )

def updateU(R,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,importance):
    m,n=R.shape
    if imp>0:
        ua=(1-weight*imp)*jit_dot(R,V0)+jit_dot( (weight*imp*np.ones((m,n))), V0)+importance*jit_dot(Lu_minus,U0)
    else:
        ua=jit_dot(R,V0)+importance*jit_dot(Lu_minus,U0)
    ub=(1-weight)*jit_dot(UVT,V0)+ (weight* (jit_dot(U0, jit_dot(V0.T,V0) ))) + importance*jit_dot(Lu_plus,U0) + reg*U0
    u=jit_sqrt(jit_divide(ua,ub))
    U1=jit_mult(U0,u)
    return U1

def fill_na(A,val=0):
    A[np.isinf(A)]=val
    A[np.isnan(A)]=val
    return A

@jit(nopython=True)
def jit_dot(A,B):
    return np.dot(A,B)

@jit(nopython=True, parallel=True)
def jit_mult(A,B):
    return np.multiply(A,B)

@jit(nopython=True, parallel=True)
def jit_divide(A,B):
    return np.divide(A,B)

@jit(nopython=True, parallel=True)
def jit_sqrt(A):
    return np.sqrt(A)

@jit(nopython=True, parallel=True)
def jit_rownorm(A,ord=2):
    norm=[]
    for i in range(A.shape[0]):
        norm.append(np.linalg.norm(A[i,:],ord=ord))
    return np.array(norm,dtype=np.float32)

if __name__== '__main__':
    args=parse_args()
    lowrank=args.low_rank
    maxite=args.max_iter
    weight=args.weight
    imp=args.imp
    reg=args.reg

    chemprot=args.path+args.R #chem-prot matrix
    chemchem=args.path+args.chemsim 
    protprot=args.path+args.protsim
    U,V=REMAP(chemprot,chemchem,protprot,args) #run remap
    P=jit_dot(U,V.T) #prediction matrix
    np.savetxt('REMAP_pred.csv',P,delimiter=',') #save prediction matrix
