import numpy as np
from scipy import sparse
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

    lowrank=args.low_rank
    maxite=args.max_iter
    weight=args.weight
    imp=args.imp
    reg=args.reg
    pchem=args.weight_chem
    pprot=args.weight_prot
    m,n=R.shape

    U0=np.random.rand(m,lowrank)
    V0=np.random.rand(n,lowrank)
    Lu_plus=0.5*(np.abs(Lu)+Lu)
    Lu_minus=0.5*(np.abs(Lu)-Lu)

    Lv_plus=0.5*(np.abs(Lv)+Lv)
    Lv_minus=0.5*(np.abs(Lv)-Lv)

    for ite in range(0,maxite):
        UVT=get_UVT(R,U0,V0)
        U0=updateU(R,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,pchem)
        V0=updateU(R.transpose(),UVT.transpose(),weight,imp,Lv_plus,Lv_minus,V0,U0,reg,pprot)

    return (U0,V0)

def get_UVT(R, U, V):
    m,n=R.shape
    
    rows,cols=R.nonzero()
    UVT=np.zeros((m,n),dtype=float)
    for i in range(len(rows)):
        row=rows[i]
        col=cols[i]
        UVT[row,col]=np.dot(U[row],V[col].T)
    return UVT

def updateU(R,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,importance):

    m,n=R.shape

    ua=(1-weight*imp)*np.dot(R,V0)+np.dot( (weight*imp*np.ones((m,n))), V0)+importance*np.dot(Lu_minus,U0)
    ub=(1-weight)*np.dot(UVT,V0)+ (weight* (np.dot(U0, np.dot(V0.T,V0) ))) + importance*np.dot(Lu_plus,U0) + reg*U0
    u=np.sqrt(np.divide(ua,ub))
    U1=np.multiply(U0,u)
    U1[np.isinf(U0)]=0
    U1[np.isnan(U0)]=0
    return U1

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
    P=np.dot(U,V.T) #prediction matrix
    np.savetxt('REMAP_pred.csv',P,delimiter=',') #save prediction matrix
