
# guess k-value
# scipy.optimize.minimize

# simulate RNA with kfold for:
# time = k*e^beta*dG_complex

# get dG_mRNA, calculate mean dG_total, TIR, & K (proportionality constant)

# dG_complex = dG_SD:aSD + dG_spacing + dG_standby + dG_start
# dG_refold  = dG_mRNA - dG_pre_SD - dG_post_footprint

# dG_final   = dG_mRNA_rRNA + dG_spacing + dG_standby + dG_start
# dG_initial = dG_mRNA

# calculate errors and total error

import sys
import shelve
import numpy as np
import pandas as pd
import scipy  as sp
from scipy.stats import pearsonr
from scipy.optimize import minimize, fmin_bfgs, fmin_l_bfgs_b

sys.path.append('..')
import KFOLDWrapper
from PyVRNA import PyVRNA
ViennaRNA = PyVRNA(parameter_file='rna_andronescu2007.par', pyindex=True)

df = pd.read_csv('JACS_2017.csv')

beta = 0.45

def custom_options(seq, initial_structure, maxtime):
    fold=ViennaRNA.RNAfold(seq)
    return dict(
        fold0=initial_structure,
        foldf=fold.structure,
        ef=0.95*fold.energy,
        nsim=1,
        tmax=maxtime,
        trange=KFOLDWrapper.get_trange(maxtime),
        pynsim=50
        )

calc_y = lambda a0,a1,x: a1*x+a0

def calc_total_error(x,y):
    LSQ = lambda b: np.sum( (y-(beta*x+b))**2.0 )
    res = minimize(LSQ,x0=1,bounds=None)
    residuals = y - calc_y(res.x[0],beta,x)
    return np.sum(residuals**2.0)

def fun(x):
    print "Initializing new iter with k={}".format(x[0])
    seqs      = list(df['used_mRNA_sequence'])
    folds     = df['final_mRNA_structure']
    dG_final  = np.array(df['dG_total']) + np.array(df['dG_mRNA'])
    taus      = x[0]*np.exp(beta*dG_final)
    options   = [custom_options(seq,fold0,tau) for seq,fold0,tau in zip(seqs,folds,taus)]
    dG_mRNAs  = []
    for output in KFOLDWrapper.run(seqs,options):
        seq = output['sequence']
        foldsf = [filter(None,traj)[-1] for traj in output['structures']]
        dGfs   = [ViennaRNA.RNAeval([seq],[fold]) for fold in foldsf]
        dG_mRNAs.append( np.mean(dGfs) )
    dG_totals = dG_final - np.array(dG_mRNAs)
    y = np.log(df['PROT.MEAN'])
    SQE = calc_total_error(dG_totals,y)
    (r,pvalue) = pearsonr(dG_totals,y)
    print "-"*50
    print "k={0:.1f}, Sum square error={1:.2f}, R^2={2:.2f}".format(x[0],SQE,r**2.0)
    print "-"*50
    return r**2.0

def main():
    # minimize(fun, x0=[650.0], bounds=[(500.0,10000.0)])
    # x = [15000.0]
    # rsqlist = []
    # while x[0] < 25000.0:
    #     rsqlist.append((x[0],fun(x)))
    #     x[0] += 1000.0
    #     print rsqlist

    # 50000.0 > R**2.0 = 0.72
    # 75000.0 > R**2.0 = 0.724
    # 100000.0 > R**2.0 = 0.71
    print "R^2={}".format(fun([150000.0]))

if __name__ == "__main__":
    main()
    