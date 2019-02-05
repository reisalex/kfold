
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
        pynsim=1000
        )

calc_y = lambda a0,a1,x: a1*x+a0

def calc_total_error(x,y):
    LSQ = lambda b: np.sum( (y-(beta*x+b))**2.0 )
    res = minimize(LSQ,x0=1,bounds=None)
    residuals = y - calc_y(res.x[0],beta,x)
    return np.sum(residuals**2.0), np.exp(res.x[0])

def fun(x):
    print "Initializing new iter with k={}".format(x[0])
    seqs      = list(df['used_mRNA_sequence'])
    folds     = df['final_mRNA_structure']
    dG_final  = np.array(df['dG_total']) + np.array(df['dG_mRNA'])
    taus      = x[0]*np.exp(beta*dG_final)
    taus      = [t if t<1000.0 else 1000.0 for t in taus[:]]
    options   = [custom_options(seq,fold0,tau) for seq,fold0,tau in zip(seqs,folds,taus)]
    dG_mRNAs  = []
    for output in KFOLDWrapper.run(seqs,options):
        seq = output['sequence']
        foldsf = [filter(None,traj)[-1] for traj in output['structures']]
        dGfs   = [ViennaRNA.RNAeval([seq],[fold]) for fold in foldsf]
        dG_mRNAs.append( np.mean(dGfs) )
    dG_totals = dG_final - np.array(dG_mRNAs)
    y = np.log(df['PROT.MEAN'])
    SQE, K = calc_total_error(dG_totals,y)
    (r,pvalue) = pearsonr(dG_totals,y)
    return r**2.0, SQE, dG_mRNAs, dG_totals, K

def simulate():
    seqs    = list(df['used_mRNA_sequence'])
    folds   = df['final_mRNA_structure']
    options = [custom_options(seq,fold0,1000.0) for seq,fold0 in zip(seqs,folds)]
    all_mean_dGs = list()
    for output in KFOLDWrapper.run(seqs,options):
        seq = output['sequence']
        dGs = [[ViennaRNA.RNAeval([seq],[fold]) for fold in filter(None,traj)] for traj in output['structures']]
        assert all(len(dG_list)==len(dGs[0]) for dG_list in dGs[1:])
        mean_dGs = []
        for i in xrange(len(dGs[0])):
            mean_dGs.append( np.mean( [dG_list[i] for dG_list in dGs] ) )
        all_mean_dGs.append(mean_dGs)

    df2 = pd.DataFrame(all_mean_dGs, columns=options[0]['trange'])
    df2.insert(0, column='sequence', value=seqs)
    df2.insert(1, column='initial_structure', value=folds)
    df2.to_csv('full_time_kfold_simulations_JACS.csv')


def main():

    kvals = map(float,[
    100,
    250,
    500,
    750,
    1000,
    2500,
    5000,
    7500,
    10000,
    25000,
    50000,
    75000,
    100000,
    250000,
    500000,
    1e6,
    2500000
    ])

    stats_table     = []
    dG_mRNAs_table  = []
    dG_totals_table = []
    for k in kvals:
        R2, RSS, dG_mRNAs, dG_totals, K = fun([k])
        stats_table.append([k,R2,RSS,K])
        dG_mRNAs_table.append(dG_mRNAs)
        dG_totals_table.append(dG_totals)

    df1 = pd.DataFrame(stats_table, columns=['k1','R^2','RSS','K'])

    df2 = pd.DataFrame(dG_mRNAs_table).transpose()
    df2.columns = kvals

    df3 = pd.DataFrame(dG_totals_table).transpose()
    df3.columns = kvals

    writer = pd.ExcelWriter('identify_k.xlsx',engine='xlsxwriter')
    df1.to_excel(writer,sheet_name='Sheet1')
    df2.to_excel(writer,sheet_name='Sheet2')
    df3.to_excel(writer,sheet_name='Sheet3')
    writer.save()

if __name__ == "__main__":
    main()
    simulate()