
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
        pynsim=10
        )

calc_y = lambda a0,a1,x: a1*x+a0

def calc_total_error(x,y):
    LSQ = lambda b: np.sum( (y-(beta*x+b))**2.0 )
    res = minimize(LSQ,x0=1,bounds=None)
    residuals = y - calc_y(res.x[0],beta,x)
    return np.sum(residuals**2.0), np.exp(res.x[0])

def fun(x, fulltime_mean, fulltime_std):
    print "Initializing new iter with k={}".format(x[0])
    seqs      = list(df['used_mRNA_sequence'])
    folds     = df['final_mRNA_structure']
    dG_final  = np.array(df['dG_total']) + np.array(df['dG_mRNA'])
    taus      = [t if t<10.0 else 10.0 for t in x[0]*np.exp(beta*dG_final)]
    # options   = [custom_options(seq,fold0,tau) for seq,fold0,tau in zip(seqs,folds,taus)]

    dG_mRNAs_avg = []
    dG_mRNAs_std = []

    for i, (seq, fold0, tau) in enumerate(zip(seqs,folds,taus)):

        # extract existing values from fulltime simulations (completed first)
        if tau == 10.0:
            dG_mRNAs_avg.append( fulltime_mean.iloc[i][-1] )
            dG_mRNAs_std.append( fulltime_std.iloc[i][-1]  )
            continue

        # only executing for one sequence
        for output in KFOLDWrapper.run([seq],[custom_options(seq,fold0,tau)]):
            seq = output['sequence']
            foldsf = [filter(None,traj)[-1] for traj in output['structures']]
            dGfs   = [ViennaRNA.RNAeval([seq],[fold]) for fold in foldsf]
            dG_mRNAs_avg.append( np.mean(dGfs) )
            dG_mRNAs_std.append( np.std(dGfs) )

    # calculate totals and statistics
    dG_totals = dG_final - np.array(dG_mRNAs_avg)
    y = np.log(df['PROT.MEAN'])
    SQE, K = calc_total_error(dG_totals,y)
    (r,pvalue) = pearsonr(dG_totals,y)

    return r**2.0, SQE, dG_mRNAs_avg, dG_mRNAs_std, dG_totals, K

def simulate():
    seqs    = list(df['used_mRNA_sequence'])
    folds   = df['final_mRNA_structure']
    options = [custom_options(seq,fold0,10.0) for seq,fold0 in zip(seqs,folds)]
    all_mean_dGs = list()
    all_std_dGs  = list()
    for output in KFOLDWrapper.run(seqs,options):
        seq = output['sequence']
        dGs = [[ViennaRNA.RNAeval([seq],[fold]) for fold in filter(None,traj)] for traj in output['structures']]
        assert all(len(dG_list)==len(dGs[0]) for dG_list in dGs[1:])
        mean_dGs = list()
        std_dGs  = list()
        for i in xrange(len(dGs[0])):
            mean_dGs.append( np.mean( [dG_list[i] for dG_list in dGs] ) )
            std_dGs.append( np.std( [dG_list[i] for dG_list in dGs] ) )

        all_mean_dGs.append(mean_dGs)
        all_std_dGs.append(std_dGs)

    df1 = pd.DataFrame(all_mean_dGs, columns=options[0]['trange'])
    df1.insert(0, column='sequence', value=seqs)
    df1.insert(1, column='initial_structure', value=folds)

    df2 = pd.DataFrame(all_std_dGs, columns=options[0]['trange'])
    df2.insert(0, column='sequence', value=seqs)
    df2.insert(1, column='initial_structure', value=folds)

    writer = pd.ExcelWriter('full_time_kfold_simulations_JACS.xlsx',engine='xlsxwriter')
    df1.to_excel(writer,sheet_name='Sheet1')
    df2.to_excel(writer,sheet_name='Sheet2')
    writer.save()

    return df1, df2

def main(fulltime_mean, fulltime_std):

    kvals = map(float,[100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000, 250000, 500000, 1e6, 2500000])
    kvals = [100.0, 250.0]
    stats_table     = []
    dG_totals_table = []
    dG_mRNAs_avg_table = []
    dG_mRNAs_std_table = []

    for i,k in enumerate(kvals):
        R2, RSS, dG_mRNAs_avg, dG_mRNAs_std, dG_totals, K = fun([k], fulltime_mean, fulltime_std)
        print k, R2, RSS, K
        stats_table.append([k,R2,RSS,K])
        dG_mRNAs_avg_table.append(dG_mRNAs_avg)
        dG_mRNAs_std_table.append(dG_mRNAs_std)
        dG_totals_table.append(dG_totals)

        df1 = pd.DataFrame(stats_table, columns=['k1','R^2','RSS','K'])

        df2 = pd.DataFrame(dG_mRNAs_avg_table).transpose()
        df2.columns = kvals[:i+1]

        df3 = pd.DataFrame(dG_mRNAs_std_table).transpose()
        df3.columns = kvals[:i+1]

        df4 = pd.DataFrame(dG_totals_table).transpose()
        df4.columns = kvals[:i+1]

        writer = pd.ExcelWriter('identify_k.xlsx',engine='xlsxwriter')
        df1.to_excel(writer,sheet_name='Sheet1')
        df2.to_excel(writer,sheet_name='Sheet2')
        df3.to_excel(writer,sheet_name='Sheet3')
        df4.to_excel(writer,sheet_name='Sheet4')
        writer.save()

if __name__ == "__main__":
    fulltime_mean, fulltime_std = simulate()
    main(fulltime_mean, fulltime_std)