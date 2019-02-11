
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
ViennaRNA = PyVRNA(parameter_file='rna_andronescu2007.par', dangles=0, pyindex=True)

df = pd.read_csv('JACS_2017.csv')
beta = 0.45

def custom_options(seq, initial_structure, kvals, dG_final):
    fold=ViennaRNA.RNAfold(seq)
    times = [t for t in kvals/np.exp(-beta*dG_final) if t < 1000.0]
    times.append(1000.0)
    tmax = times[-1]
    times += [-1.0]*(100-len(times[:]))
    assert len(times) == 100
    return dict(
        fold0=initial_structure,
        foldf=fold.structure,
        ef=0.95*fold.energy,
        nsim=1,
        tmax=tmax,
        trange=np.array(times),#np.array(times),
        pynsim=1000
        )

def get_trange(tmax):
    dt = 0.01
    io = 1
    trange = [dt]
    while (dt < tmax):
        trange.append(trange[-1]+dt)
        io += 1
        if (io > 9):
            io = 1
            dt *= 10.0
    return trange

def custom_options2(seq, initial_structure):
    fold=ViennaRNA.RNAfold(seq)
    maxtime=1000.0
    return dict(
        fold0=initial_structure,
        foldf=fold.structure,
        ef=0.95*fold.energy,
        nsim=1,
        tmax=maxtime,
        trange=get_trange(maxtime),
        pynsim=1000
        )

def get_default_options(seq):
    maxtime=1000.0
    fold=ViennaRNA.RNAfold(seq)
    return dict(
        fold0='.'*len(seq),               # starting structure
        foldf=fold.structure,             # first passage time search structure
        ef=0.95*fold.energy,              # energy first passage time threshold value
        nsim=1,                           # number of simulations internal to kfold (serial)
        tmax=maxtime,                     # total simulation time (kfold units)
        trange=get_trange(maxtime),       # time points at which to collect data
        pynsim=1000                       # number of simulations through multiprocessing (parallel)
        )                                 # total number of simulations ends up being nsim*pynsim


calc_y = lambda a0,a1,x: a1*x+a0

def calc_total_error(x,y):
    LSQ = lambda b: np.sum( (y-(beta*x+b))**2.0 )
    res = minimize(LSQ,x0=1,bounds=None)
    residuals = y - calc_y(res.x[0],beta,x)
    return np.sum(residuals**2.0), np.exp(res.x[0])

def main(kvals):
    kvals     = np.array(kvals)
    seqs      = list(df['used_mRNA_sequence'])
    folds     = df['final_mRNA_structure']
    dG_finals = np.array(df['dG_total']) + np.array(df['dG_mRNA'])

    dG_mRNAs_avg = []; dG_mRNAs_std = []; all_times = [];

    options = [custom_options(seq,fold0,kvals,dG_final) for seq, fold0, dG_final in zip(seqs,folds,dG_finals)]

    # run simulations and iterate over sequences to collect kfold data
    for output in KFOLDWrapper.run(seqs,options):
        seq   = output['sequence']
        times = [val for val in output['options']['trange'] if val != -1.0]
        times += [times[-1]]*(len(kvals)+1-len(times))
        all_times.append(times)
        avg = []; std = [];

        # get time slices, and mean/std values at each time
        for i,time in enumerate(output['options']['trange']):
            if time == -1.0: break
            folds = [traj[i] for traj in output['structures']] # ensemble at given time
            dGs   = [ViennaRNA.RNAeval([seq],[fold]) for fold in folds] # recalculated dGs
            avg.append( np.mean(dGs) )
            std.append( np.std(dGs)  )

        # extend dG_mRNAs_avg and dG_mRNAs_std by end point to "backfill" kvals where the
        # simulation tau exceeded the max time (1000.0)
        avg   += [avg[-1]]*(len(kvals)+1-len(avg))
        std   += [std[-1]]*(len(kvals)+1-len(std))

        # append to master list
        dG_mRNAs_avg.append(avg)
        dG_mRNAs_std.append(std)

    assert all(len(a)==len(kvals)+1 for a in dG_mRNAs_avg)
    assert all(len(a)==len(kvals)+1 for a in dG_mRNAs_std)
    assert all(len(a)==len(kvals)+1 for a in all_times)

    # dG_mRNAs_avg = [[seq1...],[seq2...],[seq3...]]
    # calculate totals and statistics
    stats_table   = []
    all_dG_totals = []
    for i,k in enumerate(kvals):
        dG_totals = dG_finals - np.array([avg[i] for avg in dG_mRNAs_avg])
        y = np.log(df['PROT.MEAN'])
        RSS, K = calc_total_error(dG_totals,y)
        (r,pvalue) = pearsonr(dG_totals,y)
        stats_table.append([k,r**2.0,RSS,K])
        all_dG_totals.append(dG_totals)

    # do data-point beyond last k (equilibrium or tmax value)
    dG_totals = dG_finals - np.array([avg[i+1] for avg in dG_mRNAs_avg])
    y = np.log(df['PROT.MEAN'])
    RSS, K = calc_total_error(dG_totals,y)
    (r,pvalue) = pearsonr(dG_totals,y)
    stats_table.append(['inf',r**2.0,RSS,K])
    all_dG_totals.append(dG_totals)


    df1 = pd.DataFrame(stats_table, columns=['k1','R^2','RSS','K'])

    headers = list(kvals) + ['inf']

    df2 = pd.DataFrame(dG_mRNAs_avg)#.transpose()
    df2.columns = headers

    df3 = pd.DataFrame(dG_mRNAs_std)#.transpose()
    df3.columns = headers

    df4 = pd.DataFrame(all_times)
    df4.columns = headers

    df5 = pd.DataFrame(all_dG_totals).transpose()
    df5.columns = headers

    writer = pd.ExcelWriter('identify_k.xlsx',engine='xlsxwriter')
    df1.to_excel(writer,sheet_name='stats')
    df2.to_excel(writer,sheet_name='dG_mRNA_mean')
    df3.to_excel(writer,sheet_name='dG_mRNA_std')
    df4.to_excel(writer,sheet_name='taus')
    df5.to_excel(writer,sheet_name='dG_total')
    writer.save()


    return stats_table, dG_mRNAs_avg, dG_mRNAs_std, all_dG_totals

def simulate():
    seqs    = list(df['used_mRNA_sequence'])
    folds   = df['final_mRNA_structure']
    options = [custom_options(seq,fold0) for seq,fold0 in zip(seqs,folds)]
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


if __name__ == "__main__":
    simulate()
    # kvals = map(float,[100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000, 250000, 500000, 1e6, 2500000, 5000000, 7500000, 1e7])
    # main(kvals)