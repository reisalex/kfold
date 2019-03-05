import sys
import shelve
import numpy as np
import pandas as pd
import scipy  as sp

sys.path.append('..')
import KFOLDWrapper
from PyVRNA import PyVRNA
ViennaRNA = PyVRNA(parameter_file='rna_andronescu2007.par', dangles=0, pyindex=True)

beta = 0.45
k = 31.5

def get_options(seq, maxtime):
    fold=ViennaRNA.RNAfold(seq)
    return dict(
        fold0='.'*len(seq),
        foldf=fold.structure,
        nsim=1,
        tmax=maxtime,
        trange=KFOLDWrapper.get_trange(maxtime),
        pynsim=1000
        )

def simulate():
    df = pd.read_csv('rbscalc21_ICdatabase.csv')
    used_seqs  = list(df['used_mRNA_sequence'])
    dG_16S_SDs = list(df['dG_SD_16S_hybrid'])
    cutoffs    = list(df['most_5p_paired_SD_mRNA'])
    seqs = [s[c:-21] for s,c in zip(used_seqs,cutoffs)]
    options = [get_options(seq, k*np.exp(beta * dG_16S_SD) ) for seq,dG_16S_SD in zip(seqs,dG_16S_SDs)]

    mean_dGs = []; std_dGs = [];
    for output in KFOLDWrapper.run(seqs,options):
        seq = output['sequence']
        dGs = [ViennaRNA.RNAeval([seq],[fld]) for fld in output['foldsf']]
        mean = np.mean(dGs)
        std  = np.std(dGs)
        mean_dGs.append(mean)
        std_dGs.append(std)

    dfo = pd.DataFrame()
    dfo['dG_mean'] = mean_dGs
    dfo['dG_stdev'] = std_dGs
    dfo.to_csv('all_simulation_output_at_tau_seq2footprint.csv')

if __name__ == "__main__":
    simulate()