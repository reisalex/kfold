
import time
import multiprocessing
from contextlib import contextmanager
from itertools import repeat

import kfold
import tqdm

# load PyVRNA interface to ViennaRNA
from PyVRNA import PyVRNA
ViennaRNA = PyVRNA(parameter_file='rna_turner1999.par',pyindex=True)

# default options
NCPUS = multiprocessing.cpu_count()

alphabet=set('ATCGUatcgu')

# get log-spaced time steps from 0.01 to tmax
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

# return a dictionary containing the default options for simulation
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

# check options dictionary passed by a user (improve later)
def check_options(seq,optsdict):
    assert len(optsdict['fold0'])==len(seq)
    assert len(optsdict['foldf'])==len(seq)
    assert isinstance(optsdict['ef'],float)
    assert isinstance(optsdict['nsim'],int)
    assert isinstance(optsdict['tmax'],float)
    # assert isinstance(optsdict['trange'],list)
    assert isinstance(optsdict['pynsim'],int)

# generate tuple for pool.map from dictionary
def make_input_tuple(seq,optsdict):
    return tuple([seq]+[optsdict[x] for x in ['fold0','foldf','ef','nsim','tmax','trange']])

# helper function to call kfold with pool.map
def kfold_unpack(args):
    return kfold_wrap(*args)

# wrap the kfold lib
# can be called to simulate RNA without multiprocessing
# Inputs
#   seq    :: RNA sequence (string)
#   fold0  :: starting structure to initialize trajectory from (string)
#   foldf  :: structure to monitor for when computing the first passage time distribution (string)
#   ef     :: free energy to monitor for when computing the energy first passage time distribution (float, kcal/mol)
#   nsim   :: number of simulations to run serially WITHIN kfold subroutine (int)
#   tmax   :: total kfold simulation time (float, required to be a multiple of 10.0)
#   trange :: 
# Returns
#   folds  :: list of lists with dot-parantheses folds
#   dGs    :: list of lists with free energies (kcal/mol)
#   fpt    :: list of first passage times (kfold a.u.)
#   efpt   :: list of energy first passage times (kfold a.u.)
def kfold_wrap(seq,fold0,foldf,ef,nsim,tmax,trange):
    traj,e,fpt,efpt = kfold.run(seq,fold0,foldf,ef,nsim,tmax,trange)
    folds = [["".join(traj[i][j][:len(seq)]) for j in xrange(len(trange))] for i in xrange(nsim)]
    dGs   = [[e[i][j] for j in xrange(len(trange))] for i in xrange(nsim)]
    return folds,dGs,fpt,efpt

# Python 3.3 this can be all replaced with pool.starmap() which accepts a sequence of
# argument tuples, and then automatically unpacks the arguments from each tuple; see:
# https://docs.python.org/dev/library/multiprocessing.html#multiprocessing.pool.Pool.starmap
@contextmanager
def poolcontext(*args,**kwargs):
    print "Initializing multiprocessing pool...",
    pool = multiprocessing.Pool(*args, **kwargs)
    print "Ready."
    yield pool
    pool.terminate()

# main function to simulate list of sequences
# Arguments
#   sequences :: list of RNA sequences to simulate
#   options   :: list of dictionary of options for all simulations, see function `get_default_options`
#   optsfxn   :: optional function to return options, by default uses `get_default_options` if <options> argument is None
#   N         :: number of processes to instance for multiprocessing (uses default CPU count otherwise)
def run(sequences,options=None,optsfxn=get_default_options,N=NCPUS):
    
    # check inputs
    assert isinstance(sequences,list), 'Input argument <sequences> for function `simulate_rna` should be a list of sequence(s).'
    if not options is None: assert len(sequences)==len(options), 'Invalid input or size of argument options for KFOLDWrapper. See notes in module.'

    # simulate using python's multiprocessing
    with poolcontext(processes=N) as pool:

        # iterate over list of sequences
        for i,seq in enumerate(sequences):
            assert set(seq)<alphabet, 'Invalid sequence given for kfold: {}. Allowed characters: [ATCGUatcgu].'.format(seq)

            # Set up and validate options input dictionary
            if options is None:
                opts=optsfxn(seq)
            else:
                opts = options[i]
                check_options(seq,opts)

            # Be verbose!
            print """\nSimulating: {} (first 50 nucleotides shown) up to {} kfold time,
            for {} total trajectories, with {} cores...""".format(seq[:50],opts['tmax'],opts['pynsim']*opts['nsim'],N)
            inputs = make_input_tuple(seq,opts)

            # using pool.map with no tqdm
            # results = pool.map(kfold_unpack,repeat(input_tuple,options['pynsim']))
            
            # using tqdm to track progress
            output = dict(structures=[],energies=[],fpts=[],efpts=[])
            output['sequence'] = seq
            output['options']  = opts
            for folds,dGs,fpts,efpts in tqdm.tqdm(pool.imap_unordered(kfold_unpack,repeat(inputs,opts['pynsim'])), total=opts['pynsim']):
                output['structures'].extend(folds)
                output['energies'].extend(dGs)
                output['fpts'].extend(fpts)
                output['efpts'].extend(efpts)

            yield output


def test_wrap():
    seq='ATTCTTAGGGGCGGAGCGGCGCGGCGCCCCTAAGAATTTTT'
    folds,dGs,fpt,efpt = kfold_unpack(make_input_tuple(seq,get_default_options(seq)))
    print folds
    print dGs

def test():
    sequences = [
    'ATTCTTAGGGGCGGAGCGGCGCGGCGCCCCTAAGAATTTTT',
    'AGAGGCGCCAAACAGGGGCCCGGAAACGGGCGCTGTTTTTT'
    ]
    for o in run(sequences):
        pass

if __name__ == "__main__":
    # test_wrap()
    test()