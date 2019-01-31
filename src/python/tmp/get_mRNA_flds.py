import sys
sys.path.append('..')
from PyVRNA import PyVRNA
ViennaRNA = PyVRNA(pyindex=True)

import pandas as pd

df = pd.read_csv('JACS_2017.csv')
folds = df['final_structure']

final_mRNA_stuctures = []
for f in folds:
    print f, len(f)
    bp_tuple = ViennaRNA.vienna2bp(f)
    # print bp_tuple

    bpx = []; bpy = [];
    for i,y in enumerate(bp_tuple.bpy):
        if not y > bp_tuple.length[0]:
            bpx.append(bp_tuple.bpx[i])
            bpy.append(y)

    vienna = ViennaRNA.bp2vienna(length=[bp_tuple.length[0]], bpx=bpx, bpy=bpy)
    print vienna, len(vienna)
    assert len(vienna) == bp_tuple.length[0]
    final_mRNA_stuctures.append(vienna)

df['final_mRNA_structure'] = final_mRNA_stuctures
df.to_csv('JACS_2017_2.csv')