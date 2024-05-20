#!/usr/bin/env python

import sys
import numpy as np
#from glob import glob
from collections import defaultdict
from scipy import stats

if len(sys.argv) <2:
    print('Usage: compare_count_invariants.py <list of count_invariants files>')

def compute_p_complete(d):
    # probability that the regeion is not destroyed
    s = np.sum(d, axis=1)
    idx = np.where(s != 0)[0]
    return d[idx, 10]/s[idx]

data = defaultdict()

families = set()

for f in sys.argv[1:]:
    c = f.split('/')[-2].split('_')
    fam = f"{c[0]}_{c[1]}_{c[2]}"
    repeat =  f.split('/')[-1].split('.')[0]
    families.add((fam,repeat))
    sample = c[3]
    d = []
    with open(f, 'r') as F:
        d = np.array([l.strip('\n\r').split('\t')[1:]
                      for l in F.readlines()]).astype(int)

    data[((fam,repeat),sample)] = compute_p_complete(d)

print('\t'.join("famId type statistic pvalue df".split(' ')))

for fam in families:
    res = stats.ttest_ind(data[(fam,'blood')],
                                      data[(fam,'tumor')],
                                      equal_var=False)
    out=f'{fam[0]}\t{fam[1]}\t%.4f\t%g\t%.1f'
    print(out % (res.statistic, res.pvalue, res.df))
    
