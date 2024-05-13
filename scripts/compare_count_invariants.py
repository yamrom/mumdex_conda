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

families = []

for f in sys.argv[1:]:
    c = f.split('/')[-2].split('_')
    fam = f"{c[0]}_{c[2]}"
    families.append(fam)
    fam_t = c[3]
    d = []
    with open(f, 'r') as F:
        d = np.array([l.strip('\n\r').split('\t')[1:]
                      for l in F.readlines()]).astype(int)

    data[(fam, fam_t)] = compute_p_complete(d)

results = defaultdict()


for fam in families:
    results[(fam)] = (stats.ttest_ind(data[(fam,'blood')],
                                      data[(fam,'tumor')],
                                      equal_var=False),
                      np.mean(data[(fam,'blood')]),
                      np.std(data[(fam,'blood')]),
                      np.mean(data[(fam,'tumor')]),
                      np.std(data[(fam,'tumor')]))

print('\t'.join("famId statistic pvalue df normal_mean normal_std tumor_mean tumor_std".split(' ')))

for k,v in results.items():
    out=f'{k}\t%.4f\t%g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t'
    print(out % (v[0].statistic, v[0].pvalue, v[0].df, v[1], v[2], v[3], v[4]))

    
