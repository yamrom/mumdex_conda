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
    s0 = np.sum(d[0], axis=1)
    s1 = np.sum(d[1], axis=1)
    idx = np.where(s0*s1 != 0)[0]
    return d[0][idx, 10]/s0[idx],d[1][idx, 10]/s1[idx]

data = defaultdict()
families = set()

for f in sys.argv[1:]:
    c = f.split('/')[-2].split('_')
    fam = f"{c[0]}_{c[1]}_{c[2]}"
    repeat =  f.split('/')[-1].split('.')[0]
    families.add((fam,repeat))
    sample = c[3]
    with open(f, 'r') as F:
        d = np.array([l.strip('\n\r').split('\t')[1:]
                      for l in F.readlines()]).astype(int)

    data[((fam,repeat),sample)] = d

print('\t'.join("famId type statistic ttest_ind_pvalue df ttest_rel_pvalue".split(' ')))

for fam in sorted(families):
    a,b = compute_p_complete([data[(fam,'blood')], data[(fam,'tumor')]])
    res1 = stats.ttest_ind(a, b, equal_var=False)
    res2 = stats.ttest_rel(a, b)    
    out=f'{fam[0]}\t{fam[1]}\t%.4f\t%g\t%.1f\t%g'
    print(out % (res1.statistic, res1.pvalue, res1.df, res2.pvalue))
    
