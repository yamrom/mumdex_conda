#!/usr/bin/env python

import sys
import numpy as np
#from glob import glob
from collections import defaultdict
from scipy import stats

if len(sys.argv) <2:
    print('Usage: compare_count_invariants.py <list of count_invariants files>')


def compute_p_complete(d, msl):
    # probability that the region is not destroyed
    s0 = np.sum(d[0], axis=1)
    s1 = np.sum(d[1], axis=1)
    idx = np.where(s0*s1 != 0)[0]
    return d[0][idx, 10]/s0[idx],d[1][idx, 10]/s1[idx],msl[idx]

data = defaultdict()
ms = defaultdict() 
families = set()

def diff (x):
    return x[1]-x[0]+1

for f in sys.argv[1:]:
    c = f.split('/')[-2].split('_')
    fam = f"{c[0]}_{c[1]}_{c[2]}"
    repeat =  f.split('/')[-1].split('.')[0]
    families.add((fam,repeat))
    sample = c[3]
    d = []
    
    with open(f, 'r') as F:
        d = [l.strip('\n\r').split('\t')
                      for l in F.readlines()]
        
    msl = [diff(list(map(int,y[0].split(':')[1].split('-')))) for y in d]
    msl = np.array(msl)
    assert min(msl) == 5
    assert max(msl) == 10
    d = np.array([x[1:] for x in d]).astype(int)
    
    data[((fam,repeat),sample)] = d
    ms[(fam,repeat)] = msl
    
print('\t'.join("famId type int_pv_all, ind_pv_5 ind_pv_6 ind_pv_7 ind_pv_8 ind_pv_9 ind_pv_10 rel_pv_all, rel_pv_5 rel_pv_6 rel_pv_7 rel_pv_8 rel_pv_9 rel_pv_10".split(' ')))

for fam in sorted(families):
    res = np.zeros((14))
    d_blood, d_tumor, msl = compute_p_complete([data[(fam,'blood')],
                                                data[(fam,'tumor')]],
                                               ms[(fam)])
    
    res[0] = stats.ttest_ind(d_blood, d_tumor,
                                      equal_var=False).pvalue
    for k in range(5,11):
        id = np.where(msl == k)[0]
        res[k-4] = stats.ttest_ind(d_blood[id], d_tumor[id],
                                      equal_var=False).pvalue

    res[7] = stats.ttest_rel(d_blood, d_tumor).pvalue

    for k in range(5,11):
        id = np.where(msl == k)[0]
        res[k+3] = stats.ttest_rel(d_blood[id], d_tumor[id]).pvalue
        
    out=f'{fam[0]}\t{fam[1]}\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g'
    print(out % tuple(res))
    
