#!/usr/bin/env python

import sys
import numpy as np
#from glob import glob
from collections import defaultdict
from scipy import stats

if len(sys.argv) <2:
    print('Usage: compare_count_invariants.py <list of count_invariants files>')

def compute_p_complete(d,msl):
    # probability that the regeion is not destroyed
    s = np.sum(d, axis=1)
    idx = np.where(s != 0)[0]
    return d[idx, 10]/s[idx], msl[idx]

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

    d,msl = compute_p_complete(d,msl)
    
    data[((fam,repeat),sample)] = d
    ms[(fam,repeat),sample] = msl
    
print('\t'.join("famId type pv_all, pv_5 pv_6 pv_7 pv_8 pv_9 pv_10".split(' ')))

for fam in families:
    res = np.zeros((7))
    res[0] = stats.ttest_ind(data[(fam,'blood')],
                                      data[(fam,'tumor')],
                                      equal_var=False).pvalue
    for k in range(5,11):
        id1 = np.where(ms[(fam,'blood')] == k)[0]
        id2 = np.where(ms[(fam,'tumor')] == k)[0]        
        res[k-4] = stats.ttest_ind(data[(fam,'blood')][id1],
                                      data[(fam,'tumor')][id2],
                                      equal_var=False).pvalue
        
    out=f'{fam[0]}\t{fam[1]}\t%g\t%g\t%g\t%g\t%g\t%g\t%g'
    print(out % tuple(res))
    
