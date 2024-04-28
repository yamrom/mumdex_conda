#! /usr/bin/env python
import sys
import mumdex
import numpy as np
from colorama import Fore
from collections import defaultdict

def complete_overlap(window_size,
                     pos1,
                     pos2,
                     pos,
                     offset,
                     length):
    # position in mums is reported as 1-based
    if pos1+window_size >= pos-1 and pos2-window_size <= pos-1 + length:
            return [True, False, False]
    else:
        low_in = pos1 < pos-1 and pos-1 < pos2
        high_in = pos-1 + length  > pos1 and pos-1  + length < pos2
        return [False, low_in, high_in]


def compute_inv(m1,m2):
    return m1[0]-m1[1]-(m2[0]-m2[1])
    
        
# Parse command line arguments
mums = mumdex.MUMdex(sys.argv[1])
ref = mums.Reference()
bed = []
window_size = 10

if len(sys.argv) > 2:
    with open(sys.argv[2], 'r') as F:
        for l in F:
            chrom, pos1, pos2 = l.strip('\n\r').split('\t')
            bed.append((chrom, int(pos1), int(pos2)))

data = np.zeros((len(bed), 21),dtype=int)
pattern = ""
if len(sys.argv) > 3:
    pattern = sys.argv[3]
    
buffer = []
for ind, v in enumerate(bed):    
    chrom = v[0]
    pos1 = v[1]-window_size
    pos2 = v[2]+window_size
    r_string=f"{chrom}:{str(pos1)}-{str(pos2)}"

    region = mums.range_str(r_string)
    
    # Avoid double-viewing pairs
    pairs = set()

    # Loop over mums in range
    for i in range(region[0], region[1]):

        # Look at pair for mum
        pm = mums.index(i)
        p = pm[0]

        # Only if pair not yet seen
        if not p in pairs:
            pairs.add(p)

            pair = mums.Pair(p)

            dupe = pair.dupe()
            bad_1 = pair.read_1_bad()
            bad_2 = pair.read_2_bad()
            if dupe: continue
        
            read_mums = defaultdict(list)
        
            # Loop over mums in pair
            for m in range(pair.mums_start(), pair.mums_stop()):
                mum = mums.MUM(m)

                if ref.name(mum.chromosome()) == chrom:
                    pos = mum.position()
                    read_2 = mum.read_2()
                    offset = mum.offset()
                    length = mum.length()
                    flipped = mum.flipped()
                    last_hit = mum.last_hit()
                    touches_end = mum.touches_end()
                    read_mums[(read_2,flipped)].append((pos,
                                                    offset,
                                                    length,
                                                    last_hit,
                                                    touches_end))

            for i in range(2):
                if (not i and bad_1) or (i and bad_2):
                    continue
                for j in range(2):
                    read_mums[i,j] = sorted(read_mums[i,j], key=lambda x: x[1])
                    n_mums = len(read_mums[i,j])

                    if n_mums:
                        for n,m in enumerate(read_mums[i,j]):
                            over, low, high = complete_overlap(window_size,
                                                               pos1,
                                                               pos2,
                                                               m[0],
                                                               m[1],
                                                               m[2])

                            if over:
                                data[ind, 10] += 1

                            if low and n > 0:
                                inv = compute_inv(read_mums[i,j][n-1],m)
                                if inv > -11 and inv < 11:
                                    data[ind,10+inv] += 1

                            if high and n < n_mums-1:
                                inv = compute_inv(m,read_mums[i,j][n+1])
                                if inv > -11 and inv < 11:
                                    data[ind, 10+inv] += 1

np.save(f"{pattern}.npy", data)

