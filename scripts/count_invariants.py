#! /usr/bin/env python
import sys
import mumdex
import numpy as np
from collections import defaultdict
import gzip

"""
The gist of the algorithm.

For each site get list of mum's overlaying the site.
Each mum has associated pair
and an indicator if it is dupe. 
Each read has indicator that it is bad.
We ignore duped and bad reads.
Each read is processed once.
Each read has interval of indices of mums that are on the read.
We collect mums sorted by there offset on the read into list with the key:
with key= (read, flipped).
We go through this list of mums and compute over, low_in, and high_in.
If over, we accumulate invariant 0, 
if lower_in and it is not the first mum on read, compute and update invariants.
if high_in and it is not the last mum on read, compute and update invariants.
Ignore invariants that are outside [-10,10] range.

"""

# Parse command line arguments
if len(sys.argv) > 2:
    mums = mumdex.MUMdex(sys.argv[1])
    ref = mums.Reference()
    bed = []
    window_size = 10
    bed_file = sys.argv[2]
    if bed_file.endswith("gz"):
        with gzip.open(bed_file, 'rt') as F:
            for l in F:
                chrom, pos1, pos2 = l.strip('\n\r').split('\t')
                bed.append((chrom, int(pos1), int(pos2)))
    else:
        with open(bed_file, 'r') as F:
            for l in F:
                chrom, pos1, pos2 = l.strip('\n\r').split('\t')
                bed.append((chrom, int(pos1), int(pos2)))        
else:
    print("Usage: count_invariants_bed.py <mumdex dir> <bed file> [<pattern>, <log, default = 0>]")
    
data = np.zeros((len(bed), 21)).astype(int)
pattern = ""
if len(sys.argv) > 3:
    pattern = sys.argv[3]
log = 0
if len(sys.argv) >4:
    log = sys.argv[4]
    
def complete_overlap(window_size:int,
                     pos1:int,
                     pos2:int,
                     pos:int,
                     offset:int,
                     length: int) -> [bool, bool, bool]:
    """ Input: window_size - extension of the site to the left and right,
    pos1,pos2 - the start and the end after extension
    pos - mum's position, length - mum's length
    return: over - mum completelys overlay the site, 
            low_in - mum has the left end in [pos1,pos2], 
            high_in - mum has the right end in [pos1,pos2]
    comment: mum's position is 1-based
    """
    
    # position in mums is reported as 1-based
    if pos1+window_size >= pos-1 and pos2-window_size <= pos-1 + length:
            return [True, False, False]
    else:
        low_in = pos1 < pos-1 and pos-1 < pos2
        high_in = pos-1 + length  > pos1 and pos-1  + length < pos2
        return [False, low_in, high_in]


def compute_inv(m1,m2) -> int:
    """ Returns the bridge invariant of two mums equally oriented on one read
    """
    return m1[0]-m1[1]-(m2[0]-m2[1])
    
for ind, v in enumerate(bed):    
    chrom = v[0]
    pos1 = v[1]-window_size
    pos2 = v[2]+window_size
    r_string=f"{chrom}:{str(pos1)}-{str(pos2)}"
    if log:
        print(f"r_string: {r_string}", file=sys.stderr)
        
    region = mums.range_str(r_string)
    
    # Avoid double-viewing pairs

    pairs = set([mums.index(i)[0] for i in range(region[0], region[1])])

    for p in pairs:
        pair = mums.Pair(p)

        dupe = pair.dupe()
        bad_1 = pair.read_1_bad()
        bad_2 = pair.read_2_bad()
        if dupe:
            continue
        
        read_mums = defaultdict(list)
        
        # Loop over mums in pair and
        # accumulate them in array indexed by (read,flipped)
        
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
                # sort mums according their offset in read 
                read_mums[i,j] = sorted(read_mums[i,j], key=lambda x: x[1])
                n_mums = len(read_mums[i,j])
                
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

for b,d in zip(bed,data):
    r_string=f"{b[0]}:{str(b[1])}-{str(b[2])}"
    print('\t'.join([r_string]+list(map(str,d))))
