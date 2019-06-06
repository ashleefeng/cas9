#! /usr/bin/env python

"""
Find Cas9 binding sites for HeliCas-FISH. Output sequences between binding sites for subsequent probe design.

X. Ashlee Feng xfeng17@jhu.edu June 5, 2019

"""

import pandas as pd 
import numpy as np 
import argparse
import fasta

def reverse_complement(seq):
    wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'N':'N', 'n':'n'}
    rc = ''
    for i in seq:
        rc = wc[i] + rc
    return rc
    
BD_LEN = 23 # length of PAM + protospacer
    

parser = argparse.ArgumentParser()
parser.add_argument('sequence.fasta', help='Target genomic sequence(s).')
parser.add_argument('N', type=int, help='Number of binding sites.')
parser.add_argument('output', help='Output filename.')
parser.add_argument('-r', action='store_true', help='Search the reverse complementary strand.')
parser.add_argument('-g1', type=float, default=0, help='Range of GC content. Constrain the GC content of Cas9 binding site to be within [g1 g2] Default to [0, 1].')
parser.add_argument('-g2', type=float, default=1)
parser.add_argument('-d1', type=int, default=100, help='Constrain the spacing between adjacent Cas9 binding sites to be within [d1, d2) bp. Default to [100, 200).')
parser.add_argument('-d2', type=int, default=200)
args = vars(parser.parse_args())

seq_fname = args['sequence.fasta']
N = args['N']
out_fname = args['output']
d1 = args['d1']
d2 = args['d2']
g1 = args['g1']
g2 = args['g2']
use_revcomp = args['r']

# add optional inputs here

seq_file = open(seq_fname)
fasta_reader = fasta.FASTAReader(seq_file)

out_file = open(out_fname, 'w')
out_prefix = out_fname.split('.')[0]
sgrna_fname = out_prefix + '_sgRNA.fasta'
sgrna_file = open(sgrna_fname, 'w')

while fasta_reader.isnext():
    seq_name, seq = fasta_reader.next()
    seq_len = len(seq)
    if use_revcomp:
        seq = reverse_complement(seq)
    
    ptr = 0
    prev_ptr = 0
    count = 0
    while (ptr < (len(seq)-BD_LEN)) and (count < 1e7):
        two_mer = seq[ptr:(ptr+2)]
        if (two_mer == 'CC') or (two_mer == 'cc'):
            interval_start = (prev_ptr + BD_LEN)
            interval_end = ptr
            interval = seq[interval_start:interval_end]
            if count > 0:
                if not use_revcomp:
                    out_file.write('>%s_interval_%d_%d_%d\n%s\n'%(seq_name, count, interval_start, interval_end, interval))
                else:
                    out_file.write('>%s_interval_%d_%d_%d\n%s\n'%(seq_name, count, seq_len - interval_start, seq_len - interval_end, interval))
                # print(len(interval), interval)
            sgrna_file.write('>%s_sgRNA_%d_%d\n' %(seq_name, count, ptr))
            sgrna_file.write('%s\n' %seq[ptr:(ptr+23)])
            # print(ptr, ptr - prev_ptr - BD_LEN, seq[ptr:(ptr+BD_LEN)])
                
            prev_ptr = ptr
            ptr = ptr + BD_LEN - 1 + d1
            count += 1
            
        ptr += 1
    
seq_file.close()
sgrna_file.close()


