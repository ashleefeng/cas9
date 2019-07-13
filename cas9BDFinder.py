#! /usr/bin/env python

"""
Find Cas9 binding sites for HeliCas-FISH. 

Output: file.fasta - sequences between binding sites for subsequent probe design.
        file_sgRNA.fasta - sgRNA sequences

X. Ashlee Feng xfeng17@jhu.edu July 7, 2019

"""

import pandas as pd 
import numpy as np 
import argparse
import fasta
import unittest
import sys

# Utility functions

"""
Find reverse complement of a given DNA sequence

input: string, DNA sequence
output: string, reverse complement of given sequence

"""
def reverse_complement(seq):
    wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'N':'N', 'n':'n'}
    rc = ''
    for i in seq:
        rc = wc[i] + rc
    return rc


"""
Determine GC content of a given DNA sequence

input: string, DNA sequence
output: float, fraction of GC in the given sequence

"""   
def gc_content(seq):
    count = 0
    for n in seq:
        if n in 'GgCc':
            count += 1
    return float(count)/len(seq)

"""
Determine whether a protospacer is bad

input: protospacer - string, DNA sequence of protospacer
       g1 - float, min gc content
       g2 - float, max gc content

output: boolean, True if the protospacer is bad

"""
def is_bad(protospacer, g1=0, g2=1):
    
    # Check affinity

    if ('N' in protospacer) or ('n' in protospacer):
        return True
    # 3rd nucleotide is C on guide RNA
    if protospacer[19-2] in 'Gg':
        return True
    # 17th nucleotide is U on guide RNA
    if protospacer[19-16] in 'Aa':
        return True
    # 18th nucleotide is U on guide RNA
    # if protospacer[19-17] in 'Aa':
    #     return True
    # 19th nucleotide is U on guide RNA
    if protospacer[19-18] in 'Aa':
        return True
    
    # 20th is C and U 
    # if protospacer[19-19] in 'AaGg':
    #     return True
    
    # 20th U
    if protospacer[19-19] in 'Aa':
        return True

    # Check gc content

    if (gc_content(protospacer) < g1) or (gc_content(protospacer) > g2):
        print("gc content")
        return True
    
    else:
        return False

    return False


"""
Find a cluster of cas9 binding sites

input: seq - string, DNA sequence
       start - integer, where search should start from
       clusterMinSize - integer, minimum size of cluster
       g1 - float, min gc content
       g2 - float, max gc content
       d1 - integer, min distance between two binding sites
       d2 - interger, max distance between two binding sites
       n - max number of binding sites to find

output: list of integer - starting indices for cas9 binding sites
        integer - pointer to after the last binding site, or where the next search should start

"""
def find_cluster(seq, start, clusterMinSize, g1, g2, d1, d2, n):

  # initialize parameters

    buff = []
    curr_d = 0

    ptr = start
    prev_ptr = start

    while (curr_d <= d2) and (len(buff) < n) and (ptr <= (len(seq)-BD_LEN)):

        # Find PAM site

        two_mer = seq[ptr:(ptr+2)]
        
        if (two_mer == 'CC') or (two_mer == 'cc'):
            
            protospacer = seq[(ptr+3):(ptr+23)]
            
            # Skip this PAM site if this protospacer is bad
            
            if is_bad(protospacer):
                
                ptr += 1
                print("Protospacer %d is bad: %s, gc content %.2f" %(ptr, protospacer, gc_content(protospacer)))
                continue

            buff.append(ptr)
                
            # Update parameters
            
            curr_d = ptr - prev_ptr - BD_LEN
            prev_ptr = ptr
            ptr = ptr + BD_LEN + d1
            
        else:
            
            curr_d += 1
            ptr += 1

    # Check if cluster is big enough

    if len(buff) < clusterMinSize:

       return [], -1

    return buff, ptr

"""
Main function for identifying cas9 binding sites

input: seq - string, DNA sequence
       g1 - float, min gc content
       g2 - float, max gc content
       d1 - integer, min distance between two binding sites
       d2 - interger, max distance between two binding sites
       n - max number of binding sites to find
       clusterMinSize - interger, min number of elements in the cluster

output: list of integers, starting indices for cas9 binding sites

"""     
def get_sgRNA(seq, g1, g2, d1, d2, n, clusterMinSize):
    
    # Initialize parameters
    
    ptr = 0
    sgRNA_list = []
    printed = False
    
    while (ptr <= (len(seq)-BD_LEN)) and (len(sgRNA_list) <= n):

    	if ptr > 2000 and not printed:
    		print("Searching at %dth nucleotide." %ptr)
    		printed = True
        
        cluster, ptr = find_cluster(seq, ptr, clusterMinSize, g1, g2, d1, d2, n - len(sgRNA_list))
        print(ptr)
        
        # found cluster

        if ptr != -1:

            sgRNA_list = sgRNA_list + cluster

        # no cluster found in current iteration

        else:
            
            break
    
    return sgRNA_list


"""
Save output

input: seq - string, DNA sequence
       sgRNA_indices - list of intergers, output from get_sgRNA
       out_file - file, where DNA sequence between cas9 binding sites will be stored
       use_revcomp - boolean, True if search was performed on reverse complementary sequence
       seq_name - string, parsed from fasta file
       sgrna_file - file, where guide RNA sequences will be stored

"""
def output_writer(seq, sgRNA_indices, out_file, use_revcomp, seq_name, sgrna_file):
    
    count = 0
    prev_sg_end = -1
    
    for sg_start in sgRNA_indices:
        
        sg_end = sg_start + BD_LEN
        
        if count > 0:
            
            interval_start = prev_sg_end
            interval_end = sg_start
            interval = seq[interval_start:interval_end]
            
            if not use_revcomp:
                
                out_file.write('>%s_interval_%d_%d_%d\n%s\n'%(seq_name, count, interval_start, interval_end, interval))
            
            else:
                out_file.write('>%s_interval_%d_%d_%d\n%s\n'%(seq_name, count, seq_len - interval_start, seq_len - interval_end, interval))

            sgrna_file.write('>%s_sgRNA_%d_%d\n' %(seq_name, count + 1, sg_start))
            sgrna_file.write('%s\n' %seq[sg_start : sg_end])
        
        count += 1
        prev_sg_end = sg_end
    

"""
Unit tests

"""
class TestMethods(unittest.TestCase):

    dna1 = "ATCGNatcgn"
    dna2 = "CCATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAACCATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAACCATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAACCATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAACCATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAACCATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAACCATTTTTTTTTTAAAAAAAAAAGG"
    dna3 = "AAAAATCGncgAagCGATNA"
    dna4 = "CCCTATCGAcgAagCGATaA"
    dna5 = "CCCCCCCGAcgAagCGATaA"
    dna6 = "CCATTTTTTTTTTAAAAAAAAAAGGGGGCCATTTTTTTTTTAAAAAAAAAAGGGGGGGGGGCCATTTTTTTTTTAAAAAAAAAA"


    def test_reverse_complement(self):
        self.assertEqual(reverse_complement(self.dna1), "ncgatNCGAT")

    def test_gc_content(self):
        self.assertEqual(gc_content(self.dna1), 0.4)

    def test_is_bad(self):
        self.assertTrue(is_bad(self.dna3))
        self.assertFalse(is_bad(self.dna4))
        self.assertTrue(is_bad(self.dna5, g1=0, g2=0.5))

    def test_find_cluster(self):
        self.assertEqual(find_cluster(seq=self.dna2, start=0, clusterMinSize=2, g1=0, g2=1, d1=100, d2=10000, n=3), ([0, 123, 246], 369))
        self.assertEqual(find_cluster(seq=self.dna2, start=0, clusterMinSize=2, g1=0, g2=1, d1=100, d2=10000, n=7), ([0, 123, 123*2, 123*3, 123*4, 123*5, 123*6], 123*7))
        
        self.assertEqual(find_cluster(seq=self.dna6, start=0, clusterMinSize=3, g1=0, g2=1, d1=5, d2=10, n=3), ([0, 28, 61], 89))
        self.assertEqual(find_cluster(seq=self.dna6, start=0, clusterMinSize=1, g1=0, g2=1, d1=5, d2=9, n=3), ([0, 28], 61))
        self.assertEqual(find_cluster(seq=self.dna6, start=61, clusterMinSize=1, g1=0, g2=1, d1=5, d2=9, n=3), ([61], 89))
        self.assertEqual(find_cluster(seq=self.dna6, start=0, clusterMinSize=3, g1=0, g2=1, d1=5, d2=9, n=3), ([], -1))
        self.assertEqual(find_cluster(seq=self.dna6, start=23, clusterMinSize=3, g1=0, g2=1, d1=5, d2=10, n=8), ([], -1))

    def test_get_sgRNA(self):
        self.assertEqual(get_sgRNA(seq=self.dna2, g1=0, g2=1, d1=100, d2=1000, n=3, clusterMinSize=2), [0, 123, 246])
        self.assertEqual(get_sgRNA(seq=self.dna6, g1=0, g2=1, d1=5, d2=9, n=3, clusterMinSize=3), [])
        self.assertEqual(get_sgRNA(seq=self.dna6, g1=0, g2=1, d1=5, d2=9, n=3, clusterMinSize=1), [0, 28, 61])
        self.assertEqual(get_sgRNA(seq=self.dna6, g1=0, g2=1, d1=5, d2=9, n=3, clusterMinSize=2), [0, 28])

# Constants
    
BD_LEN = 23 # length of PAM + protospacer
CLUSTER_MIN_SIZE = 10 # minimum size of binding site cluster


if __name__ == '__main__':

    # Run unit tests if no input is provided

    if len(sys.argv) == 1:
        unittest.main()
        quit()
        
    # Parse arguments

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(description='Find Cas9 binding sites for HeliCas-FISH. Run tests if no input is provided.')

    # Required arguments

    parser.add_argument('sequence.fasta', help='target genomic sequence(s).')
    parser.add_argument('N', type=int, help='max number of binding sites needed.')
    parser.add_argument('output', help='output filename.')

    # Optional arguments

    parser.add_argument('-r', action='store_true', help='search the reverse complementary strand.')
    parser.add_argument('-g1', type=float, default=0.35, help='min GC content of protospacers. Default to 0.35.')
    parser.add_argument('-g2', type=float, default=0.75, help='max GC content of protospacers. Default to 0.75.')
    parser.add_argument('-d1', type=int, default=50, help='min distance between binding sites. Default to 50.')
    parser.add_argument('-d2', type=int, default=200, help='max distance between binding sites. Default to 200.')
    parser.add_argument('-clusterMinSize', type=int, default=1, help='minimum number of Cas9 binding sites in one cluster. Default to 1.')

    args = vars(parser.parse_args())

    seq_fname = args['sequence.fasta']
    n = args['N']
    out_fname = args['output']
    d1 = args['d1']
    d2 = args['d2']
    g1 = args['g1']
    g2 = args['g2']
    use_revcomp = args['r']
    clusterMinSize = args['clusterMinSize']

    # Read input files

    seq_file = open(seq_fname)
    fasta_reader = fasta.FASTAReader(seq_file)

    out_file = open(out_fname, 'w')
    out_prefix = out_fname.split('.')[0]
    sgrna_fname = out_prefix + '_sgRNA.fasta'
    sgrna_file = open(sgrna_fname, 'w')

    while fasta_reader.isnext():
        
        # Get name and sequence of fasta entry
        
        seq_name, seq = fasta_reader.next()
        seq_len = len(seq)

        print("Working on %s" %seq_name)
        
        # Reverse complementary strand option
        
        if use_revcomp:
            seq = reverse_complement(seq)
            
        sgRNA_IDs = get_sgRNA(seq, g1, g2, d1, d2, n, clusterMinSize)

        if len(sgRNA_IDs) == 0:
            print("No good cas9 binding site was found for %s" %seq_name)
            continue

        else:
            print("Found %d cas9 binding sites" %len(sgRNA_IDs))

        output_writer(seq, sgRNA_IDs, out_file, use_revcomp, seq_name, sgrna_file)
        
    seq_file.close()
    sgrna_file.close()

