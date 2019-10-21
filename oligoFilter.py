#! /usr/bin/env python

import argparse
from fasta import FASTAReader

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Throw out bad sgRNAs based on OligoArray output.')

# Required arguments

parser.add_argument('oligo.txt', help='Output from OligoArray.')
parser.add_argument('sgRNA.fasta', help='sgRNA candidates.')
parser.add_argument('interval.fasta', help='Sequence between sgRNAs that belong to one cluster.')
parser.add_argument('max_nonspecific', type=int, help='')
parser.add_argument('num_probes', type=int, help='')

args = vars(parser.parse_args())

oligo_filename = args['oligo.txt']
sgRNA_filename = args['sgRNA.fasta']
interval_filename = args['interval.fasta']
max_nonspecific = args['max_nonspecific']
num_probes = args['num_probes']

oligo_filtered_fname = oligo_filename.split('.')[0] + '_filtered.tsv'
sgRNA_filtered_fname = sgRNA_filename.split('.')[0] + '_filtered.fasta'
interval_filtered_fname = interval_filename.split('.')[0] + '_filtered.fasta'

oligo_file = open(oligo_filename)
sgRNA_file = open(sgRNA_filename)
interval_file = open(interval_filename)
oligo_filtered_file = open(oligo_filtered_fname, 'w')
sgRNA_filtered_file = open(sgRNA_filtered_fname, 'w')
interval_filtered_file = open(interval_filtered_fname, 'w')



tad2good = {}
tad2bad = {}

good_sgRNA = set()
good_interval = set()
line_num = 0
prob_num = 0

for line in oligo_file:
    line_num += 1
    line = line.rstrip('\n')
    cols = line.split('\t')
    interval_name = cols[0]
    if len(cols) < 8:
        print line_num
        print cols
        continue
    binding_sites = cols[7].split(';')

    n_bd = 0
    bind2rna = False
    
    for bd in binding_sites:
        bd_homolog = bd.split(',')
        for bdh in bd_homolog:
            n_bd += 1

            if ("ncrna" in bdh) or ("Escherichia" in bdh):

                bind2rna = True

                break

            if n_bd > max_nonspecific:

                break
        
        if bind2rna or (n_bd > max_nonspecific):

            break



    # print("Oligo %s has %d binding sites" %(interval_name, n_bd))
    

    tokens = interval_name.split('_')

    tad_name = tokens[0]

    if (n_bd <= max_nonspecific) and not bind2rna:

        if (tad_name not in tad2good) or (tad2good[tad_name] < num_probes):

            prob_num += 1
            TAD_ID = tokens[0]
            interval_ID = int(tokens[-3])
            probe_seq = cols[-1]

            oligo_filtered_file.write('>probe_%d_%s_int_%d\n' %(prob_num, TAD_ID, interval_ID))
            oligo_filtered_file.write('%s\n' %probe_seq)
            
            sgRNA_start = interval_ID
            sgRNA_end = interval_ID + 1

            sgRNA_start_name = '_'.join(interval_name.split('_')[:-4]) + '_sgRNA_' + str(sgRNA_start)
            sgRNA_end_name = '_'.join(interval_name.split('_')[:-4]) + '_sgRNA_' + str(sgRNA_end)
            # print sgRNA_start_name
            
            if interval_name not in good_interval:
                good_interval.add(interval_name)

            if sgRNA_start_name not in good_sgRNA:

                good_sgRNA.add(sgRNA_start_name)
                
            if sgRNA_end_name not in good_sgRNA:

                good_sgRNA.add(sgRNA_end_name)

            if tad_name not in tad2good:

                tad2good[tad_name] = 1

            else:

                tad2good[tad_name] += 1

    else:

        if tad_name not in tad2bad:

            tad2bad[tad_name] = 1

        else:

            tad2bad[tad_name] += 1

for tad in sorted(tad2good.keys()):

    print("%s has %d good probes." %(tad, tad2good[tad]))

    if tad in tad2bad.keys():

        print("%s has %d bad probes." %(tad, tad2bad[tad]))

    print("")

oligo_file.close()
oligo_filtered_file.close()

iswriting = False
good_count = 0
bad_count = 0

for l in sgRNA_file:

    if iswriting:

        sgRNA_filtered_file.write(l)
        iswriting = False
        continue

    line = l.rstrip('\n')

    if line.startswith('>'):

        seq_name = line.lstrip('>')
        sgRNA_ID = '_'.join(seq_name.split('_')[:-1])
        # print sgRNA_ID
        if sgRNA_ID in good_sgRNA:
            
            good_count += 1
            sgRNA_filtered_file.write(l)
            iswriting = True

        else:

            bad_count += 1

print("Keeping %d sgRNAs" %good_count)
print("Throwing out %d sgRNAs" %bad_count)

sgRNA_filtered_file.close()
sgRNA_file.close()

iswriting = False
good_count = 0
bad_count = 0

for l in interval_file:

    if iswriting:
        interval_filtered_file.write(l)
        iswriting = False
        continue

    line = l.rstrip('\n')

    if line.startswith('>'):

        interval_name = line.lstrip('>')

        if interval_name in good_interval:
            
            good_count += 1
            interval_filtered_file.write(l)
            iswriting = True

        else:

            bad_count += 1

print("\nKeeping %d intervals" %good_count)
print("Throwing out %d intervals" %bad_count)

interval_filtered_file.close()
interval_file.close()


# previously TADcounter.py

sgRNA_file = open(sgRNA_filtered_fname)
interval_file = open(interval_filtered_fname)

print('sgRNA:')
sgRNA_reader = FASTAReader(sgRNA_file)

tad2sgRNA = {}

while sgRNA_reader.isnext:

    name, seq = sgRNA_reader.next()

    tadID = name.split('_')[0]

    if tadID not in tad2sgRNA:

        tad2sgRNA[tadID] = 1

    else:

        tad2sgRNA[tadID] += 1

for i in sorted(tad2sgRNA.keys()):

    print("TAD %s has %d sgRNAs" %(i, tad2sgRNA[i]))

print('\ninterval:')
interval_reader = FASTAReader(interval_file)

tad2interval = {}

while interval_reader.isnext:

    name, seq = interval_reader.next()

    tadID = name.split('_')[0]

    if tadID not in tad2interval:

        tad2interval[tadID] = 1

    else:

        tad2interval[tadID] += 1

for i in sorted(tad2interval.keys()):

    print("TAD %s has %d intervals" %(i, tad2interval[i]))

sgRNA_file.close()
interval_file.close()






















