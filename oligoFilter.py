#! /usr/bin/env python


import argparse

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Throw out bad sgRNAs based on OligoArray output.')

# Required arguments

parser.add_argument('oligo.txt', help='Output from OligoArray.')
parser.add_argument('sgRNA.fasta', help='sgRNA candidates.')
parser.add_argument('interval.fasta', help='Sequence between sgRNAs that belong to one cluster.')
parser.add_argument('max_nonspecific', type=int, help='')

args = vars(parser.parse_args())

oligo_filename = args['oligo.txt']
sgRNA_filename = args['sgRNA.fasta']
interval_filename = args['interval.fasta']
max_nonspecific = args['max_nonspecific'];

sgRNA_filtered_fname = sgRNA_filename.split('.')[0] + '_filtered.fasta'
interval_filtered_fname = interval_filename.split('.')[0] + '_filtered.fasta'

oligo_file = open(oligo_filename)
sgRNA_file = open(sgRNA_filename)
interval_file = open(interval_filename)
sgRNA_filtered_file = open(sgRNA_filtered_fname, 'w')
interval_filtered_file = open(interval_filtered_fname, 'w')


count1 = 0
count2 = 0

good_sgRNA = set()
good_interval = set()

for line in oligo_file:
    line = line.rstrip('\n')
    cols = line.split('\t')
    interval_name = cols[0]
    binding_sites = cols[7].split(';')
    n_bd = 0
    
    for bd in binding_sites:
        bd_homolog = bd.split(',')
        for bdh in bd_homolog:
            n_bd += 1

    print("Oligo %s has %d binding sites" %(interval_name, n_bd))
    
    if n_bd <= max_nonspecific:

        count1 += 1

        interval_ID = int(interval_name.split('_')[3])
        sgRNA_start = interval_ID
        sgRNA_end = interval_ID + 1

        sgRNA_start_name = '_'.join(interval_name.split('_')[:-4]) + '_sgRNA_' + str(sgRNA_start)
        sgRNA_end_name = '_'.join(interval_name.split('_')[:-4]) + '_sgRNA_' + str(sgRNA_end)
        print sgRNA_start_name
        
        if interval_name not in good_interval:
            good_interval.add(interval_name)

        if sgRNA_start_name not in good_sgRNA:

            good_sgRNA.add(sgRNA_start_name)
            
        if sgRNA_end_name not in good_sgRNA:

            good_sgRNA.add(sgRNA_end_name)

    else:
        count2 += 1

print("\n%d oligos are good" %count1)
print("%d oligos are bad\n" %count2)

oligo_file.close()

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
print("Thowing out %d sgRNAs" %bad_count)

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
print("Thowing out %d intervals" %bad_count)

interval_filtered_file.close()
interval_file.close()



















