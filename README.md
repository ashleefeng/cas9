# Cas9 tools

## Find Cas9 binding sites

Find Cas9 binding sites for HeliCas-FISH. Outputs `input_intervals.fa` which contains sequences between binding sites for subsequent probe design, and `input_sgRNAs.fa` contains guide RNA sequences.

```
usage: cas9BDFinder.py [-h] [-r] [-g1 G1] [-g2 G2] [-d1 D1] [-d2 D2] [-clusterMinSize CLUSTERMINSIZE] sequence.fasta N

Find Cas9 binding sites for HeliCas-FISH. Run tests if no input is provided.

positional arguments:
  sequence.fasta        target genomic sequence(s).
  N                     max number of binding sites needed.

optional arguments:
  -h, --help            show this help message and exit
  -r                    search the reverse complementary strand.
  -g1 G1                min GC content of protospacers. Default to 0.35.
  -g2 G2                max GC content of protospacers. Default to 0.75.
  -d1 D1                min distance between binding sites. Default to 50.
  -d2 D2                max distance between binding sites. Default to 200.
  -clusterMinSize CLUSTERMINSIZE
                        minimum number of Cas9 binding sites in one cluster.
                        Default to 1.
```

## Design FISH probes with OligoArray

Run OligoArray on regions between identified Cas9 binding sites.

## Finalize Cas9 sgRNA sequences

Filtered the FISH probes for good quality and discard Cas9 sites flanking bad probes.

```
usage: oligoFilter.py [-h] oligo.txt sgRNA.fasta interval.fasta max_nonspecific num_probes

Throw out bad sgRNAs and FISH probes based on OligoArray output. Outputs
oligos_filtered.tsv, sgRNAs_filtered.fasta, and intervals_filtered.fasta

positional arguments:
  oligo.txt        Output from OligoArray.
  sgRNA.fasta      sgRNA candidates.
  interval.fasta   Sequence between sgRNAs that belong to one cluster.
  max_nonspecific  Maximum number of nonspecific sites allowed for one FISH
                   probe
  num_probes       Maximum number of probes needed for one region

optional arguments:
  -h, --help       show this help message and exit
```
