# Cas9 tools

## Find Cas9 binding sites

Find Cas9 binding sites for HeliCas-FISH. Output sequences between binding sites for subsequent probe design.

```
usage: cas9BDFinder.py [-h] [-r] [-g1 G1] [-g2 G2] [-d1 D1] [-d2 D2] sequence.fasta N output

positional arguments:
  sequence.fasta  target genomic sequence(s).
  N               number of binding sites.
  output          output filename.

optional arguments:
  -h, --help      show this help message and exit
  -r              search the reverse complementary strand.
  -g1 G1          range of GC content. Constrain the GC content of Cas9 binding site to be within [g1 g2]. Default to [0, 1].
  -g2 G2
  -d1 D1          constrain the spacing between adjacent Cas9 binding sites to be within [d1, d2) bp. Default to [100, 200).
  -d2 D2

```

## Design FISH probes with OligoArray

Run OligoArray on regions between Cas9 binding sites.

Filter probes for quality.

## Finalize Cas9 sgRNA sequences

Match Cas9 binding sites with filtered probes