# Cas9 tools

## Find Cas9 binding sites

Find Cas9 binding sites for HeliCas-FISH. Output sequences between binding sites for subsequent probe design.

```
./cas9BDFinder.py [-rs] [-d d1 d2] [-g g1 g2] seq.fasta N output

Required:

seq.fasta Target genomic sequence(s).
N         Number of binding sites.
output    Output filename.

Optional:

-r Search the complementary strand. Default to forward strand.
-s Consider sequence dependence of Cas9 activity. Default does not consider sequence dependence.

-d d1 d2  Constrain the spacing between adjacent Cas9 binding sites to be within [x1, x2) bp. Default to [100, 200).
-g g1 g2  Range of GC content. Two floats within [0, 1]. Default to [0, 1].

```

## Design FISH probes with OligoArray

Run OligoArray on regions between Cas9 binding sites.

Filter probes for quality.

## Finalize Cas9 sgRNA sequences

Match Cas9 binding sites with filtered probes