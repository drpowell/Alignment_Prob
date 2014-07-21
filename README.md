FuzzyLZ and AlignCompress (M-alignment)
=========================

This repository contains two distinct but related programs

AlignCompress (M-alignment)
--------

AlignCompress aligns two sequences using the M-alignment algorithm.  The key
contribution of this alignment method is that it takes into account the
information content of the sequences to be aligned.

    L. Allison, D. R. Powell and T. I. Dix
    "Compression and Approximate Matching"
    The Computer Journal, 1999 42:1 pp1-10


There are 2 implementations of this
  * [alignCompress](README.alignCompress) written in java
  * [m-align](alignCompress/C.version/m-align/README) written in C


FuzzyLZ
-------

Compress a sequence (of DNA) using an inexact repeats model.

Please read [README.FuzzyLZ](README.FuzzyLZ).
