# Fast-re-implementation-of-DMRCaller-s-readONTbam-and-selectCytosine
A fast C++/htslib re-implementation of [DMRcaller](https://bioconductor.org/packages/DMRcaller/)'s `readONTbam()` (and the `selectCytosine()` step before it), which decodes Oxford Nanopore base-modification calls (`MM`/`ML` tags) from BAM files into per-cytosine methylation counts.
