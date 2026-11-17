# readontbam_fast

A fast C++/htslib re-implementation of
[DMRcaller](https://bioconductor.org/packages/DMRcaller/)'s `readONTbam()` (and the
`selectCytosine()` step before it), which decodes Oxford Nanopore base-modification
calls (`MM`/`ML` tags) from BAM files into per-cytosine methylation counts.

* **Same output as DMRcaller.** The R helper turns the decoder output into a GRanges that
  is `identical()` to what `readONTbam()` returns, including its quirks (see
  [Behaviour](#behaviour-reproduced-from-readontbam)), so it can be dropped into an
  existing DMRcaller workflow.
* **Fast.** About 0.8 s per BAM per Mb of genome at ~37× coverage on one core, versus
  ~3 s *per read* for `readONTbam()` in R, i.e. roughly 10,000× faster per core.
  Whole chromosomes or whole genomes run in one go; no windowing and stitching.
* **Light.** Results are streamed to disk; memory stays at ~100–150 MB per BAM being
  decoded, plus 4 bytes per base of the current chromosome.

```
readontbam_fast/
├── src/readontbam_fast.cpp    the decoder (single file, needs htslib)
├── Makefile
├── R/readONTbam_fast.R        TSV -> GRanges, an R wrapper, and a multi-chromosome loader
├── test/selftest.R            identical() check against DMRcaller::readONTbam()
├── test/make_test_data.sh     cut a small test dataset from your own BAM
└── LICENSE                    GPL-3
```

## Install

You need a C++17 compiler and htslib (≥ 1.10).

```bash
# htslib from conda/mamba
mamba create -n readontbam -c conda-forge -c bioconda htslib gxx_linux-64 make
mamba activate readontbam
make HTSLIB_PREFIX=$CONDA_PREFIX

# htslib from your system (e.g. apt install libhts-dev) or a module
make                                   # or: make HTSLIB_PREFIX=/path/to/htslib

# htslib from source, no root needed
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar xjf htslib-1.21.tar.bz2 && cd htslib-1.21
./configure --prefix=$HOME/htslib --disable-libcurl && make -j8 && make install && cd ..
make HTSLIB_PREFIX=$HOME/htslib
```

The binary is `bin/readontbam_fast`; `make install PREFIX=$HOME/.local` copies it to
`$HOME/.local/bin`. The R part needs GenomicRanges (and uses data.table when installed);
DMRcaller itself is only needed for the self-test.

## Usage

### Command line

```bash
# whole genome (every sequence in ref.fa.fai)  -> out/<chrom>/<label>_CG.tsv
readontbam_fast -f ref.fa -x CG -m C+m. -p 0.7 -t 16 -o out sample1.bam sample2.bam

# some chromosomes / one chromosome          -> out/<chrom>/<label>_CG.tsv
readontbam_fast -f ref.fa -c chr1,chr2 -x CG,CHG,CHH -m C+m. -p 0.7 -t 16 -o out *.bam

# one window (1-based, inclusive)            -> out/<label>_CG.tsv
readontbam_fast -f ref.fa -c chr1 -s 1000001 -e 2000000 -x CG -m C+m. -p 0.7 -o out s1=sample1.bam
```

| option | meaning (`readONTbam()` argument) |
|---|---|
| `-f` | reference FASTA with an up-to-date `.fai`; the same sequence as the BSgenome you would give `readONTbam()` |
| `-c` | `chr`: chromosome(s), comma separated. Omit or `all` = every sequence in the `.fai` (chromosomes missing from a BAM header are skipped) |
| `-s` / `-e` | optional window on one chromosome: keeps sites with `start <= pos <= end` |
| `-F` | reads to decode (one chromosome only). Default: the whole chromosome, or `chrom:(start+1)-end` for a window. `-F CHROM` with a window decodes every read of a pre-sliced BAM |
| `-x` | `context`: CG, CHG, CHH (comma list = one output per context) |
| `-m` | `modif`: MM code(s) counted as the modification, e.g. `C+m.` or `C+m?` (default `C+m?`) |
| `-p` | `prob_thresh`: a call is modified if `ML/255 >= p` (default 0.5) |
| `--synonymous` | `synonymous = TRUE` (CHG/CHH) |
| `--implicit-unmod` | **not** `readONTbam()` behaviour; see [Implicit calls](#implicit--mode-calls) |
| `-t` / `-j` | total threads / BAMs decoded concurrently |
| `-o`, `--overwrite` | output directory; existing outputs are skipped unless `--overwrite` |
| `--check-literal` | testing: re-run a literal transcription of the R position arithmetic for every call and abort on any difference |

BAMs must be coordinate-sorted and indexed. `LABEL=path.bam` sets the sample label;
by default it is the file name up to the first `.bam`. Output files are byte-identical
across thread counts and platforms.

### Output

One TSV per BAM, chromosome and context, with one row per reference cytosine (also
zero-coverage ones), sorted by position:

| column | content |
|---|---|
| `pos`, `strand` | position (1-based) and strand of the cytosine |
| `trinuc` | trinucleotide context (CHG/CHH; `NA` for CG) |
| `readsM`, `readsN` | reads calling the C modified / reads covering it with a call |
| `ONT_Cm`, `ONT_C` | comma-joined indices of the reads calling it modified / unmodified |
| `tri_named` | reproduces a naming quirk of `selectCytosine()` (used by the R helper only) |

Read indices number the decoded reads per chromosome (per window in window mode),
exactly as `readONTbam(chr = ...)` does.

### R

```r
source("R/readONTbam_fast.R")

## run the decoder and get GRanges back (one per BAM)
gr <- readONTbamFast(c(s1 = "sample1.bam", s2 = "sample2.bam"),
                     ref_fasta = "ref.fa", chrom = "chr1",     # add start=, end= for a window
                     context = "CG", modif = "C+m.", prob_thresh = 0.7, threads = 8,
                     genome = NULL,                             # or your BSgenome object
                     binary = "bin/readontbam_fast")

## or convert existing decoder output
gr1 <- readONTbamFastTsv("out/chr1/s1_CG.tsv", chrom = "chr1", context = "CG",
                         seqinfo = ontSeqinfo("chr1", fai = "ref.fa.fai"))

## one sample across several chromosomes of a whole-genome run, sorted like DMRcaller data
s1 <- loadONTSample("out", chroms = c("chr1", "chr2"), label = "s1", context = "CG",
                    fai = "ref.fa.fai")
```

The GRanges have the columns `context`, `trinucleotide_context`, `ONT_Cm`, `ONT_C`
(CharacterLists), `readsM`, `readsN`, as from `readONTbam()`. Pass the BSgenome you
would give `readONTbam()` as `genome` to get `identical()` objects; without it the
Seqinfo comes from the `.fai` (`genome`/`isCircular` are then `NA`).

## Self-test

`test/selftest.R` runs `DMRcaller::readONTbam()` and `readontbam_fast` on the same small
dataset and compares the results with `identical()`. It needs a DMRcaller version that
has `readONTbam()`, plus BSgenome (or BSgenomeForge) and rtracklayer. It forges and
installs a temporary BSgenome from the FASTA.

```bash
# a few dozen reads and the whole chromosome they map to, cut from your own data
test/make_test_data.sh sample1.bam ref.fa chr1:1000000-1010000 mytest
Rscript test/selftest.R CG,CHG mytest.bam mytest.fa
#   ...
#   CG   sites  ...  identical: TRUE
#   SELF-TEST PASSED
```

Keep the test small: `readONTbam()` takes about 3 s per read and context.
`SELFTEST_MODIF` (default `C+m.`) and `SELFTEST_PROB` (default 0.7) set `modif` and
`prob_thresh`.

## Validation

During development the decoder was compared with DMRcaller's `readONTbam()` /
`selectCytosine()` (source from the Bioconductor git, the version with the `cores`
argument, run in R 4.3; the decoding logic is unchanged in later releases) on Dorado
`5mC_5hmC` R10.4.1 data. All comparisons gave `identical()` results:

* 12 decoded datasets: CG, CHG and CHH; both strands, indels, soft clips; reads running
  past window ends; read numbering with `-F`; random read subsets spread over 1 Mb;
  `synonymous = TRUE`; window and whole-chromosome site sets.
* 12 cytosine enumerations, including motifs at chromosome ends, lowercase bases and Ns.
* The R wrapper end to end, for a window and a whole chromosome.
* `--check-literal` on 16 BAMs (~90k primary reads): no difference between the fast position
  arithmetic and a literal transcription of the R code.

## Behaviour reproduced from readONTbam

These are `readONTbam()`/`selectCytosine()` behaviours, kept so the output matches:

1. **Read filtering.** Primary alignments are used; secondary/supplementary ones only
   if they carry ML, have no hard clip and have SEQ.
2. **Reverse-strand Cs directly after a deletion** are shifted 1 bp onto the deleted
   base and then dropped by the reference-mismatch filter. Forward-strand insertion
   checks use `[pos, pos+len)`, reverse-strand ones `[pos, pos+len]`.
3. **Leading hard clips** are subtracted as if they were in SEQ (only matters for
   primary reads with a leading `H`).
4. **CHG minus-strand sites** are derived only from forward `CAG`/`CTG`, so minus-strand
   `CCG` sites (forward `CGG`) never appear.
5. **Thresholds.** Only calls whose code exactly matches `modif` count; a call is
   modified if `ML/255 >= prob_thresh` (0.7 → ML ≥ 179).
6. **Read indices** in `ONT_Cm`/`ONT_C` are per decoded set: per chromosome in
   whole-chromosome runs, per window in window runs. A read spanning two windows gets
   different indices in each, so read-level analyses can't be joined across windows;
   whole-chromosome runs avoid this. `readsM`/`readsN` are unaffected.

### Implicit (`.`) mode calls

In `.` (implicit) mode, e.g. `C+m.` as written by Dorado, Cs that the MM tag skips are
implicitly *unmodified*. `readONTbam()` only counts the Cs listed explicitly, so
coverage (`readsN`) is undercounted and methylation (`readsM/readsN`) overestimated. In
one CpG test set only ~73% of read Cs carried an explicit call, and mean methylation
was 0.57 with `readONTbam()` semantics versus 0.38 when skipped Cs count as unmodified
(as in the SAM specification and modkit).

`--implicit-unmod` (R: `implicit_unmod = TRUE`) counts skipped Cs as unmodified calls.
Its output then deliberately **differs** from `readONTbam()`; the default is off.

## Caveats

* **Stale FASTA index.** The decoder reads the reference through its `.fai`. If the
  FASTA was rewritten (e.g. headers shortened) after indexing, sequences are read from
  wrong offsets without any error. Re-run `samtools faidx ref.fa` when in doubt.
* The FASTA must hold the same sequence as the BSgenome used with DMRcaller, or the
  outputs are not comparable.

## Citation

If you use this tool, please also cite DMRcaller, whose methods it re-implements:

Catoni M, Tsang JMF, Greco AP, Zabet NR (2018). DMRcaller: a versatile R/Bioconductor
package for detection and visualization of differentially methylated regions in CpG and
non-CpG contexts. *Nucleic Acids Research* 46(19):e114.

## License

GPL-3 (see `LICENSE`), the same licence as DMRcaller, whose `readONTbam()` and
`selectCytosine()` logic this program is a translation of.
