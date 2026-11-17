#!/usr/bin/env Rscript
# readontbam_fast self-test: runs DMRcaller::readONTbam() (the DMRcaller installed in
# this R) and readontbam_fast on the same small dataset and checks that the two GRanges
# are identical().
#
# Usage:
#   Rscript test/selftest.R [CONTEXTS] [BAM FASTA]
#     CONTEXTS  CG (default), or a comma list such as CG,CHG,CHH
#     BAM       small coordinate-sorted, indexed BAM with MM/ML tags (default test/test.bam)
#     FASTA     reference of ONE chromosome, the one the reads map to (default test/test.fa)
#   Options via environment variables: SELFTEST_MODIF (default C+m.),
#   SELFTEST_PROB (default 0.7), READONTBAM_FAST (path to the binary).
#
# readONTbam() needs about 3 s per read and context, so keep the BAM to a few dozen
# reads (see test/make_test_data.sh). The whole chromosome in FASTA is enumerated, as
# readONTbam(chr = ...) does. Exit status 0 = PASS.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome)
})

args     <- commandArgs(trailingOnly = TRUE)
this_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
                    if (length(f)) dirname(normalizePath(f)) else getwd() })
tool_dir <- dirname(this_dir)
contexts <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else "CG"
bam      <- if (length(args) >= 3) normalizePath(args[2]) else file.path(this_dir, "test.bam")
fa_in    <- if (length(args) >= 3) normalizePath(args[3]) else file.path(this_dir, "test.fa")
modif    <- Sys.getenv("SELFTEST_MODIF", "C+m.")
prob     <- as.numeric(Sys.getenv("SELFTEST_PROB", "0.7"))
if (!file.exists(bam) || !file.exists(fa_in))
  stop("test data not found (", bam, ", ", fa_in, "). Make some with test/make_test_data.sh ",
       "and run: Rscript test/selftest.R CG my_test.bam my_test.fa")
source(file.path(tool_dir, "R", "readONTbam_fast.R"))

## DMRcaller: the installed package (SELFTEST_DMRCALLER_SRC = ";"-separated R files to source
## instead — only used to test against a DMRcaller source checkout)
src <- Sys.getenv("SELFTEST_DMRCALLER_SRC", "")
if (nzchar(src)) {
  suppressPackageStartupMessages({ library(Rsamtools); library(GenomicAlignments); library(BiocParallel) })
  for (f in strsplit(src, ";")[[1]]) source(f)
  dmr_version <- paste("source", src)
} else {
  suppressPackageStartupMessages(library(DMRcaller))
  if (!exists("readONTbam", asNamespace("DMRcaller")))
    stop("the installed DMRcaller has no readONTbam(); install a newer DMRcaller from Bioconductor")
  readONTbam <- DMRcaller::readONTbam
  selectCytosine <- DMRcaller::selectCytosine
  dmr_version <- as.character(packageVersion("DMRcaller"))
}

binary <- Sys.getenv("READONTBAM_FAST", "")
if (!nzchar(binary)) {
  binary <- file.path(tool_dir, "bin", "readontbam_fast")
  if (!file.exists(binary)) binary <- unname(Sys.which("readontbam_fast"))
}
stopifnot("readontbam_fast binary not found (run make)" = nzchar(binary) && file.exists(binary))

work <- tempfile("readontbam_selftest_")
dir.create(file.path(work, "lib"), recursive = TRUE)
dir.create(file.path(work, "seqs"))
.libPaths(c(file.path(work, "lib"), .libPaths()))

## ── reference: one chromosome; IUPAC codes other than N -> N (2bit cannot store them),
##    the same cleaned copy is given to both implementations ──
seqs <- readDNAStringSet(fa_in)
if (length(seqs) != 1) stop("FASTA must hold exactly one sequence (the chromosome of the reads)")
chrom <- sub(" .*", "", names(seqs)); names(seqs) <- chrom
seqs <- DNAStringSet(gsub("[^ACGTNacgtn]", "N", as.character(seqs)))
names(seqs) <- chrom
fa <- file.path(work, "ref.fa")
writeXStringSet(seqs, fa)
Rsamtools::indexFa(fa)

## ── forge + install a tiny BSgenome (readONTbam needs a BSgenome package) ──
pkg <- "BSgenome.Readontbam.selftest.v1"
rtracklayer::export(seqs, file.path(work, "seqs", "ref.2bit"))
writeLines(c(
  paste0("Package: ", pkg),
  "Title: readontbam_fast self-test genome",
  "Description: One reference chromosome used by the readontbam_fast self-test.",
  "Version: 0.0.1", "organism: Unknown", "common_name: Unknown", "provider: selftest",
  "genome: selftest", "release_date: 2026", "source_url: local", "organism_biocview: Unknown",
  "BSgenomeObjname: Selftest", paste0('seqnames: "', chrom, '"'), "circ_seqs: character(0)",
  "SrcDataFiles: local", "PkgExamples: genome",
  paste0("seqs_srcdir: ", file.path(work, "seqs")), "seqfile_name: ref.2bit"),
  file.path(work, "seed.dcf"))
forge <- if (requireNamespace("BSgenomeForge", quietly = TRUE)) BSgenomeForge::forgeBSgenomeDataPkg else
           BSgenome::forgeBSgenomeDataPkg
invisible(capture.output(suppressMessages(forge(file.path(work, "seed.dcf"), destdir = work))))
st <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load", "-l",
              shQuote(file.path(work, "lib")), shQuote(file.path(work, pkg))), stdout = FALSE, stderr = FALSE)
if (st != 0) stop("could not install the self-test BSgenome")
suppressPackageStartupMessages(library(pkg, character.only = TRUE, lib.loc = file.path(work, "lib")))
genome <- get(pkg)

## readONTbam() refuses BSgenome packages that are not in the online Bioconductor list when
## the internet is reachable; make it use the installed-packages fallback it already has.
invisible(try(utils::assignInNamespace("available.genomes", function(...) stop("offline"),
                                       ns = "BSgenome"), silent = TRUE))

cat(sprintf("readontbam_fast self-test | DMRcaller %s | binary %s\n  %s on %s (%s bp), modif %s, prob %s\n",
            dmr_version, binary, basename(bam), chrom, format(width(seqs), big.mark = ","), modif, prob))
ok <- TRUE
for (ctx in contexts) {
  t0 <- proc.time()
  invisible(capture.output(ref_gr <- suppressMessages(selectCytosine(genome = genome, context = ctx, chr = chrom))))
  seqinfo(ref_gr) <- seqinfo(genome)[seqlevels(ref_gr)]
  r_out <- suppressWarnings(suppressMessages({
    capture.output(res <- readONTbam(bamfile = bam, ref_gr = ref_gr, genome = genome, context = ctx,
                                     chr = chrom, modif = modif, region = NULL, prob_thresh = prob,
                                     parallel = FALSE))
    res }))
  t_r <- (proc.time() - t0)[["elapsed"]]
  t0 <- proc.time()
  cpp <- readONTbamFast(c(test = bam), ref_fasta = fa, chrom = chrom, context = ctx,
                        modif = modif, prob_thresh = prob, genome = genome, binary = binary)$test
  t_c <- (proc.time() - t0)[["elapsed"]]
  same <- identical(r_out, cpp)
  ok <- ok && same
  cat(sprintf("  %-3s  sites %8d  covered %6d  readsN %7d | R %7.1fs  C++ %5.1fs | identical: %s\n",
              ctx, length(r_out), sum(r_out$readsN > 0), sum(r_out$readsN), t_r, t_c, same))
  if (!same) print(all.equal(r_out, cpp))
}
unlink(work, recursive = TRUE)
cat(if (ok) "SELF-TEST PASSED\n" else "SELF-TEST FAILED\n")
quit(save = "no", status = if (ok) 0 else 1)
