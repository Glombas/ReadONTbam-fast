## readONTbam_fast.R — R side of readontbam_fast, a C++ port of DMRcaller::readONTbam()
##
##   source("R/readONTbam_fast.R")
##
##   * readONTbamFast()     : run the binary on BAM(s) and return a named list of GRanges
##   * readONTbamFastTsv()  : one readontbam_fast TSV  -> GRanges identical to readONTbam()
##   * loadONTSample()      : one sample over several chromosomes from an output directory
##   * ontSeqinfo()         : Seqinfo from a BSgenome or a FASTA .fai
##
## Copyright (C) 2026 Marek Glombik. GPL-3 (see LICENSE).
##
## The GRanges is assembled with the same R operations readONTbam()/selectCytosine()
## use (same column order/types, CompressedCharacterList ONT_Cm/ONT_C holding one
## comma-joined string per covered site, integer readsM/readsN), so identical() holds
## against the original output when the same Seqinfo is supplied.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

.ctxLevels <- c("CG", "CHG", "CHH")

## Seqinfo for the output. Pass the BSgenome object you would give readONTbam() to get
## objects identical to its output;
## otherwise it is built from the FASTA .fai (isCircular/genome then NA).
ontSeqinfo <- function(chrom, genome = NULL, fai = NULL) {
  if (!is.null(genome)) return(seqinfo(genome)[chrom])
  if (is.null(fai)) stop("give either a BSgenome (genome=) or the FASTA .fai path (fai=)")
  f <- utils::read.delim(fai, header = FALSE, stringsAsFactors = FALSE,
                         quote = "", comment.char = "")[, 1:2]
  names(f) <- c("name", "len")
  if (!chrom %in% f$name) stop(chrom, " not found in ", fai)
  Seqinfo(seqnames = chrom, seqlengths = as.integer(f$len[f$name == chrom]))
}

## Convert one readontbam_fast TSV into the GRanges readONTbam() would return.
readONTbamFastTsv <- function(tsv, chrom, context, seqinfo) {
  cls <- c("integer", "character", "character", "integer", "integer",
           "character", "character", "integer")
  d <- if (requireNamespace("data.table", quietly = TRUE)) {
    ## ~10x faster for whole-chromosome files; same columns/types as read.delim
    as.list(data.table::fread(tsv, sep = "\t", colClasses = cls, na.strings = "NA",
                              quote = "", showProgress = FALSE, nThread = 1))
  } else {
    utils::read.delim(tsv, colClasses = cls, na.strings = "NA", quote = "", comment.char = "")
  }
  n <- length(d$pos)

  ## --- reference cytosines, as built by selectCytosine() + window subset ---
  ref_gr <- GRanges(chrom, IRanges(start = d$pos, width = 1L), strand = d$strand)
  tri <- if (context == "CG") rep(NA_character_, n) else
           ## selectCytosine() leaves an NA-filled names attribute on CHG/CHH contexts
           ## (only the chromosome's first-matched trinucleotides are named <chrom>)
           stats::setNames(d$trinuc, ifelse(d$tri_named == 1L, chrom, NA_character_))
  mcols(ref_gr) <- DataFrame(context = factor(rep(context, n), levels = .ctxLevels),
                             trinucleotide_context = tri)
  seqinfo(ref_gr) <- seqinfo

  ## --- annotation, same steps as readONTbam() ---
  mcols(ref_gr)$ONT_Cm <- CharacterList(vector("list", n))
  mcols(ref_gr)$ONT_C  <- CharacterList(vector("list", n))
  ref_gr$readsM <- integer(n)
  ref_gr$readsN <- integer(n)

  idx_mod <- which(nzchar(d$ONT_Cm))
  mod_char_list <- stats::setNames(as.list(d$ONT_Cm[idx_mod]), idx_mod)
  ref_gr$ONT_Cm[idx_mod] <- CharacterList(mod_char_list)
  ref_gr$readsM[idx_mod] <- stats::setNames(d$readsM[idx_mod], idx_mod)

  idx_unmod <- which(nzchar(d$ONT_C))
  unmod_char_list <- stats::setNames(as.list(d$ONT_C[idx_unmod]), idx_unmod)
  idx_cov <- which(d$readsN > 0L)
  ref_gr$ONT_C[idx_unmod] <- CharacterList(unmod_char_list)
  ref_gr$readsN[idx_cov] <- stats::setNames(d$readsN[idx_cov], idx_cov)

  ref_gr
}

## Run readontbam_fast on one or more BAMs and return a named list of GRanges
## (one per BAM for a single context; a nested list [[bam]][[context]] otherwise).
##   bamfiles : named character vector (names = labels); unnamed -> labels from basenames
##   chrom    : one chromosome
##   start/end: sites with start <= pos <= end (same as start(gr) >= start & <= end);
##              leave both NULL for the whole chromosome
##   fetch    : reads to decode; default the whole chromosome, or chrom:(start+1)-end for a
##              window (== your samtools slicing)
##   implicit_unmod: TRUE deviates from readONTbam — Cs skipped in '.'-mode MM entries
##              (e.g. C+m.) are counted as unmodified reads (SAM-spec / modkit semantics)
readONTbamFast <- function(bamfiles, ref_fasta, chrom, start = NULL, end = NULL,
                           context = "CG", modif = "C+m?", prob_thresh = 0.5,
                           synonymous = FALSE, threads = 1L, fetch = NULL,
                           implicit_unmod = FALSE,
                           genome = NULL, binary = Sys.getenv("READONTBAM_FAST", "readontbam_fast"),
                           tsv_dir = tempfile("readontbam_fast_"), keep_tsv = FALSE) {
  if (is.null(names(bamfiles)) || any(!nzchar(names(bamfiles))))
    names(bamfiles) <- sub("\\.bam.*", "", basename(bamfiles))
  dir.create(tsv_dir, showWarnings = FALSE, recursive = TRUE)
  if (length(chrom) != 1) stop("readONTbamFast() handles one chromosome per call")
  if (is.null(start) != is.null(end)) stop("give both start and end, or neither")
  args <- c("-f", ref_fasta, "-c", chrom,
            if (!is.null(start)) c("-s", format(start, scientific = FALSE),
                                   "-e", format(end, scientific = FALSE)),
            "-x", paste(context, collapse = ","), "-m", paste(modif, collapse = ","),
            "-p", format(prob_thresh, digits = 17), "-t", threads, "-o", tsv_dir, "--overwrite",
            if (synonymous) "--synonymous",
            if (implicit_unmod) "--implicit-unmod",
            if (!is.null(fetch)) c("-F", fetch),
            paste0(names(bamfiles), "=", bamfiles))
  st <- system2(binary, shQuote(args))
  if (!identical(st, 0L)) stop("readontbam_fast failed (exit ", st, ")")
  si <- ontSeqinfo(chrom, genome = genome, fai = paste0(ref_fasta, ".fai"))
  out <- lapply(names(bamfiles), function(lab) {
    r <- lapply(context, function(ctx)
      readONTbamFastTsv(file.path(tsv_dir, if (is.null(start)) chrom else ".",
                                  paste0(lab, "_", ctx, ".tsv")), chrom, ctx, si))
    names(r) <- context
    if (length(context) == 1) r[[1]] else r
  })
  names(out) <- names(bamfiles)
  if (!keep_tsv) unlink(tsv_dir, recursive = TRUE)
  out
}

## ─────────────────────────────────────────────────────────────────────────────
## Loading whole-chromosome / whole-genome output
## ─────────────────────────────────────────────────────────────────────────────

## Load one sample (one context) for one or more chromosomes from a readontbam_fast
## output directory laid out as <dir>/<chrom>/<label>_<context>.tsv (decoder output for
## several chromosomes / the whole genome) or <dir>/<chrom>/<label>_<context>.rds
## (a GRanges you saved from readONTbamFastTsv()). Pieces are concatenated and sort()ed.
## Returns NULL (with a message) if any chromosome is missing.
loadONTSample <- function(input_dir, chroms, label, context, genome = NULL, fai = NULL) {
  pieces <- vector("list", length(chroms))
  for (k in seq_along(chroms)) {
    ch   <- chroms[k]
    base <- file.path(input_dir, ch, paste0(label, "_", context))
    if (file.exists(paste0(base, ".tsv"))) {
      pieces[[k]] <- readONTbamFastTsv(paste0(base, ".tsv"), ch, context,
                                       ontSeqinfo(ch, genome = genome, fai = fai))
    } else if (file.exists(paste0(base, ".rds"))) {
      pieces[[k]] <- readRDS(paste0(base, ".rds"))
    } else {
      message("MISSING input for ", label, " ", ch, " ", context, ": ", base, ".tsv/.rds")
      return(NULL)
    }
  }
  sort(do.call(c, pieces))
}
