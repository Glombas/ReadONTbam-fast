// readontbam_fast — C++/htslib re-implementation of DMRcaller::readONTbam()
//
// Reproduces the per-read MM/ML decoding, read filtering, genomic projection,
// insertion/diff/nocall filtering and per-cytosine aggregation of
// DMRcaller::readONTbam() (Bioconductor, versions with the `cores` argument;
// logic unchanged up to the cigarillo-based 1.45.x), *including its quirks*, so
// that the resulting GRanges are identical to the R output.
//
// It also enumerates the reference cytosines exactly like
// DMRcaller::selectCytosine() (followed by a start(gr) >= start & <= end
// window subset when a window is requested).
//
// Output: one TSV per (BAM, chromosome, context):
//   pos  strand  trinuc  readsM  readsN  ONT_Cm  ONT_C  tri_named
// with one row per reference cytosine (zero-coverage rows too), sorted by position.
// Works on a window, a whole chromosome, a list of chromosomes or the whole genome;
// per-site results are streamed to disk so memory does not grow with chromosome length. The companion R code (R/readONTbam_fast.R) turns this into
// the same GRanges object readONTbam() returns and saves it as .rds.
//
// Build:  make            (needs htslib; see README.md)
//
// Copyright (C) 2026 Marek Glombik
// Based on DMRcaller (Zabet, Catoni, Becker et al.; GPL-3), whose readONTbam()
// and selectCytosine() logic this program re-implements.
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version. It is distributed WITHOUT ANY WARRANTY; see the LICENSE file.

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <sys/stat.h>

static const char* VERSION = "0.2.0";

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------
struct Options {
    std::string ref_fa;
    std::vector<std::string> chroms; // empty = every sequence in the FASTA index
    int64_t win_start = -1;  // sites kept if win_start <= pos <= win_end (1-based); -1 = whole chromosome
    int64_t win_end = -1;
    std::string fetch;       // htslib region for reads (single chromosome only)
    std::vector<std::string> contexts{"CG"};
    std::vector<std::string> modif{"C+m?"};
    double prob_thresh = 0.50;
    bool synonymous = false;
    bool implicit_unmod = false;  // NOT readONTbam behaviour: count '.'-mode skipped Cs as unmodified
    int threads = 1;
    int jobs = 0;            // BAMs processed concurrently (0 = auto)
    std::string outdir = ".";
    bool check_literal = false;
    bool overwrite = false;
    std::vector<std::string> labels, bams;
};

static void usage() {
    fprintf(stderr,
"readontbam_fast %s — fast C++ port of DMRcaller::readONTbam()\n"
"\n"
"Usage:\n"
"  readontbam_fast -f ref.fa [-c CHROM[,CHROM...]] [-s START -e END] [options] [LABEL=]in.bam ...\n"
"\n"
"  whole genome      : readontbam_fast -f ref.fa ...            (every sequence in ref.fa.fai)\n"
"  some chromosomes  : readontbam_fast -f ref.fa -c chr1,chr2 ...\n"
"  one chromosome    : readontbam_fast -f ref.fa -c chr1 ...\n"
"  one window        : readontbam_fast -f ref.fa -c chr1 -s 1000001 -e 2000000 ...\n"
"\n"
"Required:\n"
"  -f, --ref FILE          reference FASTA with an up-to-date .fai (same sequence as the BSgenome\n"
"                          you would give readONTbam)\n"
"\n"
"Options:\n"
"  -c, --chrom LIST        chromosome(s), comma separated; default/'all' = every FASTA sequence\n"
"  -s, --start INT         window start (single chromosome only); sites with START <= pos <= END\n"
"  -e, --end INT           window end   (same semantics as start(gr) >= start & <= end).\n"
"                          Without -s/-e the whole chromosome is used.\n"
"  -F, --fetch REGION      reads to decode (single chromosome only). Default: the whole\n"
"                          chromosome, or CHROM:(START+1)-END when a window is given (= the\n"
"                          reads `samtools view in.bam CHROM:(START+1)-END` returns).\n"
"  -x, --context LIST      comma list of CG,CHG,CHH (default CG)\n"
"  -m, --modif LIST        MM codes counted as the modification (default C+m?)\n"
"  -p, --prob-thresh NUM   ML probability threshold, prob >= thr is modified (default 0.5)\n"
"      --synonymous        readONTbam(synonymous=TRUE) for CHG/CHH\n"
"      --implicit-unmod    DEVIATES FROM readONTbam: for MM entries in implicit mode\n"
"                          (e.g. C+m.) also count the skipped Cs as unmodified calls, as\n"
"                          the SAM spec / modkit do. Default off (exact readONTbam output).\n"
"  -t, --threads INT       total worker threads (default 1)\n"
"  -j, --jobs INT          BAMs decoded concurrently (default min(threads, #BAMs))\n"
"  -o, --outdir DIR        output directory (default .)\n"
"      --overwrite         overwrite existing outputs (default: skip outputs already present)\n"
"      --check-literal     also run a literal (slow) transcription of the R per-C\n"
"                          arithmetic and abort on any disagreement (testing)\n"
"\n"
"Output: whole chromosome(s) -> OUTDIR/CHROM/LABEL_CONTEXT.tsv; a window (-s/-e) -> OUTDIR/LABEL_CONTEXT.tsv.\n"
"LABEL defaults to the BAM basename with everything from the first '.bam' removed, like\n"
"gsub(\"\\\\.bam.*\", \"\", basename(bam)). Read indices (ONT_Cm/ONT_C) are numbered per\n"
"chromosome (or per window), exactly as readONTbam(chr = CHROM) numbers them.\n"
"Input BAMs must be coordinate-sorted and indexed. Results are streamed to disk, so memory\n"
"does not grow with chromosome length (apart from the reference and a site index).\n",
    VERSION);
}

static std::vector<std::string> split(const std::string& s, char d) {
    std::vector<std::string> out;
    size_t b = 0;
    while (true) {
        size_t e = s.find(d, b);
        out.push_back(s.substr(b, e == std::string::npos ? std::string::npos : e - b));
        if (e == std::string::npos) break;
        b = e + 1;
    }
    return out;
}

static std::string default_label(const std::string& path) {
    std::string base = path.substr(path.find_last_of('/') == std::string::npos ? 0 : path.find_last_of('/') + 1);
    size_t p = base.find(".bam");
    return p == std::string::npos ? base : base.substr(0, p);
}

static Options parse_args(int argc, char** argv) {
    Options o;
    auto need = [&](int& i) -> std::string {
        if (i + 1 >= argc) { fprintf(stderr, "missing value for %s\n", argv[i]); exit(2); }
        return std::string(argv[++i]);
    };
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(); exit(0); }
        else if (a == "--version") { printf("%s\n", VERSION); exit(0); }
        else if (a == "-f" || a == "--ref") o.ref_fa = need(i);
        else if (a == "-c" || a == "--chrom") { std::string v = need(i); if (v != "all") o.chroms = split(v, ','); }
        else if (a == "-s" || a == "--start") o.win_start = std::stoll(need(i));
        else if (a == "-e" || a == "--end") o.win_end = std::stoll(need(i));
        else if (a == "-F" || a == "--fetch") o.fetch = need(i);
        else if (a == "-x" || a == "--context") o.contexts = split(need(i), ',');
        else if (a == "-m" || a == "--modif") o.modif = split(need(i), ',');
        else if (a == "-p" || a == "--prob-thresh") o.prob_thresh = std::stod(need(i));
        else if (a == "--synonymous") o.synonymous = true;
        else if (a == "--implicit-unmod") o.implicit_unmod = true;
        else if (a == "-t" || a == "--threads") o.threads = std::stoi(need(i));
        else if (a == "-j" || a == "--jobs") o.jobs = std::stoi(need(i));
        else if (a == "-o" || a == "--outdir") o.outdir = need(i);
        else if (a == "--overwrite") o.overwrite = true;
        else if (a == "--check-literal") o.check_literal = true;
        else if (!a.empty() && a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); usage(); exit(2); }
        else {
            size_t eq = a.find('=');
            if (eq != std::string::npos && eq > 0) { o.labels.push_back(a.substr(0, eq)); o.bams.push_back(a.substr(eq + 1)); }
            else { o.labels.push_back(default_label(a)); o.bams.push_back(a); }
        }
    }
    if (o.ref_fa.empty() || o.bams.empty()) { usage(); exit(2); }
    const bool has_s = o.win_start >= 0, has_e = o.win_end >= 0;
    if (has_s != has_e) { fprintf(stderr, "give both -s and -e (or neither for whole chromosomes)\n"); exit(2); }
    if (has_s && (o.chroms.size() != 1 || o.win_end < o.win_start)) {
        fprintf(stderr, "-s/-e need exactly one chromosome (-c) and START <= END\n"); exit(2);
    }
    if (!o.fetch.empty() && o.chroms.size() != 1) { fprintf(stderr, "-F needs exactly one chromosome (-c)\n"); exit(2); }
    for (auto& c : o.contexts)
        if (c != "CG" && c != "CHG" && c != "CHH") { fprintf(stderr, "unknown context %s\n", c.c_str()); exit(2); }
    for (size_t i = 0; i < o.contexts.size(); ++i)
        for (size_t j = 0; j < i; ++j)
            if (o.contexts[i] == o.contexts[j]) { fprintf(stderr, "context %s given twice\n", o.contexts[i].c_str()); exit(2); }
    if (!(o.prob_thresh > 0 && o.prob_thresh < 1)) { fprintf(stderr, "prob-thresh must be in (0,1)\n"); exit(2); }
    if (o.synonymous) {
        for (auto& c : o.contexts)
            if (c == "CG") { fprintf(stderr, "--synonymous requires CHG or CHH context (as in readONTbam)\n"); exit(2); }
    }
    if (o.threads < 1) o.threads = 1;
    return o;
}

// ---------------------------------------------------------------------------
// Reference + cytosine enumeration (DMRcaller::selectCytosine semantics)
// ---------------------------------------------------------------------------
struct Sites {                       // one context on one chromosome window
    std::string ctx;
    int kind = 0;                    // 0=CG 1=CHG 2=CHH
    std::vector<int32_t> pos;        // 1-based, ascending
    std::vector<char> strand;        // '+' / '-'
    // selectCytosine() quirk, reproduced for identical(): getSeq(genome, chr, pos, width=3)
    // puts the (length-1) name only on its FIRST element, so per chromosome exactly one
    // forward trinucleotide (first match of the first CHG/CHH pattern that occurs) and one
    // reverse one carry names = <chrom>; every other element has an NA name.
    int64_t named_plus = -1, named_minus = -1;
};

struct ChromSites {
    std::string chrom;
    int64_t ws = 1, we = 0;          // window (1-based, inclusive)
    std::vector<Sites> ctx;          // one per requested context
    // at[p - ws] = (site index << 2) | (context slot + 1); 0 = no site. A position is a
    // site of at most one context (C followed by G / H+G / H+H are mutually exclusive).
    std::vector<uint32_t> at;
    inline bool lookup(int64_t gp, char strand, int& c, int32_t& si) const {
        if (gp < ws || gp > we) return false;
        const uint32_t v = at[gp - ws];
        if (!v) return false;
        c = (int)(v & 3u) - 1; si = (int32_t)(v >> 2);
        return ctx[c].strand[si] == strand;
    }
};

static inline bool isH_fwd(char c) { return c == 'A' || c == 'C' || c == 'T'; }
static inline bool isH_rev(char c) { return c == 'A' || c == 'G' || c == 'T'; }
static inline char comp(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        case 'M': return 'K'; case 'K': return 'M'; case 'R': return 'Y'; case 'Y': return 'R';
        case 'S': return 'S'; case 'W': return 'W'; case 'B': return 'V'; case 'V': return 'B';
        case 'D': return 'H'; case 'H': return 'D'; case 'N': return 'N';
        default: return c;
    }
}

// ref: uppercase chromosome string, ref[i] = base at 1-based position i+1.
// Same sites as selectCytosine(genome, context, chr) restricted to ws <= pos <= we:
//   CG : + at C of "CG", - at its G
//   CHG: + at C of CAG/CCG/CTG, - at the G of forward CAG/CTG only (selectCytosine derives
//        reverse CHG sites from forward hits, so minus-strand CCG = forward CGG is missing)
//   CHH: + at C of C[ACT][ACT], - at the G of [AGT][AGT]G
static Sites enumerate_sites(const std::string& ref, const std::string& ctx, int64_t ws, int64_t we) {
    Sites S; S.ctx = ctx; S.kind = ctx == "CG" ? 0 : ctx == "CHG" ? 1 : 2;
    const int64_t L = (int64_t)ref.size();
    auto first_match = [&](const char* pat) -> int64_t {
        size_t f = ref.find(pat); return f == std::string::npos ? -1 : (int64_t)f + 1; };
    if (S.kind == 1) {
        for (const char* pat : {"CAG", "CCG", "CTG"}) {
            int64_t q = first_match(pat);
            if (q > 0) { S.named_plus = q; if (pat[1] != 'C') S.named_minus = q + 2; break; }
        }
    } else if (S.kind == 2) {
        for (const char* pat : {"CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"}) {
            int64_t q = first_match(pat); if (q > 0) { S.named_plus = q; break; }
        }
        for (const char* pat : {"TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"}) {
            int64_t q = first_match(pat); if (q > 0) { S.named_minus = q + 2; break; }
        }
    }
    auto B = [&](int64_t q) -> char { return (q >= 1 && q <= L) ? ref[q - 1] : 'N'; };
    for (int64_t q = std::max<int64_t>(1, ws); q <= std::min<int64_t>(L, we); ++q) {
        const char b = ref[q - 1];
        bool plus = false, minus = false;
        if (b == 'C') {
            if (S.kind == 0)      plus = q + 1 <= L && B(q + 1) == 'G';
            else if (S.kind == 1) plus = q + 2 <= L && isH_fwd(B(q + 1)) && B(q + 2) == 'G';
            else                  plus = q + 2 <= L && isH_fwd(B(q + 1)) && isH_fwd(B(q + 2));
        } else if (b == 'G') {
            if (S.kind == 0)      minus = q - 1 >= 1 && B(q - 1) == 'C';
            else if (S.kind == 1) minus = q - 2 >= 1 && B(q - 2) == 'C' && (B(q - 1) == 'A' || B(q - 1) == 'T');
            else                  minus = q - 2 >= 1 && isH_rev(B(q - 1)) && isH_rev(B(q - 2));
        }
        if (plus || minus) { S.pos.push_back((int32_t)q); S.strand.push_back(plus ? '+' : '-'); }
    }
    return S;
}

static ChromSites build_sites(const std::string& chrom, const std::string& ref,
                              const std::vector<std::string>& contexts, int64_t ws, int64_t we) {
    ChromSites CS; CS.chrom = chrom; CS.ws = ws; CS.we = we;
    CS.at.assign((size_t)(we - ws + 1), 0u);
    for (size_t c = 0; c < contexts.size(); ++c) {
        CS.ctx.push_back(enumerate_sites(ref, contexts[c], ws, we));
        const Sites& S = CS.ctx.back();
        if (S.pos.size() >= (1u << 30)) throw std::runtime_error("too many sites for the site index");
        for (size_t i = 0; i < S.pos.size(); ++i) {
            uint32_t& slot = CS.at[S.pos[i] - ws];
            if (slot) throw std::runtime_error("position is a site in two contexts (unexpected)");
            slot = ((uint32_t)i << 2) | (uint32_t)(c + 1);
        }
    }
    return CS;
}

// ---------------------------------------------------------------------------
// Per-read decoding (DMRcaller:::.computeModifPos semantics)
// ---------------------------------------------------------------------------
struct Call { int32_t ctx; int32_t site; bool mod; };  // mod: prob >= thresh

struct Decoder {
    const Options* opt;
    const std::string* ref;          // uppercase chromosome
    const ChromSites* cs;            // reference cytosines of this chromosome/window

    // scratch (per thread)
    std::vector<char> seq, aln;
    std::vector<int64_t> ins_pos, ins_len, ins_pref, del_pos, del_len, del_pref;
    std::vector<char> ins_align;     // index = ref offset (0..ref_len), 1 if in ins_align_pos
    std::vector<int32_t> cyt;
    std::vector<uint8_t> ml;
    int64_t literal_checks = 0;

    // Appends this read's reference-cytosine calls to out.
    void decode(const bam1_t* b, std::vector<Call>& out);
};

static inline bool code_in(const std::vector<std::string>& v, const char* s, size_t n) {
    for (auto& x : v) if (x.size() == n && memcmp(x.data(), s, n) == 0) return true;
    return false;
}

void Decoder::decode(const bam1_t* b, std::vector<Call>& out) {
    const uint32_t* cig = bam_get_cigar(b);
    const uint32_t ncig = b->core.n_cigar;
    const bool rev = (b->core.flag & BAM_FREVERSE) != 0;
    const int64_t start = (int64_t)b->core.pos + 1;

    uint8_t* mmp = bam_aux_get(b, "MM");
    uint8_t* mlp = bam_aux_get(b, "ML");
    if (!mmp || !mlp) return;                           // readONTbam: no MM/ML -> no calls
    const char* mm = bam_aux2Z(mmp);
    if (!mm) return;
    if (mlp[0] != 'B') return;
    uint32_t nml = bam_auxB_len(mlp);
    if (nml == 0) return;
    ml.resize(nml);
    for (uint32_t i = 0; i < nml; ++i) ml[i] = (uint8_t)bam_auxB2i(mlp, i);

    // --- CIGAR-derived quantities ---
    int64_t lead_clip = 0;                               // .getClipping(): first op S or H
    if (ncig > 0) {
        int op0 = bam_cigar_op(cig[0]);
        if (op0 == BAM_CSOFT_CLIP || op0 == BAM_CHARD_CLIP) lead_clip = bam_cigar_oplen(cig[0]);
    }
    int64_t read_len = 0, ref_len = 0;                   // extents along query / reference
    ins_pos.clear(); ins_len.clear(); del_pos.clear(); del_len.clear();
    for (uint32_t k = 0; k < ncig; ++k) {
        int op = bam_cigar_op(cig[k]); int64_t len = bam_cigar_oplen(cig[k]);
        switch (op) {
            case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF: read_len += len; ref_len += len; break;
            case BAM_CINS: ins_pos.push_back(read_len + 1); ins_len.push_back(len); read_len += len; break;
            case BAM_CSOFT_CLIP: read_len += len; break;
            case BAM_CDEL: del_pos.push_back(read_len + 1); del_len.push_back(len); ref_len += len; break;
            case BAM_CREF_SKIP: ref_len += len; break;
            default: break;                               // H, P
        }
    }
    ins_pref.assign(ins_len.size() + 1, 0);
    for (size_t k = 0; k < ins_len.size(); ++k) ins_pref[k + 1] = ins_pref[k] + ins_len[k];
    del_pref.assign(del_len.size() + 1, 0);
    for (size_t k = 0; k < del_len.size(); ++k) del_pref[k + 1] = del_pref[k] + del_len[k];

    // ins_align_pos = unique(ref offset (1-based) of the base following each I) - 1
    ins_align.assign((size_t)ref_len + 2, 0);
    {
        int64_t rp = 1;
        for (uint32_t k = 0; k < ncig; ++k) {
            int op = bam_cigar_op(cig[k]); int64_t len = bam_cigar_oplen(cig[k]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) rp += len;
            else if (op == BAM_CINS) { if (rp - 1 >= 0 && rp - 1 <= ref_len + 1) ins_align[rp - 1] = 1; }
            else if (op == BAM_CDEL || op == BAM_CREF_SKIP) rp += len;
        }
    }

    // read sequence (SEQ, reference orientation, incl. soft clips)
    const int32_t lq = b->core.l_qseq;
    const uint8_t* s4 = bam_get_seq(b);
    seq.resize(lq);
    for (int32_t i = 0; i < lq; ++i) seq[i] = seq_nt16_str[bam_seqi(s4, i)];

    // aligned_seq = sequenceLayer(seq, cigar, from="query", to="reference"): D -> '-', N -> '.'
    aln.assign((size_t)ref_len, '-');
    {
        int64_t q = 0, r = 0;
        for (uint32_t k = 0; k < ncig; ++k) {
            int op = bam_cigar_op(cig[k]); int64_t len = bam_cigar_oplen(cig[k]);
            switch (op) {
                case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                    for (int64_t j = 0; j < len; ++j) aln[r + j] = (q + j < lq) ? seq[q + j] : 'N';
                    q += len; r += len; break;
                case BAM_CINS: case BAM_CSOFT_CLIP: q += len; break;
                case BAM_CDEL: for (int64_t j = 0; j < len; ++j) aln[r + j] = '-'; r += len; break;
                case BAM_CREF_SKIP: for (int64_t j = 0; j < len; ++j) aln[r + j] = '.'; r += len; break;
                default: break;
            }
        }
    }

    // cytosine positions in the read's original orientation (1-based)
    cyt.clear();
    if (!rev) { for (int32_t i = 0; i < lq; ++i) if (seq[i] == 'C') cyt.push_back(i + 1); }
    else      { for (int32_t q = lq - 1; q >= 0; --q) if (seq[q] == 'G') cyt.push_back(lq - q); }
    const int64_t ncyt = (int64_t)cyt.size();

    const std::string& R = *ref;
    const int64_t chrL = (int64_t)R.size();
    auto refb = [&](int64_t gpos) -> char { return (gpos >= 1 && gpos <= chrL) ? R[gpos - 1] : 'N'; };

    // --- per-cytosine call: readONTbam's position arithmetic + insert/diff/nocall filters ---
    auto handle = [&](const int64_t rp, const bool is_mod) {
        int64_t gp;
        int insert = 0;
        if (!rev) {
            const int64_t x = rp;
            // insert: any(x >= ip & x < ip + il)
            size_t j = std::upper_bound(ins_pos.begin(), ins_pos.end(), x) - ins_pos.begin();
            if (j > 0 && x < ins_pos[j - 1] + ins_len[j - 1]) insert = 1;
            size_t ki = std::lower_bound(ins_pos.begin(), ins_pos.end(), x) - ins_pos.begin();
            int64_t insert_off = ins_pref[ki];
            size_t kd = std::lower_bound(del_pos.begin(), del_pos.end(), x) - del_pos.begin();
            int64_t delete_off = del_pref[kd];
            bool del_hit = kd < del_pos.size() && del_pos[kd] == x;
            gp = start - lead_clip + (x - insert_off + delete_off) - 1;
            if (del_hit) gp += del_len[kd];
        } else {
            const int64_t x = read_len - rp + 1;
            // insert: any(x >= ip & x <= ip + il)
            size_t j = std::upper_bound(ins_pos.begin(), ins_pos.end(), x) - ins_pos.begin();
            if (j > 0 && x <= ins_pos[j - 1] + ins_len[j - 1]) insert = 1;
            size_t ki = std::lower_bound(ins_pos.begin(), ins_pos.end(), x) - ins_pos.begin();
            int64_t insert_off = ins_pref[ki];
            size_t kd = std::lower_bound(del_pos.begin(), del_pos.end(), x) - del_pos.begin();
            int64_t delete_off = del_pref[kd];
            bool del_hit = kd < del_pos.size() && del_pos[kd] == x;
            gp = start - lead_clip + (x - insert_off + delete_off) - 1;
            if (del_hit) {
                bool within = false;   // any(x >= dp & x < dp + dl) over all deletions
                for (size_t d = 0; d < del_pos.size() && !within; ++d)
                    if (x >= del_pos[d] && x < del_pos[d] + del_len[d]) within = true;
                gp += del_len[kd] + (within ? -1 : 0);
            }
        }

        if (opt->check_literal) {
            // literal transcription of the R code (vector scans) for self-checking
            int64_t lx = rev ? (read_len - rp + 1) : rp;
            int lins = 0, ldel = 0; int64_t io = 0, dof = 0;
            for (size_t k = 0; k < ins_pos.size(); ++k) {
                if (!rev ? (lx >= ins_pos[k] && lx < ins_pos[k] + ins_len[k])
                         : (lx >= ins_pos[k] && lx <= ins_pos[k] + ins_len[k])) lins = 1;
                if (ins_pos[k] < lx) io += ins_len[k];
            }
            int64_t first_dl = -1;
            for (size_t k = 0; k < del_pos.size(); ++k) {
                if (del_pos[k] == lx) { ldel = 1; if (first_dl < 0) first_dl = del_len[k]; }
                if (del_pos[k] < lx) dof += del_len[k];
            }
            int64_t lgp = start - lead_clip + (lx - io + dof) - 1;
            if (ldel) {
                if (!rev) lgp += first_dl;
                else {
                    bool w = false;
                    for (size_t k = 0; k < del_pos.size(); ++k) if (lx >= del_pos[k] && lx < del_pos[k] + del_len[k]) w = true;
                    lgp += first_dl + (w ? -1 : 0);
                }
            }
            if (lgp != gp || lins != insert) {
                fprintf(stderr, "LITERAL MISMATCH read=%s rp=%lld gp=%lld/%lld ins=%d/%d\n",
                        bam_get_qname(b), (long long)rp, (long long)gp, (long long)lgp, insert, lins);
                abort();
            }
            ++literal_checks;
        }

        // bounds
        if (gp >= start + ref_len || gp < start) return;
        const int64_t ref_idx = gp - start + 1;       // 1-based within the read's ref span
        if (ref_idx < 1 || ref_idx > ref_len) return;

        // does gp hit a requested reference cytosine on this read's strand?
        int c; int32_t si;
        if (!cs->lookup(gp, rev ? '-' : '+', c, si)) return;   // call cannot contribute

        const char ref_base = refb(gp);
        const char aln_base = aln[ref_idx - 1];
        int diff;
        if (!rev) diff = (ref_base != aln_base) || (aln_base != 'C');
        else      diff = (ref_base != aln_base) || (aln_base != 'G');

        // insertion immediately adjacent (C followed by / G preceded by an insertion)
        if (!rev) { if (ins_align[ref_idx]) insert = 1; }
        else      { if (ins_align[ref_idx - 1]) insert = 1; }

        if (insert || diff) return;

        auto A = [&](int64_t i) -> char { return aln[i - 1]; };            // 1-based
        auto Rf = [&](int64_t i) -> char { return refb(start + i - 1); };  // 1-based

        {
            int nocall = 0;
            const int kind = cs->ctx[c].kind;
            if (kind == 0) {  // CG
                int64_t i2 = rev ? ref_idx - 1 : ref_idx + 1;
                if (i2 >= 1 && i2 <= ref_len) nocall = (Rf(i2) != A(i2));
            } else {
                const bool chg = (kind == 1);
                if (!rev) {
                    if (ref_idx + 2 <= ref_len) {
                        char rc = Rf(ref_idx), r1 = Rf(ref_idx + 1), r2 = Rf(ref_idx + 2);
                        char ac = A(ref_idx), a1 = A(ref_idx + 1), a2 = A(ref_idx + 2);
                        nocall = 1;
                        if (opt->synonymous) {
                            if (chg) { if (rc == 'C' && ac == 'C' && r2 == 'G' && a2 == 'G' && isH_fwd(r1) && isH_fwd(a1)) nocall = 0; }
                            else     { if (rc == 'C' && ac == 'C' && isH_fwd(r1) && isH_fwd(a1) && isH_fwd(r2) && isH_fwd(a2)) nocall = 0; }
                        } else if (rc == ac && r1 == a1 && r2 == a2) nocall = 0;
                    }
                } else {
                    if (ref_idx - 2 >= 1) {
                        char rg = Rf(ref_idx), r1 = Rf(ref_idx - 1), r2 = Rf(ref_idx - 2);
                        char ag = A(ref_idx), a1 = A(ref_idx - 1), a2 = A(ref_idx - 2);
                        nocall = 1;
                        if (opt->synonymous) {
                            if (chg) { if (rg == 'G' && ag == 'G' && r2 == 'C' && a2 == 'C' && isH_rev(r1) && isH_rev(a1)) nocall = 0; }
                            else     { if (rg == 'G' && ag == 'G' && isH_rev(r1) && isH_rev(a1) && isH_rev(r2) && isH_rev(a2)) nocall = 0; }
                        } else if (rg == ag && r1 == a1 && r2 == a2) nocall = 0;
                    }
                }
            }
            if (nocall) return;
            out.push_back(Call{(int32_t)c, si, is_mod});
        }
    };

    // --- walk MM entries ---
    uint32_t ml_idx = 0;
    const char* p = mm;
    while (*p) {
        const char* e = strchr(p, ';');
        const char* end = e ? e : p + strlen(p);
        if (end == p) { p = e ? e + 1 : end; continue; }
        // code = up to first ','  (e.g. "C+m." / "C+m?")
        const char* comma = (const char*)memchr(p, ',', end - p);
        const char* code_end = comma ? comma : end;
        const bool wanted = code_in(opt->modif, p, code_end - p);
        // --implicit-unmod: in '.' mode, Cs skipped by MM are implicit canonical calls
        const bool implicit = wanted && opt->implicit_unmod && code_end > p && code_end[-1] == '.';
        int64_t c_count = 0;
        const char* q = comma ? comma + 1 : end;
        while (comma && q <= end) {
            char* qe;
            long v = strtol(q, &qe, 10);
            if (qe == q) break;
            const int64_t prev = c_count;
            c_count += v + 1;
            uint32_t this_ml = ml_idx++;
            q = qe; if (q < end && *q == ',') ++q; else q = end + 1;

            if (!wanted) continue;
            if (implicit) for (int64_t k = prev + 1; k < c_count && k <= ncyt; ++k) handle(cyt[k - 1], false);
            if (c_count > ncyt) continue;                  // valid <- c_counts <= length(cyt_pos)
            if (this_ml >= nml) continue;                  // ML too short (R would give NA)
            handle(cyt[c_count - 1], ((double)ml[this_ml] / 255.0) >= opt->prob_thresh);
        }
        if (implicit) for (int64_t k = c_count + 1; k <= ncyt; ++k) handle(cyt[k - 1], false);
        p = e ? e + 1 : end;
    }
}

// ---------------------------------------------------------------------------
// Per-BAM processing: streamed per-site aggregation + TSV writer
// ---------------------------------------------------------------------------
static bool file_exists(const std::string& p) { struct stat st; return stat(p.c_str(), &st) == 0; }

static inline void put_int(std::string& s, int64_t v) {
    char buf[24]; int n = 0;
    if (v == 0) { s.push_back('0'); return; }
    bool neg = v < 0; uint64_t u = neg ? (uint64_t)(-v) : (uint64_t)v;
    while (u) { buf[n++] = (char)('0' + u % 10); u /= 10; }
    if (neg) s.push_back('-');
    while (n) s.push_back(buf[--n]);
}

// Sites are written in position order as soon as no later read can reach them. Reads come
// in coordinate order and every call lands at gp >= read start, so all sites before the
// start of the most recently fetched record are final.
struct StreamOut {
    struct Entry { std::vector<int32_t> mod, unmod; };  // read indices, ascending, unique
    const ChromSites* cs = nullptr; int c = 0; const std::string* ref = nullptr;
    std::string path, tmp;
    FILE* f = nullptr;
    std::deque<Entry> win;            // entries for sites [base, base + win.size())
    int64_t base = 0;                 // next site index to write
    std::string buf;

    void open(const std::string& p) {
        path = p; tmp = p + ".tmp";
        f = fopen(tmp.c_str(), "w");
        if (!f) throw std::runtime_error("cannot write " + tmp);
        buf = "pos\tstrand\ttrinuc\treadsM\treadsN\tONT_Cm\tONT_C\ttri_named\n";
    }
    inline void add(int32_t si, bool mod, int32_t rid) {
        if (si < base) throw std::runtime_error("call before an already written site — is the BAM coordinate-sorted?");
        const size_t off = (size_t)(si - base);
        if (off >= win.size()) win.resize(off + 1);
        std::vector<int32_t>& v = mod ? win[off].mod : win[off].unmod;
        if (v.empty() || v.back() != rid) v.push_back(rid);
    }
    void write_row(int64_t i, const Entry* e) {
        const Sites& S = cs->ctx[c];
        const int64_t p = S.pos[i]; const char st = S.strand[i];
        put_int(buf, p); buf.push_back('\t'); buf.push_back(st); buf.push_back('\t');
        if (S.kind == 0) buf += "NA";
        else if (st == '+') buf.append(ref->data() + (p - 1), 3);
        else { const char* q = ref->data() + (p - 3); buf.push_back(comp(q[2])); buf.push_back(comp(q[1])); buf.push_back(comp(q[0])); }
        buf.push_back('\t');
        size_t nm = 0, nn = 0;
        if (e) {
            nm = e->mod.size();
            // readsN = unique reads over modified + unmodified calls (both sorted ascending)
            size_t a = 0, b = 0;
            while (a < e->mod.size() || b < e->unmod.size()) {
                if (b == e->unmod.size() || (a < e->mod.size() && e->mod[a] < e->unmod[b])) ++a;
                else if (a == e->mod.size() || e->unmod[b] < e->mod[a]) ++b;
                else { ++a; ++b; }
                ++nn;
            }
        }
        put_int(buf, (int64_t)nm); buf.push_back('\t'); put_int(buf, (int64_t)nn); buf.push_back('\t');
        if (e) for (size_t k = 0; k < e->mod.size(); ++k) { if (k) buf.push_back(','); put_int(buf, e->mod[k]); }
        buf.push_back('\t');
        if (e) for (size_t k = 0; k < e->unmod.size(); ++k) { if (k) buf.push_back(','); put_int(buf, e->unmod[k]); }
        const bool named = (st == '+' && p == S.named_plus) || (st == '-' && p == S.named_minus);
        buf += named ? "\t1\n" : "\t0\n";
        if (buf.size() > (1u << 20)) drain();
    }
    void drain() {
        if (!buf.empty() && fwrite(buf.data(), 1, buf.size(), f) != buf.size()) throw std::runtime_error("write failed " + tmp);
        buf.clear();
    }
    // write every site with position < pos_limit
    void flush_before(int64_t pos_limit) {
        const std::vector<int32_t>& P = cs->ctx[c].pos;
        const int64_t lim = std::lower_bound(P.begin() + base, P.end(), pos_limit) - P.begin();
        while (base < lim) {
            if (!win.empty()) { write_row(base, &win.front()); win.pop_front(); }
            else write_row(base, nullptr);
            ++base;
        }
    }
    void finish() {
        flush_before(INT64_MAX);
        drain();
        if (fclose(f) != 0) { f = nullptr; throw std::runtime_error("write failed " + tmp); }
        f = nullptr;
        if (rename(tmp.c_str(), path.c_str()) != 0) throw std::runtime_error("rename failed " + path);
    }
    void abandon() { if (f) { fclose(f); f = nullptr; } remove(tmp.c_str()); }
};

struct BamResult { std::string label; bool ok = false; bool skipped = false; bool absent = false; std::string err;
                   double secs = 0; int64_t n_records = 0, n_kept = 0, n_calls = 0; };

static BamResult process_bam(const Options& o, const std::string& label, const std::string& bam,
                             const std::string& fetch, const std::vector<std::string>& outs,
                             const std::string& ref, const ChromSites& cs, int nthreads, bool absent_ok) {
    BamResult R; R.label = label;
    auto t0 = std::chrono::steady_clock::now();
    bool all_exist = true;
    for (auto& p : outs) if (!file_exists(p)) all_exist = false;
    if (all_exist && !o.overwrite) { R.ok = true; R.skipped = true; return R; }

    samFile* fp = sam_open(bam.c_str(), "r");
    if (!fp) { R.err = "cannot open " + bam; return R; }
    if (nthreads > 1) hts_set_threads(fp, std::min(nthreads, 8));
    sam_hdr_t* hdr = sam_hdr_read(fp);
    hts_idx_t* idx = hdr ? sam_index_load(fp, bam.c_str()) : nullptr;
    if (!hdr || !idx) { R.err = "cannot read header/index of " + bam; if (hdr) sam_hdr_destroy(hdr); sam_close(fp); return R; }
    const int tid = sam_hdr_name2tid(hdr, cs.chrom.c_str());
    if (tid < 0) {
        hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
        if (absent_ok) { R.ok = true; R.absent = true; return R; }
        R.err = "chromosome " + cs.chrom + " not in BAM header"; return R;
    }
    hts_itr_t* itr = sam_itr_querys(idx, hdr, fetch.c_str());
    if (!itr) { R.err = "bad fetch region " + fetch; hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp); return R; }

    std::vector<StreamOut> so(cs.ctx.size());
    std::vector<Decoder> dec(nthreads);
    std::vector<bam1_t*> batch;
    try {
        for (size_t c = 0; c < so.size(); ++c) { so[c].cs = &cs; so[c].c = (int)c; so[c].ref = &ref; so[c].open(outs[c]); }
        for (auto& d : dec) { d.opt = &o; d.ref = &ref; d.cs = &cs; }

        // reads decoded per parallel batch; small enough to bound memory (bam1_t buffers keep
        // their largest size), large enough to amortise the per-batch thread start-up
        const size_t BATCH = std::max<size_t>(256, 64 * (size_t)nthreads);
        batch.resize(BATCH);
        for (auto& b : batch) b = bam_init1();
        std::vector<int32_t> batch_idx(BATCH);
        std::vector<std::vector<Call>> calls(BATCH);
        int32_t read_index = 0;   // readONTbam's read index: position among *kept* reads
        int64_t last_pos0 = -1;   // 0-based POS of the last record fetched
        int ret = 0;
        bool eof = false;
        while (!eof) {
            size_t n = 0;
            while (n < BATCH) {
                ret = sam_itr_next(fp, itr, batch[n]);
                if (ret < 0) { eof = true; break; }
                ++R.n_records;
                const bam1_t* b = batch[n];
                if (b->core.tid == tid) {
                    if (b->core.pos < last_pos0) throw std::runtime_error(bam + " is not coordinate-sorted");
                    last_pos0 = b->core.pos;
                }
                // --- readONTbam filtering ---
                if (b->core.n_cigar == 0) continue;                       // !is.na(cigar)
                if (b->core.tid != tid) continue;                         // rname %in% chr
                const uint16_t fl = b->core.flag;
                const bool sec = fl & BAM_FSECONDARY, sup = fl & BAM_FSUPPLEMENTARY;
                if (sec || sup) {
                    uint8_t* mlp = bam_aux_get(b, "ML");
                    bool has_tag = mlp && mlp[0] == 'B' && bam_auxB_len(mlp) > 0;
                    bool hard = false;
                    const uint32_t* cg = bam_get_cigar(b);
                    for (uint32_t k = 0; k < b->core.n_cigar; ++k) if (bam_cigar_op(cg[k]) == BAM_CHARD_CLIP) hard = true;
                    bool valid_seq = b->core.l_qseq > 0;
                    if (!(has_tag && !hard && valid_seq)) continue;
                }
                batch_idx[n] = ++read_index;
                ++n;
            }
            if (ret < -1) throw std::runtime_error("error reading " + bam + " (truncated/corrupt?)");
            if (n > 0) {
                R.n_kept += n;
                auto work = [&](int tix) {
                    size_t lo = n * tix / nthreads, hi = n * (tix + 1) / nthreads;
                    for (size_t k = lo; k < hi; ++k) { calls[k].clear(); dec[tix].decode(batch[k], calls[k]); }
                };
                if (nthreads == 1 || n < 64) { for (size_t k = 0; k < n; ++k) { calls[k].clear(); dec[0].decode(batch[k], calls[k]); } }
                else {
                    std::vector<std::thread> th;
                    for (int t = 1; t < nthreads; ++t) th.emplace_back(work, t);
                    work(0);
                    for (auto& x : th) x.join();
                }
                // merge in read order (keeps read indices ascending per site)
                for (size_t k = 0; k < n; ++k) {
                    const int32_t rid = batch_idx[k];
                    for (const Call& cl : calls[k]) { so[cl.ctx].add(cl.site, cl.mod, rid); ++R.n_calls; }
                }
            }
            // no later record starts before last_pos0, and calls land at gp >= read start
            if (last_pos0 >= 0) for (auto& s : so) s.flush_before(last_pos0 + 1);
        }
        for (auto& s : so) s.finish();
    } catch (std::exception& e) {
        R.err = e.what();
        for (auto& s : so) s.abandon();
    }
    for (auto& b : batch) if (b) bam_destroy1(b);
    hts_itr_destroy(itr); hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
    if (!R.err.empty()) return R;
    int64_t lit = 0; for (auto& d : dec) lit += d.literal_checks;
    if (o.check_literal) fprintf(stderr, "[%s %s] literal check OK on %lld calls\n", label.c_str(), cs.chrom.c_str(), (long long)lit);
    R.ok = true;
    R.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return R;
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    Options o = parse_args(argc, argv);
    auto T0 = std::chrono::steady_clock::now();
    mkdir(o.outdir.c_str(), 0755);

    faidx_t* fai = fai_load(o.ref_fa.c_str());
    if (!fai) { fprintf(stderr, "cannot load FASTA index for %s\n", o.ref_fa.c_str()); return 1; }
    const bool all_mode = o.chroms.empty();
    std::vector<std::string> chroms = o.chroms;
    if (all_mode) for (int i = 0; i < faidx_nseq(fai); ++i) chroms.push_back(faidx_iseq(fai, i));
    {   // drop repeated names, keep order
        std::vector<std::string> u;
        for (auto& c : chroms) if (std::find(u.begin(), u.end(), c) == u.end()) u.push_back(c);
        chroms.swap(u);
    }
    for (auto& c : chroms)
        if (!faidx_has_seq(fai, c.c_str())) { fprintf(stderr, "chromosome %s not found in %s\n", c.c_str(), o.ref_fa.c_str()); return 1; }
    // whole-chromosome output goes to OUTDIR/<chrom>/ (the layout loadONTSample() in R/readONTbam_fast.R reads);
    // a single window keeps the flat OUTDIR/<label>_<ctx>.tsv layout
    const bool window_mode = o.win_start >= 0;
    const bool multi = !window_mode;

    const int nb = (int)o.bams.size();
    int jobs = o.jobs > 0 ? o.jobs : std::min(o.threads, nb);
    jobs = std::max(1, std::min(jobs, nb));
    const int per = std::max(1, o.threads / jobs);
    int nfail = 0, nok = 0;

    for (const std::string& chrom : chroms) {
        auto TC = std::chrono::steady_clock::now();
        // reference chromosome (uppercase, like a BSgenome/DNAString)
        hts_pos_t len = 0;
        char* seq = faidx_fetch_seq64(fai, chrom.c_str(), 0, HTS_POS_MAX, &len);
        if (!seq || len <= 0) { fprintf(stderr, "cannot read %s from %s\n", chrom.c_str(), o.ref_fa.c_str()); return 1; }
        std::string ref(seq, (size_t)len);
        free(seq);
        for (auto& ch : ref) ch = (char)toupper((unsigned char)ch);

        const bool window = o.win_start >= 0;
        const int64_t ws = window ? o.win_start : 1;
        const int64_t we = window ? std::min<int64_t>(o.win_end, (int64_t)ref.size()) : (int64_t)ref.size();
        std::string fetch = !o.fetch.empty() ? o.fetch
                          : window ? chrom + ":" + std::to_string(ws + 1) + "-" + std::to_string(we)
                                   : chrom;
        ChromSites cs;
        try { cs = build_sites(chrom, ref, o.contexts, ws, we); }
        catch (std::exception& e) { fprintf(stderr, "%s: %s\n", chrom.c_str(), e.what()); return 1; }
        for (auto& S : cs.ctx) fprintf(stderr, "[readontbam_fast] %s sites in %s:%lld-%lld: %zu\n", S.ctx.c_str(),
                                       chrom.c_str(), (long long)ws, (long long)we, S.pos.size());
        const std::string dir = multi ? o.outdir + "/" + chrom : o.outdir;
        if (multi) mkdir(dir.c_str(), 0755);

        std::vector<BamResult> res(nb);
        std::atomic<int> next{0};
        std::mutex mu;
        auto worker = [&]() {
            while (true) {
                int i = next++;
                if (i >= nb) break;
                std::vector<std::string> outs;
                for (auto& S : cs.ctx) outs.push_back(dir + "/" + o.labels[i] + "_" + S.ctx + ".tsv");
                res[i] = process_bam(o, o.labels[i], o.bams[i], fetch, outs, ref, cs, per, all_mode);
                std::lock_guard<std::mutex> lk(mu);
                const BamResult& r = res[i];
                if (r.skipped) fprintf(stderr, "[readontbam_fast] %s %s: outputs exist, skipped\n", r.label.c_str(), chrom.c_str());
                else if (r.absent) fprintf(stderr, "[readontbam_fast] %s %s: not in BAM header, skipped\n", r.label.c_str(), chrom.c_str());
                else if (r.ok) fprintf(stderr, "[readontbam_fast] %s %s: %lld records, %lld kept reads, %lld site calls, %.1fs\n",
                                       r.label.c_str(), chrom.c_str(), (long long)r.n_records, (long long)r.n_kept,
                                       (long long)r.n_calls, r.secs);
                else fprintf(stderr, "[readontbam_fast] FAILED %s %s: %s\n", r.label.c_str(), chrom.c_str(), r.err.c_str());
            }
        };
        std::vector<std::thread> th;
        for (int j = 0; j < jobs; ++j) th.emplace_back(worker);
        for (auto& t : th) t.join();
        for (auto& r : res) { if (r.ok) ++nok; else ++nfail; }
        if (multi) fprintf(stderr, "[readontbam_fast] %s done in %.1fs\n", chrom.c_str(),
                           std::chrono::duration<double>(std::chrono::steady_clock::now() - TC).count());
    }
    fai_destroy(fai);
    fprintf(stderr, "[readontbam_fast] done: %d/%d BAM x chromosome jobs ok, %.1fs total\n", nok, nok + nfail,
            std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count());
    return nfail ? 1 : 0;
}
