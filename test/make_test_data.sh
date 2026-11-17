#!/usr/bin/env bash
# Cut a small self-test dataset from your own data.
#   test/make_test_data.sh in.bam ref.fa CHROM:START-END [out_prefix]
# Writes <out_prefix>.bam(.bai) with the primary/secondary/supplementary records overlapping
# the region, and <out_prefix>.fa(.fai) with the WHOLE chromosome CHROM (readONTbam and the
# self-test enumerate every cytosine of the chromosome). Pick a region with a few dozen
# reads: readONTbam() in R needs ~3 s per read and context.
set -euo pipefail
[ $# -ge 3 ] || { sed -n '2,7p' "$0"; exit 1; }
bam=$1; ref=$2; region=$3; out=${4:-mytest}
chrom=${region%%:*}
samtools view -b -o "$out.bam" "$bam" "$region"
samtools index "$out.bam"
samtools faidx "$ref" "$chrom" > "$out.fa"
samtools faidx "$out.fa"
echo "$(samtools view -c -F 0x900 "$out.bam") primary reads in $out.bam; reference $chrom in $out.fa"
echo "Run: Rscript test/selftest.R CG $out.bam $out.fa"
