#!/usr/bin/env bash
set -euo pipefail

# Usage: $0 <sample> <HPRC_file>
#
# sample     Sample ID (HG\d+ or NA\d+)
# HPRC_file  HPRC remote BAM file path starting with s3://human-pangenomics/
#
# Required environment variables:
#   HP_RDIR      reference directory
#   HP_MT        mitochondrial reference name
#   HP_RMT       mitochondrial chromosome
#   HP_MINLEN    minimum alignment length
#   HP_MINPC     minimum alignment percentage
#   HP_MINID     minimum alignment identity
#
# Required programs:
#   s5cmd samtools minimap2 filterAlignments.pl

S=$1
F=$2
P=2

# MT BAM already exists: nothing more to do.
if [[ -s "$S.MT.bam" ]]; then
    exit 0
fi

# ----------------------------------------------------------------------
# Download BAM
# ----------------------------------------------------------------------

if [[ ! -s "$S.bam" ]]; then
    s5cmd cp "$F" "$S.bam"
fi

# Download BAI file if it exists.
if s5cmd ls "$F.bai" >/dev/null 2>&1; then
    s5cmd cp "$F.bai" "$S.bam.bai"
fi

# ----------------------------------------------------------------------
# Get sample read count
# ----------------------------------------------------------------------

if [[ ! -s "$S.count" ]]; then
    samtools view -c -@ "$P" "$S.bam" > "$S.count"
fi

# ----------------------------------------------------------------------
# Extract candidate mitochondrial reads
# ----------------------------------------------------------------------

if [[ ! -s "$S.mt.bam" ]]; then
    if [[ -s "$S.bam.bai" ]]; then
        samtools view -@ "$P" "$S.bam" "$HP_RMT" > "$S.mt.bam"
    else
        samtools fastq "$S.bam" |
            minimap2 \
                --eqx \
                -ax map-hifi \
                -t "$P" \
                "$HP_RDIR/$HP_MT.fa" \
                /dev/stdin |
            samtools view -b -F 4 \
            > "$S.mt.bam"
    fi
fi

# ----------------------------------------------------------------------
# Select reads with sufficient mitochondrial alignment
# ----------------------------------------------------------------------

if [[ ! -s "$S.MT.ids" ]]; then
    samtools fastq "$S.mt.bam" |
        minimap2 \
            --eqx \
            -ax map-hifi \
            -t "$P" \
            "$HP_RDIR/$HP_MT2.fa" \
            /dev/stdin |
        samtools view -h -F 4 |
        filterAlignments.pl \
            --minlen "$HP_MINLEN" \
            --minpc "$HP_MINPC" \
            --minid "$HP_MINID" |
        cut -f1 \
        > "$S.MT.ids"
fi

# ----------------------------------------------------------------------
# Extract and sort selected reads
# ----------------------------------------------------------------------

if [[ ! -s "$S.MT.bam" ]]; then
    samtools view \
        -N "$S.MT.ids" \
        -b "$S.mt.bam" |
        samtools sort \
            -@ "$P" \
            -o "$S.MT.bam"

    samtools index "$S.MT.bam"
    samtools idxstats "$S.MT.bam" > "$S.MT.idxstats"
fi

# ----------------------------------------------------------------------
# Get MT read count
# ----------------------------------------------------------------------

if [[ ! -s "$S.MT.count" ]]; then
    paste "$S.count" "$S.MT.idxstats" |
        head -1 |
        perl -ane 'print "$F[0]\t$F[0]\t$F[3]\n";' |
        sed "s|^|$S\t|" \
        > "$S.MT.count"
fi

# ----------------------------------------------------------------------
# Cleanup
# ----------------------------------------------------------------------

if [[ -s "$S.MT.bam" ]]; then
    rm -f \
        "$S.bam" \
        "$S.bam.bai" \
        "$S.mt.bam" \
        "$S.MT.ids"
fi
```

