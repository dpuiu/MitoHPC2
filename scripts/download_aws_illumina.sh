#!/usr/bin/env bash
set -exuo pipefail

# Usage: $0 <sample> <HPRC_file>
#
# sample     Sample ID (HG\d+ or NA\d+)
# HPRC_file  HPRC remote CRAM/BAM path starting with s3://human-pangenomics/
#
# Required environment variables:
#   HP_RMT      mitochondrial chromosome
#   HP_RNUMT    NUMT chromosome/region
#   HP_MTLEN    mitochondrial reference length
#   HP_RDIR     reference directory
#   HP_RNAME    reference name
#
# Required programs:
#   s5cmd samtools idxstats2count.pl filterSam.pl

S=$1
F=$2
P=2

# MT BAM already exists: nothing more to do.
if [[ -s "$S.MT.bam.bai" ]]; then
    exit 0
fi

# ----------------------------------------------------------------------
# Download CRAM/BAM
# ----------------------------------------------------------------------

if [[ ! -s "$S.cram" ]]; then
    s5cmd cp "$F" "$S.cram"
fi

# ----------------------------------------------------------------------
# Download CRAI if available; otherwise create it locally
# ----------------------------------------------------------------------

if [[ ! -s "$S.cram.crai" ]]; then
  if s5cmd ls "$F.crai" >/dev/null 2>&1; then
      s5cmd cp "$F.crai" "$S.cram.crai"
  else
      samtools index -@ "$P" "$S.cram"
  fi
fi

# ----------------------------------------------------------------------
# Count reads mapped to MT
# ----------------------------------------------------------------------

if [[ ! -s "$S.MT.count" ]]; then
    samtools idxstats -@ "$P" "$S.cram" | \
        tee $S.idxstats | \
        idxstats2count.pl \
            -sample "$S" \
            -chrM "$HP_RMT" \
            > "$S.MT.count"
fi

# ----------------------------------------------------------------------
# Extract MT and NUMT reads and filter them
# ----------------------------------------------------------------------

if [[ ! -s "$S.MT.bam" ]]; then
    samtools view \
        -h \
        "$S.cram" \
        "$HP_RMT:1-$HP_MTLEN" \
        "$HP_RNUMT" \
        -F 0x90C \
        -T "$HP_RDIR/$HP_RNAME.fa" |
        filterSam.pl "$HP_RMT:1-$HP_MTLEN" "$HP_RNUMT" |
        samtools view -b \
        > "$S.MT.bam"
fi

if [[ ! -s "$S.MT.bam.bai" ]]; then
    samtools index "$S.MT.bam"
    samtools idxstats "$S.MT.bam" > "$S.MT.idxstats"
fi

# get subsampling rate
if [[ ! -s "$S.MT.r" ]]; then
  samtools stats "$S.MT.bam" | \
      tee "$S.MT.stats" | \
      grep -m  1 'bases mapped' | cut -f3 | \
      perl -ane 'print ($ENV{HP_MTLEN}*$ENV{HP_C}/$F[0]);' > "$S.MT.r"
  
  #samtools view -s `cat $S.MT.r`  $S.MT.bam -F 0x90C 
fi

# ----------------------------------------------------------------------
# Remove downloaded CRAM and index after MT extraction
# ----------------------------------------------------------------------

#if [[ -s "$S.MT.bam" ]]; then
#    rm -f "$S.cram" "$S.cram.crai"
#fi
