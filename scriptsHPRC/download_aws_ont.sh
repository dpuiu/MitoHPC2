#!/bin/bash -eux

# Usage: $0 <HPRC_file> <sample>
#
# HPRC_file  HPRC BAM file path relative to s3://human-pangenomics/
# sample     Sample ID

S=$1                 # sample id
F=$2                 # HPRC remote file name

MT=chrM
MT2=chrM2
P=2
MINLEN=6000
MINPC=0.95
MINID=0.95
RMT=chrM

# Test references exist
test -s "$MT.fa"  || exit 1
test -s "$MT2.fa" || exit 1

if [ -s $S.MT.bam ]; then
  exit 0
fi

# Download HPRC BAM
if [ ! -s "$S.bam" ]; then
  s5cmd cp $F $S.bam
fi

# Get sample stats
if [ ! -s "$S.count" ]; then
  samtools view -c -@ "$P" "$S.bam" > "$S.count"
fi

# Extract reads and align to the mitochondrial reference
if [ ! -s "$S.mt.bam" ]; then
  samtools fastq "$S.bam" | \
    minimap2 --eqx -ax map-ont "$MT.fa" /dev/stdin -t "$P" | \
    samtools view -b -F 4 > "$S.mt.bam"
fi

# Select reads with >= MIN_LENGTH bases of mitochondrial alignment
if [ ! -s "$S.MT.ids" ]; then
  samtools fastq "$S.mt.bam" | \
    minimap2 --eqx  -ax map-hifi "$MT2.fa" /dev/stdin -t "$P" | \
    samtools view -h -F 4  | \
    ./filterAlignments.pl --minlen "$MINLEN" --minpc "$MINPC" --minid "$MINID" | \
    cut -f1 > "$S.MT.ids"
fi

# Extract and sort selected reads
if [ ! -s "$S.MT.bam" ]; then
  samtools view -N "$S.MT.ids" -b "$S.mt.bam" |\
    samtools sort -@ "$P" -o "$S.MT.bam"
  samtools index $S.MT.bam
  samtools idxstats $S.MT.bam > $S.MT.idxstats
fi

if [ ! -s "$S.MT.count" ]; then
  paste $S.count $S.MT.idxstats | head -1 | perl -ane 'print "$F[0]\t$F[0]\t$F[3]\n";' | sed "s|^|\$S\t|" > $S.MT.count
fi

# Cleanup
if [ -s $S.MT.bam ]; then
  rm $S.cram $S.cram.crai $S.mt.bam $S.MT.ids
fi

