#!/bin/bash 
set -euxo pipefail

# Usage: $0 <sample> <HPRC_file> 
#
# sample     Sample ID (HG\d+ or NA\d+)
# HPRC_file  HPRC renote CRAM/BAM file path starting with s3://human-pangenomics/

S=$1		# sample id
F=$2		# HPRC remote file name
P=2

if [ -s $S.MT.bam ]; then
  exit 0
fi

# Download HPRC CRAM file
if [ ! -s $S.cram ]; then
  s5cmd cp $F $S.cram
fi

# Download/GENERATE HPRC CRAI file if exists
if s5cmd ls "$F.crai" >/dev/null 2>&1  ; then
    s5cmd cp $F.crai $S.cram.crai
else
    samtools index -@ $P $S.cram
fi

# Get the count file
if [ ! -s $S.MT.count ]; then
  samtools idxstats -@ $P $S.cram | \
    idxstats2count.pl -sample $S -chrM $HP_RMT > $S.MT.count
fi

# Extract reads and align to the mitochondrial reference
if [ ! -s $S.MT.bam ]; then
  #samtools view  $S.cram $RMT $RNUMT -b > $S.MT.bam
  samtools view -h $S.cram $HP_RMT:1-$HP_MTLEN $HP_RNUMT -F 0x90C -T $HP_RDIR/$HP_RNAME.fa | \
    filterSam.pl $HP_RMT:1-$HP_MTLEN $HP_RNUMT | samtools view -b  > $S.MT.bam
  samtools index $S.MT.bam
  samtools idxstats $S.MT.bam > $S.MT.idxstats
  #samtools depth $S.MT.bam -r $RMT > $S.MT.depth
fi

# Cleanup
if [ -s $S.MT.bam ]; then
  rm $S.cram $S.cram.crai
fi

