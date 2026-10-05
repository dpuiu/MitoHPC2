#!/bin/bash -eux

# Usage: $0 <sample> <HPRC_file> 
#
# sample     Sample ID (HG\d+ or NA\d+)
# HPRC_file  HPRC renote CRAM/BAM file path starting with s3://human-pangenomics/

S=$1		# sample id
F=$2		# HPRC remote file name

MT=chrM
P=2
RNUMT="chr1:629084-634672 chr17:22521208-22521639"
RMT=chrM
RNAME=hs38DH
MTLEN=16569

if [ -s $S.MT.bam ]; then
  exit 0
fi

# Download HPRC CRAM & CRAI files
if [ ! -s $S.cram ]; then
  s5cmd cp $F      $S.cram
  s5cmd cp $F.crai $S.cram.crai
fi

# Generate CRAI file if not avail
if [ ! -s $S.cram.crai ]; then
  samtools index -@ $P $S.cram
fi

# Get the count file
if [ ! -s $S.MT.count ]; then
  samtools idxstats -@ $P $S.cram | \
    ./idxstats2count.pl -sample $S -chrM $RMT > $S.MT.count
fi

# Extract reads and align to the mitochondrial reference
if [ ! -s $S.MT.bam ]; then
  #samtools view  $S.cram $RMT $RNUMT -b > $S.MT.bam
  samtools view -h $S.cram $RMT:1-$MTLEN $RNUMT -F 0x90C | \
    ./filterSam.pl $RMT:1-$MTLEN $RNUMT | samtools view -b  > $S.MT.bam
  samtools index $S.MT.bam
  samtools idxstats $S.MT.bam > $S.MT.idxstats
  #samtools depth $S.MT.bam -r $RMT > $S.MT.depth
fi

# Cleanup
if [ -s $S.MT.bam ]; then
  rm $S.cram $S.cram.crai
fi
