#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that aligns Illumina pair-end reads to a genome refenece using bwa mem

# Input:  genome bwa index, read FASTQ/FASTA files
# Output: alignment .bam, ,bai, .idxstats files

###############################################################################################################

I=$1  # input prefix
O=$2  # output prefix

test -s $HP_RDIR/$HP_RNAME.fa

if [ ! -s $O.bam ] ; then
  test -s $HP_RDIR/$HP_RNAME.amb
  test -s $HP_RDIR/$HP_RNAME.ann
  test -s $HP_RDIR/$HP_RNAME.bwt
  test -s $HP_RDIR/$HP_RNAME.pac
  test -s $HP_RDIR/$HP_RNAME.sa

  test -s ${I}_1.f*q*
  test -s ${I}_2.f*q*

  bwa mem $HP_RDIR/$HP_RNAME ${I}_[12].f*q* -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    head -100000000 | \
    samtools view -bu | \
    samtools sort -m $HP_MM -@ $HP_P -o $O.bam
fi

if [ ! -s $I.bam.bai ] ; then
  samtools index $O.bam -@ $HP_P
fi

if [ ! -s $I.idxstats ] ; then
  samtools idxstats $O.bam > $O.idxstats
fi
