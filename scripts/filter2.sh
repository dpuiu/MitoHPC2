#!/bin/bash -x

S="$1"
M="mutect2"
OS="out/$S/$S.$M"
OSS="$OS.$M"

if [ ! -s $OSS.00.orig.vcf ] ; then
  test -s $OS.max.vcf
  test -s $OS.fa
  test -s $OSS.vcf

  cat $OSS.vcf | bcftools norm -m-any -f $OS.fa  | fix${M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header | \
    fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -rlen $HP_MTLEN -mfile $OS.max.vcf  | \
    filterVcf.pl -sample $S -source $M -header $HP_SDIR/$M.vcf  -depth $HP_DP | bedtools sort  -header > $OSS.00.orig.vcf

  #annotateVcf.sh  $OSS.00.orig.vcf
fi


