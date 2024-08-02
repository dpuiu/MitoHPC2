#!/bin/bash -e

S="$1"
M="mutect2"
OS="out/$S/$S.$M"

if [ ! -s $OS.00.max.vcf ] ; then
  test -s $OS.vcf

  bcftools norm -m-any -f $HP_RDIR/$HP_MT.fa  $OS.vcf   | fix${M}Vcf.pl -file $HP_RDIR/$HP_MT.fa  | bedtools sort -header | \
  filterVcf.pl -sample $S -source $M -header $HP_SDIR/$M.vcf -depth $HP_DP | uniqVcf.pl | bedtools sort -header  | tee $OS.00.vcf | \
  maxVcf.pl | bedtools sort -header > $OS.00.max.vcf
fi
