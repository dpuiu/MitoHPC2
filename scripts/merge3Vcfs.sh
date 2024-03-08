#!/bin/bash -eux

###############################################################################

M1=mutect2   ; M3=freebayes  ; M2=varscan     ;M=merge3
OM1="out.$M1"; OM2="out.$M2" ; OM3="out.$M3" ; OM=$HP_ODIR
test -s $OM1/$M1.00.concat.vcf ; test -s $OM2/$M2.00.concat.vcf ; test -s $OM3/$M3.00.concat.vcf

cat $OM1/$M1.00.concat.vcf $OM2/$M2.00.concat.vcf $OM3/$M3.00.concat.vcf | egrep -v '^##bcftools|^##source' | uniqVcf.pl  -min 2 | bedtools sort -header > $OM/$M.00.concat.vcf
snpCount.sh $M $HP_T1
snpCount.sh $M $HP_T2
snpCount.sh $M $HP_T3

###############################################################################

MM1="mutect2.mutect2" ; MM3="freebayes.freebayes" ; MM2="varscan.varscan" ; MM="merge3.merge3"
test -s $OM1/$MM1.00.concat.vcf ; test -s $OM2/$MM2.00.concat.vcf ; test -s $OM3/$MM3.00.concat.vcf

cat $OM1/$MM1.00.concat.vcf $OM2/$MM2.00.concat.vcf $OM3/$MM3.00.concat.vcf | egrep -v '^##bcftools|^##source' | uniqVcf.pl  -min 2 | bedtools sort -header > $OM/$MM.00.concat.vcf
snpCount.sh $MM $HP_T1
snpCount.sh $MM $HP_T2
snpCount.sh $MM $HP_T3

