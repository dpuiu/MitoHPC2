#!/bin/bash -eux

###############################################################################

if [ $HP_I -lt "1" ] ; then exit 0 ; fi

test -s $HP_ODIR/$HP_M1.00.concat.vcf
test -s $HP_ODIR/$HP_M2.00.concat.vcf
test -s $HP_ODIR/$HP_M3.00.concat.vcf

cat $HP_ODIR/$HP_M1.00.concat.vcf $HP_ODIR/$HP_M2.00.concat.vcf $HP_ODIR/$HP_M3.00.concat.vcf | egrep -v '^##bcftools|^##source' | uniqVcf.pl  -min 2 | bedtools sort -header > $HP_ODIR/$HP_M.00.concat.vcf
snpCount.sh $HP_M $HP_T1
snpCount.sh $HP_M $HP_T2
snpCount.sh $HP_M $HP_T3

echo $HP_I
if [ $HP_I -lt "2" ] ; then exit 0 ; fi

###############################################################################

test -s $HP_ODIR/$HP_M1.$HP_M1.00.concat.vcf
test -s $HP_ODIR/$HP_M2.$HP_M2.00.concat.vcf
test -s $HP_ODIR/$HP_M3.$HP_M3.00.concat.vcf

cat $HP_ODIR/$HP_M1.$HP_M1.00.concat.vcf $HP_ODIR/$HP_M2.$HP_M2.00.concat.vcf $HP_ODIR/$HP_M3.$HP_M3.00.concat.vcf | egrep -v '^##bcftools|^##source' | uniqVcf.pl  -min 2 | bedtools sort -header > $HP_ODIR/$HP_M.$HP_M.00.concat.vcf
snpCount.sh $HP_M.$HP_M $HP_T1
snpCount.sh $HP_M.$HP_M $HP_T2
snpCount.sh $HP_M.$HP_M $HP_T3
