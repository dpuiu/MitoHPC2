#!/usr/bin/env bash
set -ex
##############################################################################################################

# Program that merges multiple MitoHPC results
#  input directories must be given as command line arguments

##############################################################################################################

if [ "$#" -lt 1 ]; then 
  echo "ERROR"
  exit 1
fi

#########################################

mkdir -p $HP_ODIR

#for F in count.tab cvg.tab $HP_M.cvg.tab $HP_M.haplocheck.tab $HP_M.haplogroup1.tab $HP_M.haplogroup.tab $HP_M.merge.bed
#do
#  find $@ -name $F | xargs cat | uniq.pl  | { sed -u 1q; sort; } > $HP_ODIR/$F
#done
#cut -f1 $HP_ODIR/cvg.tab | tail -n +2 > $HP_IN

#for F in $HP_M.fa
#do
#   find $@ -name $F | xargs cat > $HP_ODIR/$F
#   samtools faidx $HP_ODIR/$F
#done

#for S in $HP_M.$HP_M
for M in $HP_M
do
  #find $@ -name $M.00.concat.vcf | xargs cat | bedtools sort -header > $HP_ODIR/$M.00.concat.vcf
  #cat $HP_ODIR/$M.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $HP_ODIR/$M.00.AF.histo

  snpCount.sh $M $HP_T1 
  snpCount.sh $M $HP_T2 
  snpCount.sh $M $HP_T3 
done
 	

