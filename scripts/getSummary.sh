#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that generates a summary of all the SNVs in all the samples

###############################################################################################################

#initialize ODIR (output dir) variable
if  [ "$#" -lt 1 ] ; then ODIR=$HP_ODIR ; else ODIR=$1 ; fi

#get read count and mtDN-CN stats 
if [ $HP_CN ]  && [ $HP_CN -ne 0 ] || [ ! -s $ODIR/count.tab ] ; then
    grep -v "^#" $HP_IN | cut -f3 | perl -ane 'print "$F[0].count\n";' | xargs cat | cut -f1,2,3,4 | uniq.pl -i 0 | getCN.pl > $ODIR/count.tab
fi
if [ $HP_I -lt 1 ] ; then exit 0 ; fi

###########################################################
# get 1st iteration stats

M="${2:-$HP_M}" 
###########################################################

#count,cvg
grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.cvg.stat|"  | xargs cat | uniq.pl -i 0  > $ODIR/cvg.tab

#March 7th 2023
#awk '{print $3}' $HP_IN | sed "s|$|.$M.00.vcf|" | xargs cat | grep -v "^##sample="  |   uniq.pl | bedtools sort -header | eval $HP_FRULE   > $ODIR/$M.00.concat.vcf

#Sept 22nd 2023; removed " | eval $HP_FRULE"
grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$M.00.vcf|" | xargs cat | grep -v "^##sample="  |   uniq.pl | bedtools sort -header   > $ODIR/$M.00.concat.vcf

snpSort.sh $ODIR/$M.00.concat
cat $ODIR/$M.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $ODIR/$M.00.AF.histo

#haplogroups; sed "s|\.MT||" 
if [ "$HP_O" == "Human" ] ; then
  grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$M.haplogroup|" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | \
    perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' | \
    tee $ODIR/$M.haplogroup.tab | \
    perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $ODIR/$M.haplogroup1.tab

  #haplocheck
  grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$M.haplocheck|" | xargs cat  | uniq.pl | sed 's|^"Sample"|"Run"|' | sed 's|"||g' | sed 's| ||g' > $ODIR/$M.haplocheck.tab
fi

#fasta
grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$M.fa|"        | xargs cat > $ODIR/$M.fa
samtools faidx  $ODIR/$M.fa

#cvg
grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$M.merge.bed|" | xargs cat > $ODIR/$M.merge.bed

#sv's
if [[ ! -z "$HP_V" ]] ; then
  grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$HP_V.00.vcf|" | xargs cat | uniq.pl | bedtools sort -header > $ODIR/$HP_V.00.concat.vcf
fi

#snv counts
snpCount.sh $M $HP_T1
snpCount.sh $M $HP_T2
snpCount.sh $M $HP_T3

#cleanup
rm -f fastp.html fastp.json

##########################################################
# get 2nd iteration stats

if [ $HP_I -lt 2 ] ; then exit 0 ; fi
if [ $M == "mutserve" ] ; then exit 0 ; fi

MM=$M.$M
grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$M.cvg.stat|" | xargs cat | uniq.pl -i 0  > $ODIR/$S.cvg.tab

#Sept 22nd 2023 ; removed "| eval $HP_FRULE " 
#awk '{print $3}' $HP_IN | sed "s|$|.$MM.00.vcf|"  | xargs cat |  grep -v "^##sample=" |   grep -v "^##bcftools_annotateCommand" | uniq.pl | bedtools sort -header  | eval $HP_FRULE  > $ODIR/$MM.00.concat.vcf  
grep -v "^#" $HP_IN | cut -f3 | sed "s|$|.$MM.00.vcf|"  | xargs cat |  grep -v "^##sample=" |   grep -v "^##bcftools_annotateCommand" | uniq.pl | bedtools sort -header   > $ODIR/$MM.00.concat.vcf  

snpSort.sh $ODIR/$MM.00.concat
cat $ODIR/$MM.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $ODIR/$MM.00.AF.histo

#snv counts
snpCount.sh $MM $HP_T1
snpCount.sh $MM $HP_T2
snpCount.sh $MM $HP_T3
