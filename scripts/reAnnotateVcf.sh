#!/usr/bin/env bash
set -eux

##############################################################################################################

# Program that re-annotates a VCF file

# Input:  VCF file
# Output: VCF file

###############################################################################################################

I=$1
O=$2

test -s $I

cat $I | deAnnotateVcf.pl > $O
bgzip -f $O
tabix -f $O.gz

bcftools annotate -a $HP_RDIR/dbSNP.vcf.gz     -c "ID" $O.gz |\
  bcftools annotate -a $HP_RDIR/Hypervariable.bed.gz      -c "CHROM,FROM,TO,Hypervariable" -h <(echo '##INFO=<ID=Hypervariable,Number=1,Type=String,Description="Hypervariable">') |\
  bcftools annotate -a $HP_RDIR/Homopolymer.bed.gz      -c "CHROM,FROM,TO,Homopolymer"	-h <(echo '##INFO=<ID=Homopolymer,Number=1,Type=Flag,Description="Homopolymer">') |\
  bcftools annotate -a $HP_RDIR/Hotspot.bed.gz      -c "CHROM,FROM,TO,Hotspot"       -h <(echo '##INFO=<ID=Hotspot,Number=1,Type=Flag,Description="Hotspot">') |\
  bcftools annotate -a $HP_RDIR/CDS.bed.gz     -c "CHROM,FROM,TO,CDS"     -h <(echo '##INFO=<ID=CDS,Number=1,Type=String,Description="CDS">') |\
  bcftools annotate -a $HP_RDIR/COMPLEX.bed.gz -c "CHROM,FROM,TO,COMPLEX" -h <(echo '##INFO=<ID=COMPLEX,Number=1,Type=String,Description="COMPLEX">') |\
  bcftools annotate -a $HP_RDIR/RNR.bed.gz     -c "CHROM,FROM,TO,RNR"     -h <(echo '##INFO=<ID=RNR,Number=1,Type=String,Description="rRNA">') |\
  bcftools annotate -a $HP_RDIR/TRN.bed.gz     -c "CHROM,FROM,TO,TRN"     -h <(echo '##INFO=<ID=TRN,Number=1,Type=String,Description="tRNA">') |\
  bcftools annotate -a $HP_RDIR/DLOOP.bed.gz   -c "CHROM,FROM,TO,DLOOP"   -h <(echo '##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="DLOOP">') |\
  bcftools annotate -a $HP_RDIR/MCC.bed.gz -c "CHROM,FROM,TO,MCC"         -h <(echo '##INFO=<ID=MCC,Number=1,Type=String,Description="New YALE protein_gene_missense_constraint; missense_OEUF">') |\
  annotateVcf.pl - $HP_RDIR/HG.vcf.gz        |\
  annotateVcf.pl - $HP_RDIR/NUMT.vcf.gz      |\
  annotateVcf.pl - $HP_RDIR/NONSYN.vcf.gz    |\
  annotateVcf.pl - $HP_RDIR/STOP.vcf.gz      |\
  annotateVcf.pl - $HP_RDIR/MITOMAP.vcf.gz   |\
  annotateVcf.pl - $HP_RDIR/MITIMPACT.vcf.gz |\
  annotateVcf.pl - $HP_RDIR/HelixMTdb.vcf.gz |\
  annotateVcf.pl - $HP_RDIR/gnomAD31.vcf.gz  |\
  annotateVcf.pl - $HP_RDIR/UKB_mutect2_Nature_20230816.vcf.gz | \
  annotateVcf.pl - $HP_RDIR/MLC.vcf.gz       > $O
rm $O.gz
