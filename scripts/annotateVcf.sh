#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that annotates a VCF file

# Input:  VCF file
# Output: VCF file

###############################################################################################################

O=$1   # input/output vcf file

test  -s $O

bgzip -f $O ;tabix -f $O.gz

if [ -s $HP_RDIR/HV.bed.gz ]    ; then bcftools annotate -a $HP_RDIR/HV.bed.gz    $O.gz -c "CHROM,FROM,TO,Hypervariable"  -h <(echo '##INFO=<ID=Hypervariable,Number=1,Type=String,Description="Hypervariable">') > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/HP.bed.gz ]    ; then bcftools annotate -a $HP_RDIR/HP.bed.gz    $O.gz -c "CHROM,FROM,TO,Homopolymer"  -h <(echo '##INFO=<ID=Homopolymer,Number=0,Type=Flag,Description="Homoloplymer">')  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/HS.bed.gz ]    ; then bcftools annotate -a $HP_RDIR/HS.bed.gz    $O.gz -c "CHROM,FROM,TO,Hotspot"  -h <(echo '##INFO=<ID=Hotspot,Number=0,Type=Flag,Description="Hotspot">')      > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/CDS.bed.gz ]   ; then bcftools annotate -a $HP_RDIR/CDS.bed.gz   $O.gz -c "CHROM,FROM,TO,CDS" -h <(echo '##INFO=<ID=CDS,Number=1,Type=String,Description="CDS">')   > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/COMPLEX.bed.gz ]   ; then bcftools annotate -a $HP_RDIR/COMPLEX.bed.gz   $O.gz -c "CHROM,FROM,TO,COMPLEX" -h <(echo '##INFO=<ID=COMPLEX,Number=1,Type=String,Description="COMPLEX">')   > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/RNR.bed.gz ]   ; then bcftools annotate -a $HP_RDIR/RNR.bed.gz   $O.gz -c "CHROM,FROM,TO,RNR"  -h <(echo '##INFO=<ID=RNR,Number=1,Type=String,Description="rRNA">') > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/TRN.bed.gz ]   ; then bcftools annotate -a $HP_RDIR/TRN.bed.gz   $O.gz -c "CHROM,FROM,TO,TRN"  -h <(echo '##INFO=<ID=TRN,Number=1,Type=String,Description="tRNA">') > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/DLOOP.bed.gz ] ; then bcftools annotate -a $HP_RDIR/DLOOP.bed.gz $O.gz -c "CHROM,FROM,TO,DLOOP"  -h <(echo '##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="DLOOP">') > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/MCC.bed.gz ]  ; then bcftools annotate -a $HP_RDIR/MCC.bed.gz $O.gz -c "CHROM,FROM,TO,MCC"  -h <(echo '##INFO=<ID=MCC,Number=1,Type=String,Description="New YALE protein_gene_missense_constraint; missense_OEUF">') > $O ; fi

if [ -s $HP_RDIR/HG.vcf.gz ]    ; then bcftools annotate -a $HP_RDIR/HG.vcf.gz    $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/NUMT.vcf.gz ]  ; then bcftools annotate -a $HP_RDIR/NUMT.vcf.gz  $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/dbSNP.vcf.gz ] ; then bcftools annotate -a $HP_RDIR/dbSNP.vcf.gz $O.gz -c "ID"    > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/NONSYN.vcf.gz ]  ; then bcftools annotate -a $HP_RDIR/NONSYN.vcf.gz   $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/STOP.vcf.gz ]  ; then bcftools annotate -a $HP_RDIR/STOP.vcf.gz   $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/MITOMAP.vcf.gz ]  ; then bcftools annotate -a $HP_RDIR/MITOMAP.vcf.gz   $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/MITIMPACT.vcf.gz ]    ; then bcftools annotate -a $HP_RDIR/MITIMPACT.vcf.gz $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/HelixMTdb.vcf.gz ]    ; then bcftools annotate -a $HP_RDIR/HelixMTdb.vcf.gz $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/gnomAD31.vcf.gz  ]    ; then bcftools annotate -a $HP_RDIR/gnomAD31.vcf.gz $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/UKB_dragen.05.vcf.gz ] ; then bcftools annotate -a $HP_RDIR/UKB_dragen.05.vcf.gz $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi
if [ -s $HP_RDIR/MLC.vcf.gz ]  ; then bcftools annotate -a $HP_RDIR/MLC.vcf.gz $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz ; fi

if [ ! -s $O ] ; then zcat $O.gz > $O ; fi
rm -f $O.gz $O.gz.tbi
