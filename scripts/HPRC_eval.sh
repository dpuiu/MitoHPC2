#!/bin/bash -ex

M=$1    
T=".10"
T=$2
MT="$M$T"


echo "$MT:all"
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -sm 
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -snp -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -ins -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -del -sm

echo
echo "$MT:hom"
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -hom -sm 
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -snp -hom -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -ins -hom -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -del -hom -sm

echo
echo "$MT:het"
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -het -sm 
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -snp -het -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -ins -het -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/HPRC/RefSeq/chrM.39.Hh.norm.vcf $MT.concat.vcf -del -het -sm

