#!/bin/bash -e

M=$1	
T=".10"
T=$2
MT="$M$T"

echo "$MT:all"
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -snp -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -ins -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -del -sm

#perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -snp -af 0.1 -sm
#perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -ins -af 0.1 -sm
#perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -del -af 0.1 -sm

#echo
#echo "$MT:noh"
#perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -snp -noh -sm
#perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -ins -noh -sm
#perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -del -noh -sm

echo
echo "$MT:hom"
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -snp -hom -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -ins -hom -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -del -hom -sm

echo
echo "$MT:het"
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -snp -het -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -ins -het -sm
perl  ~/Homo_sapiens_mito/MitoHPC2/scripts/cmpVcf.pl ~/Homo_sapiens_mito/MitoHPC2/examples/Sim/RefSeq/chrM.30.Hh.norm.vcf $MT.concat.vcf -del -het -sm


