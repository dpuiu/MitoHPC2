#!/usr/bin/env bash
#set -e
set -euo pipefail

##############################################################################################################

# Program that checks if all dependencies are installed

##############################################################################################################

#. $HP_SDIR/init.sh

echo "########################"  > checkInstall.log
echo -n "DATE:"     >  checkInstall.log
date             >> checkInstall.log
echo "########################"  >> checkInstall.log
echo "EXECUTABLES PATHS:"  >> checkInstall.log
#test executables
which perl       >> checkInstall.log	        #usually available on Linux
which gcc        >> checkInstall.log	        #module load gcc
which java       >> checkInstall.log	        #module load java
which python     >> checkInstall.log		#module load python
#which R

which bwa        >> checkInstall.log	        #module load bwa
which samtools   >> checkInstall.log         #module load samtools
which bedtools   >> checkInstall.log         #module load bedtools
which fastp      >> checkInstall.log
which samblaster >> checkInstall.log
which bcftools   >> checkInstall.log
which tabix      >> checkInstall.log
which freebayes  >> checkInstall.log
which minimap2   >> checkInstall.log
which plink2     >> checkInstall.log
#which gridss
#which delly

which gatk       >> checkInstall.log
which mutserve   >> checkInstall.log
which haplogrep  >> checkInstall.log
which haplocheck >> checkInstall.log
which varscan    >> checkInstall.log

######################################################

echo "########################"  >> checkInstall.log
echo "EXECUTABLES VERSIONS:"  >> checkInstall.log

samtools --version   | head -1 >> checkInstall.log
bcftools --version   | head -1 >> checkInstall.log
bedtools --version   | head -1 >> checkInstall.log
plink2 --version     | head -1 >> checkInstall.log
fastp --version 2>&1 | head -1  >> checkInstall.log
echo -n "bwa "  >> checkInstall.log; (bwa 2>&1 | grep Version -m 1 >> checkInstall.log) || true
samblaster --version 2>&1 | head -1 >> checkInstall.log
echo -n "tabix " >> checkInstall.log ; (tabix 2>&1 |  grep -v ^$ | head -1 >> checkInstall.log) || true
perl --version | grep -v ^$ | head -1  >> checkInstall.log
java -version | head -1               >> checkInstall.log
#gridss
#delly


echo "########################"  >> checkInstall.log
echo "ENV VARIABLES:"  >> checkInstall.log

echo "HP_HDIR=$HP_HDIR" >> checkInstall.log
echo "HP_BDIR=$HP_BDIR" >> checkInstall.log
echo "HP_SDIR=$HP_SDIR" >> checkInstall.log
echo "HP_RDIR=$HP_RDIR" >> checkInstall.log

cat checkInstall.log

echo Success!
