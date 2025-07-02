#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that downloads and install software prerequisites

##############################################################################################################

#. $HP_SDIR/init.sh

cd $HP_HDIR/
mkdir -p prerequisites/ $HP_BDIR/ $HP_JDIR/

cd prerequisites/

wget -N -c --no-check-certificate https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
wget -N -c https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2

wget -N -c https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
wget -N -c https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
wget -N -c https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
wget -N -c https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
wget -N -c http://opengene.org/fastp/fastp
#wget -N -c https://github.com/OpenGene/fastp/archive/refs/tags/v0.24.1.tar.gz

wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc15/mutserve.zip
wget -N -c https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
wget -N -c  https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar
wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz
wget -N -c https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip
wget -N -c https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip

wget -N -c https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
wget -N -c https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip

