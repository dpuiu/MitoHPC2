#!/usr/bin/env bash
set -ux # do not add -e

##############################################################################################################

# Program that downloads and installs software prerequisites and genome reference
#  -f : opt; force reinstall

##############################################################################################################

if [ -z $HP_SDIR ] ; then echo "Variable HP_SDIR not defined. Make sure you followed the SETUP ENVIRONMENT instructions" ;  exit 0 ; fi
if [ -z $HP_HDIR ] ; then echo "Variable HP_HDIR not defined. Make sure you followed the SETUP ENVIRONMENT instructions" ;  exit 0 ; fi

#. $HP_SDIR/init.sh
cd $HP_HDIR
mkdir -p prerequisites/ $HP_BDIR/ $HP_JDIR/ $HP_RDIR/
cd prerequisites/

#compile using multiple threads

##############################################################################################################

which bwa
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c --no-check-certificate https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
  if [ ! -s $HP_BDIR/bwa ] ; then
    tar -xjvf bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make  CFLAGS="-g -Wall -Wno-unused-function -O2 -fcommon"  # compiling using gcc v10.+ fails unless "-fcommon" is added
    mv bwa $HP_BDIR/
    cd -
  fi
fi

which minimap2
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2
  tar -xjvf minimap2-2.28.tar.bz2 
  cd minimap2-2.28/
  make ;  mv minimap2 $HP_BDIR
  cd -
fi

##############################################################################################################

which htsfile
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
  if [ ! -s $HP_BDIR/tabix ] ; then
    tar -xjvf htslib-1.21.tar.bz2
    cd htslib-1.21
    ./configure --prefix=$HP_HDIR/ --with-curl # --disable-bz2
    make ; make install
    cd -
  fi
fi

which samtools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
  if [ ! -s $HP_BDIR/samtools ] ; then
    tar -xjvf samtools-1.21.tar.bz2
    cd samtools-1.21
    ./configure --prefix=$HP_HDIR/ --with-curl # --without-curses --disable-bz2
    make ;  make install
    cd -
  fi
fi

which bcftools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
  if [ ! -s $HP_BDIR/bcftools ] ; then
    tar -xjvf  bcftools-1.21.tar.bz2
    cd bcftools-1.21
    ./configure --prefix=$HP_HDIR/ # --disable-bz2
    make  ; make install
    cd -
  fi
fi


which samblaster
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
  if [ ! -s $HP_BDIR/samblaster ] ; then
    tar -xzvf samblaster-v.0.1.26.tar.gz
    cd samblaster-v.0.1.26
    make ; mv samblaster $HP_BDIR/
    cd -
  fi
fi

which bedtools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
  if [ ! -s $HP_BDIR/bedtools ] ; then
    tar -xzvf bedtools-2.31.1.tar.gz
    cd bedtools2/
    make install prefix=$HP_HDIR/
    cd -
  fi
fi

which fastp
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c http://opengene.org/fastp/fastp
  mv fastp $HP_BDIR/
  chmod a+x $HP_BDIR/fastp
  #wget -N -c https://github.com/OpenGene/fastp/archive/refs/tags/v0.24.1.tar.gz
  #tar -xzvf v0.24.1.tar.gz 
  #cd fastp-0.24.1/
  #make
  #make install prefix=$HP_HDIR/
  #cd -
fi

#########################################################################################

#if [ ! -s $HP_JDIR/gatk.jar ] ; then # 2023/04/26
if [[ ! -s $HP_JDIR/gatk.jar || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
  unzip -o gatk-4.6.0.0.zip
  mv gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar $HP_JDIR/gatk.jar
  mv gatk-4.6.0.0/gatk $HP_BDIR/
fi

#if [ ! -s $HP_JDIR/mutserve.jar ] ; then  # 2023/04/26
if [[ ! -s $HP_JDIR/mutserve.jar || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc15/mutserve.zip
  unzip -o mutserve.zip
  mv mutserve.jar $HP_JDIR/
fi

which freebayes
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
 
  wget -N -c https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
  gunzip freebayes-1.3.6-linux-amd64-static.gz  -c >  $HP_BDIR/freebayes
  chmod a+x $HP_BDIR/freebayes
fi

if [[ ! -s $HP_JDIR/VarScan.jar || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c  https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar
  mv VarScan.v2.4.6.jar $HP_JDIR/VarScan.jar
fi

which Rscript
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
  tar -xzvf R-4.3.0.tar.gz
  cd R-4.3.0
  ./configure --prefix=$HP_HDIR
  make
  make install
fi

which gridss
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz
  tar -xzvf gridss-2.13.2.tar.gz
  mv gridss $HP_BDIR/
  mv gridss-2.13.2-gridss-jar-with-dependencies.jar $HP_JDIR/gridss.jar
fi

#which delly
#if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  #wget -N -c https://github.com/dellytools/delly/releases/download/v1.3.1/delly_v1.3.1_linux_x86_64bit
  #mv delly_v1.3.1_linux_x86_64bit $HP_BDIR/delly
#fi

which plink2
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
  unzip plink2_linux_x86_64_latest.zip
  mv plink2 $HP_BDIR/
fi

#test -s $HP_BDIR/clair3_latest.sif
#test -d ~/clair3_sandbox
#if [[ $? != 0 ]] ; then
#  singularity pull docker://hkubal/clair3:latest
#  mv -i clair3_latest.sif $HP_BDIR
#  singularity build --sandbox ~/clair3_sandbox $HP_BDIR/clair3_latest.sif
#fi

#test -d ~/deepvariant_sandbox
#if [[ $? != 0 ]] ; then
#  singularity pull docker://google/deepvariant:latest
#  mv -i deepvariant_latest.sif $HP_BDIR
#  singularity build --sandbox ~/deepvariant_sandbox $HP_BDIR/deepvariant_latest.sif
#fi

####################################################################################

#if [ ! -s $HP_JDIR/haplogrep.jar ] ; then # 2023/04/26
if [[ ! -s $HP_JDIR/haplogrep.jar || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip
  unzip -o haplogrep.zip
  mv haplogrep.jar $HP_JDIR/
fi

#if [ ! -s $HP_JDIR/haplocheck.jar ] ; then # 2023/04/26
if [[ ! -s $HP_JDIR/haplocheck.jar || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip
  unzip -o haplocheck.zip
  mv haplocheck.jar $HP_JDIR/
fi

#####################################################################################

#if [ ! -s $HP_RDIR/$HP_RNAME.fa ] ; then # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_RNAME.fa.fai || $# == 1 && $1 == "-f" ]] ; then
  wget -qO- $HP_RURL | zcat -f > $HP_RDIR/$HP_RNAME.fa
  #wget -q $HP_RURL -O $HP_RDIR/$HP_RNAME.fa  # wget on fedora does not download ftp links
  #curl -L $HP_RURL -o $HP_RDIR/$HP_RNAME.fa
  samtools faidx $HP_RDIR/$HP_RNAME.fa
fi

#if [ ! -s $HP_RDIR/$HP_MT.fa ] ; then  # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_MT.dict || $# == 1 && $1 == "-f"  ]] ; then
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RMT > $HP_RDIR/$HP_MT.fa
  samtools faidx $HP_RDIR/$HP_MT.fa
  rm $HP_RDIR/$HP_MT.dict
  java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $HP_RDIR/$HP_MT.fa --OUTPUT $HP_RDIR/$HP_MT.dict
fi

#if [ ! -s $HP_RDIR/$HP_NUMT.fa ] ; then # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_NUMT.bwt || $# == 1 && $1 == "-f"  ]] ; then
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RNUMT > $HP_RDIR/$HP_NUMT.fa
  bwa index $HP_RDIR/$HP_NUMT.fa -p $HP_RDIR/$HP_NUMT
fi

#if [ ! -s $HP_RDIR/$HP_MTC.fa ] ; then  # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_MTC.dict || $# == 1 && $1 == "-f" ]] ; then
  circFasta.sh $HP_MT $HP_RDIR/$HP_MT $HP_E $HP_RDIR/$HP_MTC
fi

#if [ ! -s $HP_RDIR/$HP_MTR.fa ] ; then # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_MTR.dict || $# == 1 && $1 == "-f" ]] ; then
  rotateFasta.sh $HP_MT $HP_RDIR/$HP_MT $HP_E $HP_RDIR/$HP_MTR
fi
