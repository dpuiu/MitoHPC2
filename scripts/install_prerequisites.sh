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
mkdir -p prerequisites/ $HP_BDIR/ $HP_RDIR/
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
    cp bwa $HP_BDIR/
    cd -
  fi
fi

which minimap2
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2
  tar -xjvf minimap2-2.28.tar.bz2 
  cd minimap2-2.28/
  make ;  cp minimap2 $HP_BDIR
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
    ./configure --prefix=$HP_HDIR/ --with-curl  --without-curses # --disable-bz2
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
    make ; cp samblaster $HP_BDIR/
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
  cp fastp $HP_BDIR/
  chmod a+x $HP_BDIR/fastp
  #wget -N -c https://github.com/OpenGene/fastp/archive/refs/tags/v0.24.1.tar.gz
  #tar -xzvf v0.24.1.tar.gz 
  #cd fastp-0.24.1/
  #make
  #make install prefix=$HP_HDIR/
  #cd -
fi
#########################################################################################

which gatk
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
  unzip -o gatk-4.6.0.0.zip
  cp gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar $HP_BDIR/gatk.jar
  cp gatk-4.6.0.0/gatk $HP_BDIR/
fi


which mutserve
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc15/mutserve.zip
  unzip -o mutserve.zip
  cp mutserve mutserve.jar $HP_BDIR/
fi

which freebayes
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then 
  wget -N -c https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
  gunzip freebayes-1.3.6-linux-amd64-static.gz  -c >  $HP_BDIR/freebayes
  chmod a+x $HP_BDIR/freebayes
fi

which varscan
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c  https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar
  cp VarScan.v2.4.6.jar $HP_BDIR/VarScan.jar
  echo "#!/usr/bin/env bash" > $HP_BDIR/varscan
  echo "java -jar $HP_BDIR/VarScan.jar \$@" >> $HP_BDIR/varscan
  chmod a+x $HP_BDIR/varscan
fi

#which Rscript
#if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
#  wget -N -c https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
#  tar -xzvf R-4.3.0.tar.gz
#  cd R-4.3.0
#  ./configure --prefix=$HP_HDIR --with-readline=no
#  make
#  make install
#fi

which gridss
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz
  tar -xzvf gridss-2.13.2.tar.gz
  cp gridss $HP_BDIR/
  cp gridss-2.13.2-gridss-jar-with-dependencies.jar $HP_BDIR/gridss.jar
fi

#which delly
#if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  #wget -N -c https://github.com/dellytools/delly/releases/download/v1.3.1/delly_v1.3.1_linux_x86_64bit
  #cp delly_v1.3.1_linux_x86_64bit $HP_BDIR/delly
#fi

which plink2
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
  unzip plink2_linux_x86_64_latest.zip
  cp plink2 $HP_BDIR/
fi
####################################################################################

which haplogrep
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip
  unzip -o haplogrep.zip
  cp haplogrep.jar haplogrep $HP_BDIR
fi

which haplocheck
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip
  unzip -o haplocheck.zip
  cp haplocheck.jar haplocheck $HP_BDIR/
fi

