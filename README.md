# MitoHPC2 : Mitochondrial High Performance Caller v2 #

Pipeline for Calling Mitochondrial Homoplasmies and Heteroplasmies

# Citing #

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

# Prerequisites # 

[bwa](https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2) 
[minimap2](https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2) 
[htslib](https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2) 
[samtools](https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2) 
[bcftools](https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2) 
[bedtools](https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz) 
[samblaster](https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz) 
[gatk](https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip) 
[mutserve](https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc15/mutserve.zip) 
[freebayes](https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz) 
[VarScan](https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar) 
[clair3](docker://hkubal/clair3:latest) 
[gridss](https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz) 
[haplogrep](https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip) 
[haplocheck](https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip) 
[fastp](http://opengene.org/fastp/fastp) 
[R](https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz) 
[plink2](https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip) 

# INSTALL # 
    
     # install git
     sudo apt-get -y update                            
     sudo apt-get install -y git                          

     # clone and init MitoHPC2
     git clone https://github.com/dpuiu/MitoHPC2.git
     cd MitoHPC2/scripts
     export HP_SDIR=`pwd`
     . ./init.sh

     # install system prerequisites (if admin or checkInstall.sh below fails)
     sudo ./install_sysprerequisites.sh                

     # install prerequisites and check install
     ./install_prerequisites.sh
     ./checkInstall.sh

For additional information check https://github.com/dpuiu/MitoHPC/blob/main/README.md

## CUSTOM ANNOTATION ## 

The pipleine has been updated so that all MitoHPC2/RefSeq/*.{vcf,bed}.gz files are used for annotation. 
If you have any custom annotation files you would like to use, just copy them to MitoHPC2/RefSeq/
Make sure the VCF/BED files are gzipped and indexed.

# RUN #

Check https://github.com/dpuiu/MitoHPC/blob/main/README.md first !!!

## SINGLE SNV CALLER (Illumina)

    # copy init file to work directory
    cp $HP_SDIR/init.sh .

    # init
    . ./init.sh

    # run
    $HP_SDIR/run.sh | tee run.all.sh | bash

    # check output
    ls $HP_ODIR/mutect2.*

## MULTIPLE SNV CALLERS (Illumina)

In addition, one can run multiple(3) SNV callers and merge the results. 
Only the SNV called by at least 2 the metods make it into the final/merged set.

    # init     
    cp $HP_SDIR/init3.sh .
    cat init3.sh
      ...
      export HP_M1=mutect2   
      export HP_M2=varscan
      export HP_M3=freebayes
      export HP_M=merge3
      ...

    . ./init3.sh

    # run
    $HP_SDIR/run3.sh | tee run3.all.sh | bash       

    # output
    ls $HP_ODIR/merge3.*

## SINGLE SNV CALLER (PacBio HiFi/ONT)

    # copy init file to work directory
    cp $HP_SDIR/init.lr.sh .

    # edit init.lr.sh; set the SNV caller
    nano init.lr.sh
      HP_M=bcftools  # or
      HP_M=varscan

    # init
    . ./init.lr.sh

    # run
    $HP_SDIR/run.lr.sh | tee run.lr.all.sh | bash       

    # check output
    ls $HP_ODIR/.*

# Examples #

## 40 HPRC samples ##

Illumina vs PacBio HiFi vs ONT R10

### Illumina ###

    cd examples/HPRC/Illumina/
    cat in.url | ls2in.pl -word -out $PWD/out | sort | head -n 5 > in.txt
    cp $HP_SDIR/init.sh .
    . ./init.sh
    $HP_SDIR/run.sh | tee run.all.sh | bash
    ls $HP_ODIR/

### PacBio HiFi ###

    cd examples/HPRC/HiFi/
    cat in.url | ls2in.pl -word -out $PWD/out | sort | head -n 5 > in.txt
    cp $HP_SDIR/init.lr.sh .
    . ./init.lr.sh
    $HP_SDIR/run.lr.sh | tee run.lr.all.sh | bash
    ls $HP_ODIR/         

## 30 Simulated samples ##


