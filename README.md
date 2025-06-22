# MitoHPC2 : Mitochondrial High Performance Caller v2 #

For Calling Mitochondrial Homoplasmies and Heteroplasmies

# Citing #

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

# INSTALL # 
    
     # make sure git is installed
     $ which git

     # install system prerequisites(if admin and git not installed)
     $ sudo apt-get -y update                            
     $ sudo apt-get install git                          

     # clone and init 
     $ git clone https://github.com/dpuiu/MitoHPC2.git
     $ cd MitoHPC2/scripts
     $ export HP_SDIR=`pwd`
     $ . ./init.sh

     # install system prerequisites(if admin or checkInstall.sh below fails)
     $ sudo ./install_sysprerequisites.sh                

     # install prerequisites and check install
     $ ./install_prerequisites.sh
     $ ./checkInstall.sh

For additional information check https://github.com/dpuiu/MitoHPC/blob/main/README.md

## CUSTOM ANNOTATION ## 

The pipleine has been updated so that all MitoHPC2/RefSeq/*.{vcf,bed}.gz files are used for annotation. 
If you have any custom annotation files you would like to use, just copy them to MitoHPC2/RefSeq/
Make sure the VCF/BED files are gzipped and indexed.

# RUN #

Check https://github.com/dpuiu/MitoHPC/blob/main/README.md first !!!

## SINGLE SNV CALLER (Illumina)

    # copy init file to work directory
    $ cp $HP_SDIR/init.sh .

    # init
    $ . ./init.sh

    # run
    $ $HP_SDIR/run.sh | tee run.all.sh | bash

    # check output
    $ ls $HP_ODIR/mutect2.*

## MULTIPLE SNV CALLERS (Illumina)

In addition, one can run multiple(3) SNV callers and merge the results. 
Only the SNV called by at least 2 the metods make it into the final/merged set.

    # init     
    $ cp $HP_SDIR/init3.sh .
    $ cat init3.sh
      ...
      export HP_M1=mutect2   
      export HP_M2=varscan
      export HP_M3=freebayes
      export HP_M=merge3
      ...

    $ . ./init3.sh

    # run
    $ $HP_SDIR/run3.sh | tee run3.all.sh | bash       

    # output
    $ ls $HP_ODIR/merge3.*

## SINGLE SNV CALLER (PacBio HiFi/ONT)

    # copy init file to work directory
    $ cp $HP_SDIR/init.lr.sh .

    # edit init.lr.sh; set the SNV caller
    $ nano init.lr.sh
      HP_M=bcftools  # or
      HP_M=varscan

    # init
    $ . ./init.lr.sh

    # run
    $ $HP_SDIR/run.lr.sh | tee run.lr.all.sh | bash       

    # check output
    $ ls $HP_ODIR/.*

# Examples #

## 40 HPRC samples ##

Illumina vs PacBio HiFi 

### Illumina ###

    $ cd examples/HPRC/Illumina/
    $ cat in.url | ls2in.pl -word -out $PWD/out | sort | head -n 5 > in.txt
    $ cp $HP_SDIR/init.sh .
    $ . ./init.sh
    $ $HP_SDIR/run.sh | tee run.all.sh | bash
    $ ls $HP_ODIR/

### PacBio HiFi ###

    $ cd examples/HPRC/HiFi/
    $ cat in.url | ls2in.pl -word -out $PWD/out | sort | head -n 5 > in.txt
    $ cp $HP_SDIR/init.lr.sh .
    $ . ./init.lr.sh
    $ $HP_SDIR/run.lr.sh | tee run.lr.all.sh | bash
    $ ls $HP_ODIR/         

## 30 Simulated samples ##

Illumina vs ONT R10

