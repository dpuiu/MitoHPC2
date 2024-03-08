# MitoHPC2 : Mitochondrial High Performance Caller 2#

    For Calling Mitochondrial Homoplasmies and Heteroplasmies

# Citing #

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

## INPUT ##
 
    * Illumina paired-end reads pre-aligned to a reference: .bam or .cram files 
    * Alignment indeces: .bai or .crai files
    * Alignment indexstats: .idxstats files

## SYSTEM PREREQUISITES ##

    * OS:                UNIX/LINUX 
    * SOFTWARE PACKAGES: git,wget,tar,unzip,make,autoconf,gcc,java,perl,python 
 
## PIPELINE PREREQUISITES ##

    * SOFTWARE PACKAGES: bwa,minimap2,samtools,bedtools,fastp,samblaster,bcftools,htslib,freebayes
    * JAVA JARS:         gatk,mutserve,varscan2,haplogrep,haplocheck
    * HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    $ git clone https://github.com/dpuiu/MitoHPC2.git

### SETUP ENVIRONMENT (important) ###

    $ cd MitoHPC2/scripts
    $ export HP_SDIR=`pwd`                           # set script directory variable 
    $ echo "export HP_SDIR=`pwd`" >> ~/.bashrc       # add variable to your ~/.bashrc so it gets initialized automatically
    $ . ./init.sh                                    # or one of init.{hs38DH,hg19,mm39}.sh  correponding to  a different reference


## RUN PIPELINE  ##
 
   $ $HP_SDIR/run3.sh                                # runs mutect2,varscan,freebayes, combine results

