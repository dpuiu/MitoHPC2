# MitoHPC2 : Mitochondrial High Performance Caller v2 #

Pipeline for Calling Mitochondrial Homoplasmies and Heteroplasmies

# Citing #

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, 
[NAR 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/)

# Prerequisites # 

[bwa](https://github.com/lh3/bwa/releases) 0.7.17, 
[minimap2](https://github.com/lh3/minimap2/releases) 2.28, 
[htslib](https://github.com/samtools/htslib/releases) 1.21, 
[samtools](https://github.com/samtools/samtools/releases) 1.21, 
[bcftools](https://github.com/samtools/bcftools/releases) 1.21, 
[bedtools](https://github.com/arq5x/bedtools2/releases) 2.31.1, 
[fastp](http://opengene.org/fastp/fastp) 0.24.0, 
[samblaster](https://github.com/GregoryFaust/samblaster/releases) 0.1.26, 
[gatk Mutect2](https://github.com/broadinstitute/gatk/releases/) 4.6.0.0, 
[mutserve](https://github.com/seppinho/mutserve/releases) 2.0.0-rc15, 
[freebayes](https://github.com/freebayes/freebayes/releases) 1.3.6, 
[VarScan](https://github.com/dkoboldt/varscan/releases) 2.4.6, 
[Clair3](https://github.com/HKU-BAL/Clair3/releases) 1.1.1, 
[ClairS-TO](https://github.com/HKU-BAL/ClairS-TO) 0.4.2,
[deepvariant](https://github.com/google/deepvariant) 1.10.0-beta, 
[deepsomatic](https://github.com/google/deepsomatic/) 1.9.0,
[gridss](https://github.com/PapenfussLab/gridss/releases) 2.13.2, 
[haplogrep](https://github.com/seppinho/haplogrep-cmd/releases) v2.4.0, 
[haplocheck](https://github.com/genepi/haplocheck/releases) 1.3.3, 
[R](https://cran.r-project.org/src/base/) 4.3.0, 
[plink2](https://www.cog-genomics.org/plink/2.0/) 2.0 

# INSTALL # 
    
     # install git
     sudo apt-get -y update                            
     sudo apt-get install -y git                          

     # clone and init MitoHPC2
     git clone https://github.com/dpuiu/MitoHPC2.git
     cd MitoHPC2/scripts
     export HP_SDIR=`pwd -P`
     . ./init.sh

     # install system prerequisites (if admin or checkInstall.sh below fails)
     sudo ./install_sysprerequisites.sh                

     # install prerequisites and check install
     ./install_prerequisites.sh
     ./checkInstall.sh

For additional information check [MitoHPC README](https://github.com/dpuiu/MitoHPC/blob/main/README.md) 

## CUSTOM ANNOTATION ## 

The pipleine has been updated so that all MitoHPC2/RefSeq/*.{vcf,bed}.gz files are used for annotation. 
If you have any custom annotation files you would like to use, just copy them to MitoHPC2/RefSeq/
Make sure the VCF/BED files are gzipped and indexed.

# RUN #

Check [MitoHPC README](https://github.com/dpuiu/MitoHPC/blob/main/README.md) first !!!

Copies of the  Clair3, ClairS-To, DeepVarian, DeepSomatic SIF files can be found at [ftp](ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/MitoHPC2/bin/)

Copies of the HPRC samples BAM alignments can be found at [ftp](ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/MitoHPC2/bin/examples/HPRC/)

## SINGLE SNV CALLER (Illumina)

    # copy init file to work directory
    cp $HP_SDIR/init.sh .
  
    # edit SNV caller if needed
    nano ./init.sh 
      ...
      export HP_M=mutect2          # mutect2,mutserve,freebayes,varscan

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

    . ./init3.sh

    # run
    $HP_SDIR/run3.sh | tee run.all.sh | bash       

    # output
    ls $HP_ODIR/merge3.*

## SNV CALLER IINIT PacBio HiFi

    # copy init file to work directory
    cp $HP_SDIR/init.hifi.sh .

    # edit init.lr.sh; set the SNV caller
    cat init.hifi.sh
      ...
      export HP_M=deepsomatic                       # varscan,bcftools,clair3,clairs-to,deepvariant,deepsomatic
      export HP_PLATFORM="hifi_revio"		    # used by clairs-to
      export HP_MODELTYPE="PACBIO_TUMOR_ONLY"       # used by deepsomatic

    # init
    . ./init.hifi.sh

## SNV CALLER INIT ONT

    # copy init file to work directory
    cp $HP_SDIR/init.ont.sh .

    # edit init.lr.sh; set the SNV caller
    cat init.ont.sh
      ...
      export HP_M=clairs-to                         # varscan,bcftools,clair3,clairs-to,deepvariant,deepsomatic
      export HP_PLATFORM="ont_r10_dorado_sup_5khz"  # used by clairs-to ; ont_r10_dorado_sup_4khz, ont_r10_dorado_hac_4khz, ont_r10_dorado_sup_5khz, ont_r10_dorado_sup_5khz_ss, ont_r10_dorado_sup_5khz_ssrs
      export HP_MODELTYPE="ONT_TUMOR_ONLY"	    # used by deepsomatic

    # init
    . ./init.ont.sh

## SNV CALLER LONG READS

    # run
    $HP_SDIR/run.lr.sh | tee run.all.sh | bash

    # check output
    ls $HP_ODIR/.*

    # evaluate results(compare to Illumina mutect2 T=10)             
    $HP_SDIR/eval.sh Illumina/out/mutect2.10.concat.vcf $HP_M.10.concat.vcf | uniq.pl | column -t > $HP_M.10.eval

# Examples #

## 30 Simulated samples ##

## 39 HPRC samples ##

### HiFi10 ###

* Note: 
** 1st number: short read heteroplasmy thold  
** 2nd: long read heteroplasmy thold  

* clairs-to: 
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.10.10.eval), 
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.05.10.eval), 
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.10.05.eval), 
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.05.05.eval) 

* deepsomatic:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.10.10.eval), 
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.05.10.eval), 
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.10.05.eval), 
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.05.05.eval)

### ONT10 ###

* clairs-to:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.10.10.eval), 
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.05.10.eval), 
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.10.05.eval), 
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.05.05.eval)

* deepsomatic:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.10.10.eval), 
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.05.10.eval), 
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.10.05.eval), 
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.05.05.eval)

* clairs-to intersect deepsomatic:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to_deepsomatic.10.10.eval)

