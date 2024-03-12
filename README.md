# MitoHPC2 : Mitochondrial High Performance Caller v2 #

For Calling Mitochondrial Homoplasmies and Heteroplasmies

# Citing #

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

# INSTALL/RUN ##
    
Check https://github.com/dpuiu/MitoHPC/blob/main/README.md

The pipleine has been updated so that all MitoHPC2/RefSeq/*.{vcf,bed}.gz files are used for annotation. 
If you have any custom annotation files you would like to use, just copy them to MitoHPC2/RefSeq. 
Make sure the VCF/BED files are gzipped and indexed.

In addition, one can run multiple(3) SNV callers and merge the results.
     
    * Example: Run mutect2,varscan,freebayes    
    * Merges the calls; filter the SNV called by at least 2 of the metods
    
    $ cp $HP_SDIR/init3.sh .
    $ . ./init3.sh
    $ $HP_SDIR/run3.sh | tee run3.all.sh | bash       

    # output
    $ ls $HP_ODIR/merge3.*
