# MitoHPC2 : Mitochondrial High Performance Caller v2 #

For Calling Mitochondrial Homoplasmies and Heteroplasmies

# Citing #

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

# INSTALL/RUN ##
    
Check https://github.com/dpuiu/MitoHPC/blob/main/README.md

# RUN MUTIPLE SNV CALLERS, MERGE RESULTS #
     
    * Runs mutect2,varscan,freebayes    
    * Merges the calls; an SNV must be called by at least 2 of the metods
    
    $ cp $HP_SDIR/init3.sh .
    $ . ./init3.sh
    $ $HP_SDIR/run3.sh        

    $ ls $HP_ODIR/merge3.*
