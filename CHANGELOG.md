# 2025/08/05 #
* Added deepvariant as one of the SNV callers

# 2025/07/23 #
* Added clair3 as one of the SNV callers

# 2025/01/08 #
* MitoHPC2/examples/Sim/Illumina/ <- MitoHPC/examples2   : 30 WGS simulated Illumina samples
* MitoHPC2/examples/HPRC/{Illumina,HiFi}                 : 40 WGS HPRC Illumina/HiFi samples; direct MT filtering using "samtools view s3://"

# 2025/01/07 #
* RefSeq/{CDS,RNR,TRN}.bed.gz*         : switched to Hugo gene names (MT-*)
* scripts/install_sysprerequisites.sh  : repaced java packages with "openjdk-17-jdk"  and "java-17-openjdk"
* added Dockerfile                     : for creating docker images;

# 2024/12/03 #
* scripts/uniqVcf.pl                   : converted REF/ALT to upper case
* scripts/maxVcf.pl                    : converted REF/ALT to upper case

# 2024/10/14 #
* scripts/install_sysprerequisites.sh  : added libcurl4-gnutls-dev libssl-dev zlib1g-dev
* scripts/install_prerequisites.sh     : software version update
                                         minimap2-2.26   => 2.28
                                         samtools-1.16   => 1.21 --with-curl
                                         htslib-1.16     => 1.16 --with-curl
                                         bcftools-1.16   => 1.21
                                         bedtools-2.30.0 => 2.31.1
                                         gatk-4.3.0.0    => 4.6.0.0
                                         freebayes-1.3.6 => 1.3.8
                                         delly_v1.2.9    => 1.3.1
* scripts/ls2in.pl                     : added "-word" parameter : word only sample names
* examples/HPRC/Illumina/in.url        : URLs of 40 HPRC Illumina samples
* examples/HPRC/HiFi/in.url            : URLs of        40 HPRC HiFi     samples

# 2024/10/11 #
* scripts/install_prerequisites.sh     : updated gridss, added Rscript , delly, plink2 (uncomment)
* scripts/checkInstall.sh              : added Rscript , delly, plink2 (uncomment)
* scripts/init.hifi.sh                 : lowered min coverage to HP_DP=10
* scripts/run.sh                       : checked if init script has been run
* scripts/filter.sh                    : multiple updates
                                         updated st.pl call
                                         added varscan2Vcf.pl call
                                         allowed for other SV callers(delly)
                                         updated circFasta.sh call
                                         updated rotateFasta.h call
                                         removed "-H A" from  bcftools consensus
* scripts/filterHiFi.sh                : multiple updates
                                         miniamp2 : added --eqx
                                         support for freebayes SNV caller
* scripts/cat2concat.pl                : added FORMAT/AD
* scripts/concat2merge.pl              : updated GT
* scripts/st.pl                        : added -sample parameter
* added scripts/varscan2Vcf.pl         : process varscan output
* scripts/filterVcf.pl                 : GT:0/1;1/0;1; REF,ALT converted to uppercase if filtering thold>0
* scripts/getSummary.sh                : summarize SV's

# 2024/09/15 #
* scripts/filterHiFi.sh                : added "-x map-hifi" option to minimap2
* scripts/init.hifi.sh                 : initialized HP_MTLEN

# 2024/09/12 # 
* scripts/README.md                    : updated README.med file
* scripts/run.hifi.sh                  : added to process HiFi reads
* scripts/fixbcftoolsVcf.pl            : added the FORMAT/AD tag => can process HiFi reads
* scripts/fixfreebayesVcf.pl           : added the FORMAT/AD tag
* scripts/fixmutserveVcf.pl            : added the FORMAT/AD tag

# 2024/08/02 #
* added scripts/filter1.sh             : re-filter the original mutect2 calls and add AD tag to mutect2.00 file

# 2024/08/01 #
* scripts/fixsnpPos.pl                 : using lower cases for 2nd iteration "flipped" SNVs
* scripts/differenceVcf.pl             : allow for lower eq upper cases

# 2024/07/31 #
* scripts/filter2.sh                   : file path update

# 2024/07/30 #
* scripts/annotateVcf.sh               : added dbSNP INFO
* scripts/filter.sh                    : added dbSNP INFO
* scripts/filterVcf.pl                 : added the FORMAT/AD tag
* scripts/fixmutect2Vcf.pl             : added the FORMAT/AD tag
* scripts/*.vcf                        : added the FORMAT/AD tag
* scripts/filter2.sh                   : added new script to filter the original mutect2.mutect2 calls => mutect2.mutect2.00.orig

# 2024/05/23 #
* RefSeq/{CDS,COMPLEX,TRN}.bed.gz.*    : annotate the overlapping features using both names
* removed RefSeq/STOP.vcf.gz*          : STOP gain annotation present in MLC.vcf.gz*

# 2024/03/12 #
* added scripts/filterHiFi.sh          : process HiFi reads 

* added scripts/init3.sh               : setup parameters associated with running 3 SNV callers
* added scripts/run3.sh                : run 3 SNV callers (default: mutect2, varscan, freebayes)
* added scripts/filter3.sh             : filter the SNVs called by at least 2 out of 3 SNV callers
* added scripts/merge3Vcf.sh           : merge results (AD,AF,....)
* updated scripts/annotateVcf.sh       : annotate SNvs using the *.bed.gz and *.vcf.gz files available in the RefSeq/
                                         Existing:   
                                            dbSNP.vcf.gz
                                            HG.vcf.gz
                                            MITIMPACT.vcf.gz
                                            MLC.vcf.gz 
                                            NONSYN.vcf.gz
                                            NUMT.vcf.gz

                                            CDS.bed.gz
                                            COMPLEX.bed.gz
                                            DLOOP.bed.gz
                                            Homopolymer.bed.gz
                                            Hotspot.bed.gz
                                            Hypervariable.bed.gz
                                            MCC.bed.gz
                                            RNR.bed.gz
                                            TRN.bed.gz

                                         New:      
                                           gnomAD31.vcf.gz
                                           MITOMAP.vcf.gz
                                           UKB_mutect2_Nature_20230816.vcf.gz
                                           HelixMTdb.vcf.gz


# 2025/01/07 #
* RefSeq/{CDS,RNR,TRN}.bed.gz*         : switched to Hugo gene names (MT-*)  
* scripts/install_sysprerequisites.sh  : repaced java packages with "openjdk-17-jdk"  and "java-17-openjdk"
* added Dockerfile                     : for creating docker images; 

# 2025/01/08 #
* Translating multiple PERL scripts to PYTHON
* MitoHPC2/examples/Sim/Illumina/ <- MitoHPC/examples2   : 30 WGS simulated Illumina samples
* MitoHPC2/examples/HPRC/{Illumina,HiFi}                 : 40 WGS HPRC Illumina/HiFi samples; direct MT filtering using "samtools view s3://"
