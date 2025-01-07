# 2024/05/23 #

* updated RefSeq/{CDS,COMPLEX,TRN}.bed.gz.* : annotate the overlapping features using both names
* remove  RefSeq/STOP.vcf.gz*               : STOP gain annotation present in MLC.vcf.gz*

# 2024/07/30 #

* scripts/annotateVcf.sh               : added dbSNP INFO
* scripts/filter.sh                    : added dbSNP INFO
* scripts/filterVcf.pl                 : added the FORMAT/AD tag
* scripts/fixmutect2Vcf.pl             : added the FORMAT/AD tag
* scripts/*.vcf                        : added the FORMAT/AD tag

* scripts/filter2.sh                   : added new script to filter the original mutect2.mutect2 calls => mutect2.mutect2.00.orig

# 2024/07/31 #

* scripts/filter2.sh	               : file path update

# 2024/08/01 #

* scripts/fixsnpPos.pl	               : using lower cases for 2nd iteration "flipped" SNVs 
* scripts/differenceVcf.pl             : allow for lower eq upper cases 

# 2024/08/02 #

* scripts/filter1.sh	               : added new script to re-filter the original mutect2 calls and add AD tag to mutect2.00 file

# 2024/09/12 # 

* scripts/README.md                    : updated README.med file
* scripts/run.hifi.sh                  : added to process HiFi reads
* scripts/fixbcftoolsVcf.pl            : added the FORMAT/AD tag => can process HiFi reads
* scripts/fixfreebayesVcf.pl           : added the FORMAT/AD tag
* scripts/fixmutserveVcf.pl            : added the FORMAT/AD tag

# 2024/09/15 #

* scripts/filterHiFi.sh                : minimap2 added "-x map-hifi"
* scripts/init.hifi.sh	               : initialized HP_MTLEN

# 2024/10/11 #

* scripts/install_prerequisites.sh     : updated gridss, added Rscript , delly, plink2 (uncomment)
* scripts/checkInstall.sh              : added Rscript , delly, plink2 (uncomment)
* scripts/init.hifi.sh                 : lower min coverage to HP_DP=10
* scripts/run.sh                       : chech if init script has been run
* scripts/filter.sh                    : multiple updates
                                         updated st.pl call
                                         added varscan2Vcf.pl call
                                         allow for other SV callers(delly)
                                         update circFasta.sh call
                                         update rotateFasta.h call
                                         bcftools consensus : removed "-H A"
* scripts/filterHiFi.sh                : multiple updates
                                         miniamp2 : added --eqx
                                         support for freebayes SNV caller
* scripts/cat2concat.pl                : added FORMAT/AD
* scripts/concat2merge.pl              : updated GT
* scripts/st.pl                        : added -sample parameter
* scripts/varscan2Vcf.pl               : new script
* scripts/filterVcf.pl                 : GT:0/1;1/0;1; REF,ALT converted to uppercase if filtering thold>0
* scripts/getSummary.sh                : summarize SV's

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
* ls2in.pl                             : added "-word" parameter : word only sample names

* examples/HPRC/Illumina/in.url        : URLs of 40 HPRC Illumina samples
* examples/HPRC/HiFi/in.url            : URLs of	40 HPRC	HiFi     samples

# 2024/12/03 #

* scripts/uniqVcf.pl		       : convert REF/ALT to upper case
* scripts/maxVcf.pl		       : convert REF/ALT to upper case

# 2025/01/07 #

* RefSeq/{CDS,RNR,TRN}.bed.gz*         : use Hugo gene names (MT-*)  
* scripts/install_sysprerequisites.sh  : repace "default-jdk default-jre" with "openjdk-18-jdk"
