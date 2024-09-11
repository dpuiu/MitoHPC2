# 2024/05/23 #

* updated RefSeq/{CDS,COMPLEX,TRN}.bed.gz.* ; annotate the overlapping features using both names
* remove  RefSeq/STOP.vcf.gz*               ; STOP gain annotation present in MLC.vcf.gz*

# 2024/07/30 #

* scripts/annotateVcf.sh     : added dbSNP INFO
* scripts/filter.sh          : added dbSNP INFO
* scripts/filterVcf.pl       : added the FORMAT/AD tag
* scripts/fixmutect2Vcf.pl   : added the FORMAT/AD tag
* scripts/*.vcf              : added the FORMAT/AD tag

* scripts/filter2.sh         : added new script to filter the original mutect2.mutect2 calls => mutect2.mutect2.00.orig

# 2024/07/31 #

* scripts/filter2.sh	     : file path update

# 2024/08/01 #

* scripts/fixsnpPos.pl	     : using lower cases for 2nd iteration "flipped" SNVs 
* scripts/differenceVcf.pl   : allow for lower eq upper cases 

# 2024/08/02 #

* scripts/filter1.sh	     : added new script to re-filter the original mutect2 calls and add AD tag to mutect2.00 file

# 2024/09/22 # 

* scripts/README.md          : updated README.med file
* scripts/run.hifi.sh        : to process HiFi reads
* scripts/fixbcftoolsVcf.pl  : added the FORMAT/AD tag => can process HiFi reads
* scripts/fixfreebayesVcf.pl : added the FORMAT/AD tag
* scripts/fixmutserveVcf.pl  : added the FORMAT/AD tag
