# 2024/05/23 #

* updated RefSeq/{CDS,COMPLEX,TRN}.bed.gz.* ; annotate the overlapping features using both names
* remove  RefSeq/STOP.vcf.gz*               ; STOP gain annotation present in MLC.vcf.gz*

# 2024/07/30 #

* scripts/annotateVcf.sh   : added dbSNP INFO
* scripts/filter.sh        : added dbSNP INFO
* scripts/filterVcf.pl     : added the FORMAT/AD tag
* scripts/fixmutect2Vcf.pl : added the FORMAT/AD tag
* scripts/*.vcf            : added the FORMAT/AD tag

* scripts/filter2.sh       : added new script to filter the original mutect2.mutect2 calls => mutect2.mutect2.00.orig
