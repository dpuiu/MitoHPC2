#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that annotates a VCF file

# Input:  VCF file
# Output: VCF file

###############################################################################################################

export O=$1   # input/output vcf file
test  -s $O

bgzip -f $O ;tabix -f $O.gz

if [ -s $HP_RDIR/dbSNP.vcf.gz ] ; then bcftools annotate -a $HP_RDIR/dbSNP.vcf.gz $O.gz -c "ID,INFO"    > $O ; bgzip -f $O ; tabix -f $O.gz ; fi

ls $HP_RDIR/*.bed.gz |                        perl -ane 'chomp; /.+\/(.+).bed.gz/; print "bcftools annotate -a $_ $ENV{O}.gz -c \"CHROM,FROM,TO,$1\" -h  <(echo \"##INFO=<ID=$1,Number=1,Type=String,Description=\"$1\">\")  > $ENV{O} ; bgzip -f $ENV{O} ; tabix -f $ENV{O}.gz\n";' | bash
ls $HP_RDIR/*.vcf.gz | grep -v dbSNP.vcf.gz | perl -ane 'chomp; print "bcftools annotate -a $_ $ENV{O}.gz -c \"INFO\"  > $ENV{O} ; bgzip -f $ENV{O} ; tabix -f $ENV{O}.gz\n";' | bash

if [ ! -s $O ] ; then zcat $O.gz | grep -v "^##bcftools" > $O ; fi
rm -f $O.gz $O.gz.tbi
