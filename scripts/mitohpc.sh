#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program launches the MitoHPC pipeline with default parameters

###############################################################################################################

if [ -f ./init.sh ] ; then  . ./init.sh ; else  . $HP_SDIR/init.sh ; fi
find $HP_ADIR/  -name "*.bam" -o -name "*.cram" -readable | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN
$HP_SDIR/run.sh | tee run.all.sh | bash
