#!/bin/bash -eux

###############################################################################

M1=mutect2 ; M2=varscan ; M3=freebayes

. $HP_SDIR/init.$M1.sh ; run.sh | tee run.$M1.sh | bash
. $HP_SDIR/init.$M2.sh ; run.sh | tee run.$M2.sh | bash
. $HP_SDIR/init.$M3.sh ; run.sh | tee run.$M3.sh | bash
. $HP_SDIR/init.sh

mkdir -p $HP_ODIR
merge3Vcfs.sh
