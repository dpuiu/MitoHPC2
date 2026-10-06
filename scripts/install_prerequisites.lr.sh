#!/usr/bin/env bash
set -ux # do not add -e

#which go
which singularity

####################################################################################

#pull docker images as SIF's
singularity pull $HP_BDIR/clairs.sif	  docker://hkubal/clairs-to:v0.4.2
singularity pull $HP_BDIR/deepsomatic.sif docker://google/deepsomatic:1.10.0 # 1.9.0

#or download pre-generated SIF's
#wget ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/MitoHPC2/bin/clairs-to_v0.4.2.sif   -O $HP_BDIR/clairs-to.sif
#wget ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/MitoHPC2/bin/deepsomatic_1.10.0.sif -O $HP_BDIR/deepsomatic.sif

#rockfish
#singularity build --sandbox ~/clairs-to_sandbox/          $HP_BDIR/clairs-to.sif
#singularity build --sandbox ~/deepsomatic_sandbox/        $HP_BDIR/deepsomatic.sif
#singularity build --sandbox ~/deepsomatic_sandbox_1.10.0/ $HP_BDIR/deepsomatic_1.10.0.sif ; ln -s ~/deepsomatic_sandbox_1.10.0/ ~/deepsomatic_sandbox
