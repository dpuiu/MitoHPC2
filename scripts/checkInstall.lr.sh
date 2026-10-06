#!/usr/bin/env bash
#set -ex
set -euo pipefail

####################################################################################

echo "########################"  >> checkInstall.log
echo "EXECUTABLES PATHS LONG READS:"  >> checkInstall.log

#which go
which singularity >> checkInstall.log

####################################################################################
#clair3/deepvariant support
if [ ! -f $HP_BDIR/clairs-to.sif ] && [ ! -d "$HOME/clairs-to_sandbox" ] && [ ! -L "$HOME/clairs-to_sandbox" ]; then
    echo "ERROR: Directory $HP_BDIR/clairs-to_sandbox does not exist" >&2
    exit 1
fi

if [ ! -f $HP_BDIR/deepsomatic.sif ] && [ ! -d "$HOME/deepsomatic_sandbox" ] && [ ! -L "$HOME/deepsomatic_sandbox" ]; then
    echo "ERROR: Directory $HP_BDIR/deepsomatic_sandbox does not exist" >&2
    exit 1
fi

cat checkInstall.log
echo Success!
