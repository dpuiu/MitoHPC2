#!/usr/bin/env bash
#set -ex
set -euo pipefail

test -s $HP_RDIR/$HP_RNAME.fa
test -s $HP_RDIR/$HP_MT.fa
test -s $HP_RDIR/$HP_MT2.fa
test -s $HP_RDIR/$HP_MTC.fa
test -s $HP_RDIR/$HP_MTR.fa
test -s $HP_RDIR/$HP_NUMT.fa

test -s $HP_RDIR/$HP_MTC.bwt
test -s $HP_RDIR/$HP_NUMT.bwt

test -s $HP_RDIR/$HP_MT.dict

echo Success!
