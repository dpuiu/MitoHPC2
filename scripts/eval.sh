#!/bin/bash -e

R=$1	# reference VCF
Q=$2	# euqery VCF

test -s $R
test -s $Q

###################################
echo `basename $PWD`
echo
echo "$MT:all"
cmpVcf.pl $R $Q -sm
cmpVcf.pl $R $Q -snp -sm
cmpVcf.pl $R $Q -ins -sm
cmpVcf.pl $R $Q -del -sm

echo
echo "$MT:homoplasmies"
cmpVcf.pl $R $Q -hom -sm
cmpVcf.pl $R $Q -snp -hom -sm
cmpVcf.pl $R $Q -ins -hom -sm
cmpVcf.pl $R $Q -del -hom -sm

echo
echo "$MT:heteroplasmies"
cmpVcf.pl $R $Q -snp -sm
cmpVcf.pl $R $Q -snp -het -sm
cmpVcf.pl $R $Q -ins -het -sm
cmpVcf.pl $R $Q -del -het -sm


