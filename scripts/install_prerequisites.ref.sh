#!/usr/bin/env bash
set -ux # do not add -e

#####################################################################################

#if [ ! -s $HP_RDIR/$HP_RNAME.fa ] ; then # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_RNAME.fa.fai || $# == 1 && $1 == "-f" ]] ; then
  wget -qO- $HP_RURL | zcat -f > $HP_RDIR/$HP_RNAME.fa
  #wget -q $HP_RURL -O $HP_RDIR/$HP_RNAME.fa  # wget on fedora does not download ftp links
  #curl -L $HP_RURL -o $HP_RDIR/$HP_RNAME.fa
  samtools faidx $HP_RDIR/$HP_RNAME.fa
fi

#if [ ! -s $HP_RDIR/$HP_MT.fa ] ; then  # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_MT.dict || $# == 1 && $1 == "-f"  ]] ; then
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RMT > $HP_RDIR/$HP_MT.fa
  samtools faidx $HP_RDIR/$HP_MT.fa
  rm $HP_RDIR/$HP_MT.dict
  java $HP_JOPT -jar $HP_BDIR/gatk.jar CreateSequenceDictionary --REFERENCE $HP_RDIR/$HP_MT.fa --OUTPUT $HP_RDIR/$HP_MT.dict
fi

#if [ ! -s $HP_RDIR/$HP_NUMT.fa ] ; then # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_NUMT.bwt || $# == 1 && $1 == "-f"  ]] ; then
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RNUMT > $HP_RDIR/$HP_NUMT.fa
  bwa index $HP_RDIR/$HP_NUMT.fa -p $HP_RDIR/$HP_NUMT
fi

#if [ ! -s $HP_RDIR/$HP_MTC.fa ] ; then  # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_MTC.dict || $# == 1 && $1 == "-f" ]] ; then
  circFasta.sh $HP_MT $HP_RDIR/$HP_MT $HP_E $HP_RDIR/$HP_MTC
fi

#if [ ! -s $HP_RDIR/$HP_MTR.fa ] ; then # 2023/04/26
if [[ ! -s $HP_RDIR/$HP_MTR.dict || $# == 1 && $1 == "-f" ]] ; then
  rotateFasta.sh $HP_MT $HP_RDIR/$HP_MT $HP_E $HP_RDIR/$HP_MTR
fi
