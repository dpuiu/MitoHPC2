#!/usr/bin/env bash
set -x

##############################################################################################################

# Program that downcramples a SAM file; ~ max MCOUNT alignments starting at each position are allowed

# Examples:
#      downcrampleSam.sh in out
#      downcrampleSam.sh in out 10
#      downcrampleSam.sh in out 10 chrM
#      downcrampleSam.sh in out 10 "chrM chr1:629084-634672 chr1:76971223-76971280 chr17:22521208-22521639"

##############################################################################################################

IN=$1                                                                                   # input bam prefix(including path)
OUT=$2                                                                                  # output bam prefix(including path)
MCOUNT=${3:-10}                                                                         # max number of reads starting at a certain pos
R=${4:-"chrM chr1:629084-634672 chr17:22521208-22521639"}                               # MT+NUMT regions (MT 1st) 
                                                                                        # chr1:76971223-76971280 & chr17:22521208-22521639 align to chrM 5'/3'
                                                                                        # chr1:629084-634672 : 5Kb NUMT

test -s $IN
IP=`basename "$IN" .sam`
IP=`basename "$IP" .bam`
IP=`basename "$IP" .cram`

OP=`basename "$OUT" .bam`

test -s $HP_RDIR/$HP_RNAME.fa
test -s $HP_RDIR/$HP_RMT.fa.fai 

################################

#downcrample and index
if [ ! -s $OUT ] ; then
  samtools view -h $IN $R -F 0x90C -T $HP_RDIR/$HP_RNAME.fa | \
    downsampleSam.pl -max $MCOUNT -hg38 | grep -v "^\@" | cut -f1 | sort | uniq -d  > $OUT.ids
  cat $OUT.ids | samtools view -b -N /dev/stdin  $IN $R  -T $HP_RDIR/$HP_RNAME.fa > $OUT

  samtools index    $OUT
  samtools idxstats $OUT > $OP.idxstats

  cat $IN  | bedtools bamtobed -cigar | grep ^$HP_RMT | bedtools genomecov -d -i - -g $HP_RDIR/$HP_RMT.fa.fai > $IN.cvg #| getSummary.pl -i -1 -t $IP > $IP.cvg.summary
  cat $OUT | bedtools bamtobed -cigar | grep ^$HP_RMT | bedtools genomecov -d -i - -g $HP_RDIR/$HP_RMT.fa.fai > $OP.cvg #| getSummary.pl -i -1 -t $OP > $OP.cvg.summary
fi

