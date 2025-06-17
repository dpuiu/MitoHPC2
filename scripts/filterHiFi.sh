#!/usr/bin/env bash
set -eux

#########################################################################################################################################

# Program that runs the heteroplasmy pipeline on a single HiFi sample

# Input arguments
#  1: sample names
#  2: BAM/CRAM alignment file; full path
#  3: output prefix; full path

#########################################################################################################################################

#set variables

export S=$1             # sample name
IDIR=`dirname "$2"`     # dir path
I=${2%.*}               # bam/cram file name prefix
X=${2##*.}              #                    extension
O=$3                    # output prefix
M="${4:-$HP_M}"         # 2024/03/12 

OS=$O.$M                # output prefix + snv_caller
OSS=$OS.$M
#RG="@RG\tID:$S\tSM:$S\tPL:HiFi"
 RG="@RG\tID:$S\tSM:$S\tPL:HiFi"

export PC="0.95"

ODIR=`dirname "$3"`; mkdir -p $ODIR

MAXDP=2000
MINAF=0.01
#########################################################################################################################################
# test if count and VCF output files exist; exit if they do

if [ $HP_I -lt 1 ] && [ -s $O.count ]    ; then exit 0 ; fi
if [ $HP_I -eq 1 ] && [ -s $OS.00.vcf  ] ; then exit 0 ; fi
if [ $HP_I -ge 2 ] && [ -s $OSS.00.vcf ] ; then exit 0 ; fi

#########################################################################################################################################
# test alignment file exists, is indexed and sorted by coordinates
test -s $2
test -s $2.bai || test -s $2.crai || test -s $I.bai || test -s $I.crai

# generate index and indexstats files
if [ ! -s $O.idxstats ] ; then
  if [ ! -s $I.idxstats ] ; then samtools idxstats $2 > $O.idxstats 
                            else cp $I.idxstats $O.idxstats 
  fi
fi

# get read counts
if [ ! -s $O.count ]; then cat $O.idxstats | idxstats2count.pl -sample $S -chrM $HP_RMT > $O.count ; fi

# test if there are any MT reads
MTCOUNT=`tail -1 $O.count| cut -f4`
if [ $MTCOUNT -lt 1 ]; then echo "ERROR: There are no MT reads in $2; plese remove it from $HP_IN" ; exit 1 ; fi


if [ ! -s $O.fa ] ; then
  if [ "$X" == "cram" ] ; then T="-T $HP_RFILE.fa" ; else T="" ; fi
  #MTCOUNT=`samtools view -F 0x900 $2 $HP_RMT $T -c`
  #COUNT=`samtools view -F 0x900 $2 $T -c`
  #echo -e "$S\t$COUNT\t$COUNT\t$MTCOUNT" > $O.count
  #if [ $MTCOUNT -lt 1 ]; then echo "ERROR: There are no MT reads in $2; plese remove it from $HP_IN" ; exit 1 ; fi
  if [ $HP_L ]; then R=`tail -1 $O.count | perl -ane '$R=$ENV{HP_L}/($F[-1]+1);  if($R<1) { print "-s $R"} else {print ""} '` ; else R="" ; fi
  samtools view -F 0x900 $R $2 $HP_RMT $T    | perl -ane 'print "$F[0]\t",length($F[9]),"\n";' | sort > $O.len
  samtools view -F 0x900 $R $2 $HP_RMT $T -b | samtools fasta > $O.fa
fi

if [ ! -s $O.bam ] ; then
  cat $O.fa | minimap2 -x map-ont $HP_RDIR/$HP_MT.fa /dev/stdin -R $RG -a --eqx | samtools view -b | samtools sort | samtools view -h > $O.sam
  samtools view -b $O.sam | bedtools bamtobed | cut -f 1-4 | bed2bed.pl  | count.pl -i 3 -j 4 | sort | join $O.len -  -a 1 --nocheck-order | \
    perl -ane 'print "$F[0]\t$F[1]\t$F[2]\t",$F[2]/$F[1],"\n"' | sort -k4,4nr | tee $O.score  | perl -ane 'print  if($F[-1]>$ENV{PC});' | cut -f1 | \
    samtools view -N /dev/stdin $O.sam  -b > $O.bam
  samtools index $O.bam
  rm $O.sam
fi

# count aligned reads; compute cvg; get coverage stats; get split alignments

if [ ! -s $O.cvg ]    ; then cat $O.bam | bedtools bamtobed | bedtools genomecov -i - -g $HP_RDIR/$HP_MT.fa.fai  -d | tee $O.cvg  | cut -f3 | st.pl -sample $S > $O.cvg.stat ; fi
if [ ! -f $O.sa.bed ] ; then samtools view -h $O.bam | sam2bed.pl | sort -k1,1 -k2,2n -k3,3n | tee $O.bed | grep -P '\d\d+D\t|\d\d+I\t' > $O.sa.bed ; fi

if [ ! -s $OS.vcf ] ; then
  if [ "$M" == "bcftools" ] ; then
    bcftools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -d $MAXDP | bcftools call --ploidy 2 -mv -Ov > $OS.vcf
  elif [ "$M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction $MINAF $O.bam -f $HP_RDIR/$HP_MT.fa  > $OS.vcf
  elif [ "$M" == "varscan" ] ; then
    samtools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -r $HP_MT -B -d $MAXDP  | java -jar $HP_JDIR/VarScan.jar mpileup2snp   --min-coverage $HP_DP -B --variants --min-var-freq $MINAF --output-vcf 1 > $OS.orig.vcf
    samtools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -r $HP_MT -B -d $MAXDP  | java -jar $HP_JDIR/VarScan.jar mpileup2indel --min-coverage $HP_DP -B --variants --min-var-freq $MINAF --output-vcf 1 | grep -v "^#" >> $OS.orig.vcf
    cat $HP_SDIR/$M.vcf >  $OS.vcf
    cat $HP_RDIR/$HP_MT.fa.fai | perl -ane 'print "##contig=<ID=$F[0],length=$F[1]>\n"' >> $OS.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$S" >> $OS.vcf
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT:DP:AD:AF\t[%GT:%DP:%AD:%FREQ]\n'  $OS.orig.vcf |  perl -lane 'print "$1:",int($2*100+.5)/10000 if(/(.+):(.+)%$/);' | sort -k2,2n >> $OS.vcf
  else
    echo "Unsuported SNV caller 1"
    exit 1
  fi
fi

if [ ! -s $OS.00.vcf ] ; then
  #cat $HP_SDIR/$M.vcf  > $OS.vcf
  #fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OS.vcf
  bcftools norm -m-any  -f $HP_RDIR/$HP_MT.fa $OS.vcf|  fix${M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header  > $OS.fix.vcf
  cat $OS.fix.vcf | maxVcf.pl | bedtools sort -header |tee $OS.max.vcf | bgzip -f -c > $OS.max.vcf.gz ; tabix -f $OS.max.vcf.gz
  cat $OS.fix.vcf | filterVcf.pl -sample $S -source $M -header $HP_SDIR/$M.vcf -depth $HP_DP | uniqVcf.pl | bedtools sort -header > $OS.00.vcf
  annotateVcf.sh $OS.00.vcf
fi

if [ ! -s $OS.haplogroup ] ; then
  if [ "$HP_O" == "Human" ] ; then
    if [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "chrM" ] ; then  java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup
    elif [ "$HP_MT" == "RSRS" ] ; then                           java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup --rsrs
    fi
  fi
fi

if  [ ! -s $OS.fa ]  ; then
  bcftools consensus -f $HP_RDIR/$HP_MT.fa $OS.max.vcf.gz | perl -ane 'chomp; if($.==1) { print ">$ENV{S}\n" } else { s/N//g; print } END {print "\n"}' > $OS.fa
  samtools faidx $OS.fa
fi

#################################################

if [ ! -s $OS.bam ] ; then
   cat $O.fa |  minimap2 -x map-ont $OS.fa /dev/stdin -R $RG -a  --eqx | samtools view -b | samtools sort | samtools view -h > $OS.sam
   bedtools bamtobed -i $OS.sam | cut -f 1-4 | bed2bed.pl  | count.pl -i 3 -j 4  | sort | join $O.len -  -a 1 --nocheck-order | \
     perl -ane 'print "$F[0]\t$F[1]\t$F[2]\t",$F[2]/$F[1],"\n"' | perl -ane 'print  if($F[-1]>$ENV{PC});'  | cut -f1 | \
     samtools view -N /dev/stdin $OS.sam  -b > $OS.bam
  samtools index $OS.bam
  rm $OS.sam $O.fa $O.len
fi

if [ ! -s $OS.cvg ]   ; then bedtools genomecov -d -ibam $OS.bam | tee $OS.cvg  | cut -f3 | st.pl -sample $S > $OS.cvg.stat  ; fi
if [ ! -f $OS.sa.bed ] ; then samtools view -h $OS.bam | sam2bed.pl | sort -k1,1 -k2,2n -k3,3n  | tee $OS.bed | grep -P '\d\d+D\t|\d\d+I\t'  > $OS.sa.bed ; fi
if [ ! -f $OS.merge.bed ] ; then cat $OS.bed | grep -P '\d=\t|\d\M\t' | bedtools merge | bed2bed.pl > $OS.merge.bed ; fi


if [ ! -s $OSS.00.vcf ] ; then
  cat $HP_SDIR/$M.vcf > $OSS.00.vcf
  echo "##sample=$S" >> $OSS.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OSS.00.vcf

  if [ "$M" == "bcftools" ] ; then
    bcftools mpileup -f $OS.fa $OS.bam -d $MAXDP | bcftools call --ploidy 2 -mv -Ov > $OSS.vcf
  elif [ "$M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction $MINAF $OS.bam -f $OS.fa  > $OSS.vcf
 elif [ "$M" == "varscan" ] ; then
    samtools mpileup -f $OS.fa $OS.bam -r $S -B -d $MAXDP  | java -jar $HP_JDIR/VarScan.jar mpileup2snp   --min-coverage $HP_DP -B --variants --min-var-freq $MINAF --output-vcf 1 > $OSS.orig.vcf
    samtools mpileup -f $OS.fa $OS.bam -r $S -B -d $MAXDP  | java -jar $HP_JDIR/VarScan.jar mpileup2indel --min-coverage $HP_DP -B --variants --min-var-freq $MINAF --output-vcf 1 | grep -v "^#" >> $OSS.orig.vcf
    cat $HP_SDIR/$M.vcf >  $OSS.vcf
    echo "##sample=$S" >> $OSS.vcf
    cat $OS.fa.fai | perl -ane 'print "##contig=<ID=$F[0],length=$F[1]>\n"' >> $OSS.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$S" >> $OSS.vcf
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT:DP:AD:AF\t[%GT:%DP:%AD:%FREQ]\n'  $OSS.orig.vcf |  perl -lane 'print "$1:",int($2*100+.5)/10000 if(/(.+):(.+)%$/);' | sort -k2,2n >> $OSS.vcf
  else
    echo "Unsuported SNV caller 1"
    exit 1
  fi

  cat  $OSS.vcf |\
    bcftools norm -m-any  -f $OS.fa  | \
    fix${M}Vcf.pl  | \
    bedtools sort -header | \
    fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -rlen $HP_MTLEN -mfile $OS.max.vcf  | \
    cat $OS.max.vcf - | \
    filterVcf.pl -sample $S -source $M |  bedtools sort  >> $OSS.00.vcf

  annotateVcf.sh $OSS.00.vcf
fi
