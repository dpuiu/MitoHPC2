#!/usr/bin/env bash
#set -eux

#########################################################################################################################################

# Program that runs the heteroplasmy pipeline on a single sample

# Input arguments
#  1: sample names
#  2: BAM/CRAM alignment file; full path
#  3: output prefix; full path
#  4:  SNV caller

#########################################################################################################################################

#set variables
export S=$1		# sample name
N=`basename "$2" .bam`
export N=`basename "$N" .cram`
IDIR=`dirname "$2"`
I=$IDIR/$N
ODIR=`dirname "$3"`; mkdir -p $ODIR 

M="${4:-$HP_M}" # 2024/03/11 

O=$3
OR=${O}R

ON=$O.$HP_NUMT
OS=$O.$M

OSC=${OS}C
OSR=${OS}R

OSS=$OS.$M
OSSR=${OSS}R

MINAF=0.01
MAXDP=2000

#########################################################################################################################################
# test if count and VCF output files exist; exit if they do

if [ $HP_I -lt 1 ] && [ -s $O.bam ]    ; then exit 0 ; fi
if [ $HP_I -eq 1 ] && [ -s $OS.00.vcf  ] ; then exit 0 ; fi
if [ $HP_I -ge 2 ] && [ -s $OSS.00.vcf ] ; then exit 0 ; fi

#########################################################################################################################################
# test alignment file exists and is sorted by coordinates

# test alignment input file exists , is indexed and sorted 
test -s $2
test -s $2.bai || test -s $2.crai
samtools view -H $2 | grep -m 1 -P "^@HD.+coordinate$" > /dev/null

# test references match
#RCOUNT=`samtools view -H $2 | grep -c "^@SQ"`
#if [ $HP_RCOUNT != $RCOUNT ] ; then echo "ERROR: HP_RCOUNT=$HP_RCOUNT does not match the number of \@SQ lines=$RCOUNT in $2"; exit 1 ; fi

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

#########################################################################################################################################
# subsample reads

if [ ! -s $O.fq ] ; then
  R=""
  if [ $HP_L ]; then
    R=`tail -1 $O.count  | perl -ane '$R=$ENV{HP_L}/($F[-1]+1); print $R if($R<1)'`
    if [ $R ] ; then R="-s $R"  ; fi
  fi

  samtools view $R $2 $HP_RMT $HP_RNUMT -bu -F 0x900 -T $HP_RDIR/$HP_RNAME.fa -@ $HP_P   | \
    samtools sort -n -O SAM -m $HP_MM -@ $HP_P |  \
    perl -ane 'if(/^@/) {print} elsif($P[0] eq $F[0]) {print $p,$_}; @P=@F; $p=$_;' | \
    samblaster $HP_DOPT --addMateTags   | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout $HP_FOPT  > $O.fq
fi
#########################################################################################################################################
# realign subsampled reads

if  [ ! -s $O.bam ] ; then
  cat $O.fq | \
    bwa mem $HP_RDIR/$HP_MTC - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view  -F 0x90C -h | \
    tee $O.sam  | \
    samtools view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix | \
    count.pl -i 3 -j 4  | sort > $O.score

  cat $O.fq | \
    bwa mem $HP_RDIR/$HP_NUMT - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix  | \
    count.pl -i 3 -j 4  | sort > $ON.score

  join $O.score $ON.score -a 1 --nocheck-order | perl -ane 'next if(@F==3 and $F[2]>$F[1]);print join "\t",@F;print "\n"' | \
     intersectSam.pl $O.sam - > $O.sam2 ;  mv $O.sam2 $O.sam

  cat $O.sam | \
     circSam.pl -ref_len $HP_RDIR/$HP_MT.fa.fai -offset 0 | \
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P -o $O.bam 
  samtools index $O.bam

  cat $O.sam | \
     circSam.pl -ref_len $HP_RDIR/$HP_MT.fa.fai -offset $HP_E | \
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P -o $OR.bam 
  samtools index $OR.bam
  rm -f $O.*sam $O.score 
fi

#########################################################################################################################################
# count aligned reads; compute cvg; get coverage stats; get split alignments

#if [ ! -s $O.count ]  ; then samtools idxstats $O.bam  -@ $HP_P | idxstats2count.pl -sample $S -chrM $HP_MT > $O.count ; fi
if [ ! -s $O.cvg ]    ; then cat $O.bam  | bedtools bamtobed -cigar | grep "^$HP_MT" | bedtools genomecov -i - -g $HP_RDIR/$HP_MT.fa.fai  -d | tee $O.cvg  | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $O.cvg.stat  ; fi
if [ ! -f $O.sa.bed ] ; then samtools view -h $O.bam  -@ $HP_P | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $O.sa.bed ; fi

if [ $HP_I -lt 1 ]  ; then exit 0 ; fi

#########################################################################################################################################
# compute SNVs using mutect2/mutserve/freebayes

if [ ! -s $OS.vcf ] ; then
  if [ "$M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $HP_RDIR/$HP_MT.fa -I $O.bam       -O $OS.orig.vcf $HP_GOPT --native-pair-hmm-threads $HP_P --callable-depth 6 --max-reads-per-alignment-start 0 -min-AF $MINAF # --mitochondria-mode 
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $HP_RDIR/$HP_MT.fa -V $OS.orig.vcf -O $OS.vcf --min-reads-per-strand 2

    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $HP_RDIR/$HP_MTR.fa -I $OR.bam        -O $OSR.orig.vcf $HP_GOPT --native-pair-hmm-threads $HP_P  --callable-depth 6 --max-reads-per-alignment-start 0 -L "$HP_MT:285-315" -min-AF $MINAF # "$HP_MT:16254-16284"  
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $HP_RDIR/$HP_MTR.fa -V $OSR.orig.vcf  -O $OSR.vcf --min-reads-per-strand 2
    cat $OSR.vcf | perl -ane 'next if(/^#/) ; $F[1]=($F[1]-$ENV{HP_E})%$ENV{HP_MTLEN} ; print join "\t",@F; print "\n" ' >> $OS.vcf
    rm -f $OS.orig.* $OS.vcf.*

  elif [ "$M" == "mutserve" ] ; then
    if [ "$HP_MT" == "chrM" ] ||  [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "RSRS" ] ; then
      java $HP_JOPT -jar $HP_JDIR/mutserve.jar call --deletions --insertions --level $MINAF --output $OS.vcf --reference $HP_RDIR/$HP_MT.fa $O.bam
      #mv $O.txt $OS.txt
    else
      echo "Wrong mutserve reference"
      exit 1
    fi
  elif [ "$M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction $MINAF $O.bam -f $HP_RDIR/$HP_MT.fa  > $OS.vcf
  elif [ "$M" == "varscan" ] ; then
    #samtools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -r $HP_MT -B -d $MAXDP  | tee \
    #  >(java -jar $HP_JDIR/VarScan.jar pileup2snp   -B  --variants --min-var-freq $MINAF > $OS.snp.txt) \
    #  >(java -jar $HP_JDIR/VarScan.jar pileup2indel -B  --variants --min-var-freq $MINAF > $OS.indel.txt) > /dev/null
    samtools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -r $HP_MT -B -d $MAXDP  | java -jar $HP_JDIR/VarScan.jar pileup2snp   -B  --variants --min-var-freq $MINAF > $OS.txt
    samtools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -r $HP_MT -B -d $MAXDP  | java -jar $HP_JDIR/VarScan.jar pileup2indel -B  --variants --min-var-freq $MINAF | grep -v "^#" >> $OS.txt

    cat $HP_SDIR/$M.vcf  > $OS.vcf
    fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OS.vcf
    cat $OS.txt| \
      perl -ane 'next if($.==1 or $F[2]=~/N/i or $F[3]=~/N/i); $DP=$F[4]+$F[5];$F[6]=~/(.+)%/;$AF=$1/100; print join "\t",($F[0],$F[1],".",$F[2],$F[-1],".",".",".","GT:DP:AF","0/1:$DP:$AF");print "\n";' | \
      perl -ane 'if($F[4]=~/\+(.+)/) {($F[4],$F[7])=("$F[3]$1","INDEL")} elsif($F[4]=~/\-(.+)/) {($F[3],$F[4],$F[7])=("$F[3]$1",$F[3],"INDEL")} print join "\t",@F; print "\n";' | uniqVcf.pl | sort -k1,1 -k2,2n >> $OS.vcf
  else
    echo "Unsuported SNV caller 1"
    exit 1
  fi
  rm -f $OSR.*
fi

if [ ! -s $OS.00.vcf ] ; then
  # filter SNVs ; to update fixmuserveVcf.pl !!! 2024/03/04
  bcftools norm -m-any -f $HP_RDIR/$HP_MT.fa  $OS.vcf   | fix${M}Vcf.pl -file $HP_RDIR/$HP_MT.fa  | bedtools sort -header> $OS.fix.vcf
  cat $OS.fix.vcf | maxVcf.pl | bedtools sort -header |tee $OS.max.vcf | bgzip -f -c > $OS.max.vcf.gz ; tabix -f $OS.max.vcf.gz
  cat $OS.fix.vcf | filterVcf.pl -sample $S -source $M -header $HP_SDIR/$M.vcf -depth $HP_DP | uniqVcf.pl | bedtools sort -header > $OS.00.vcf
  annotateVcf.sh $OS.00.vcf
fi

#########################################################################################################################################
#  identify SVs

if [ $HP_V ] && [ "$HP_V" == "gridss" ] ; then
  OV=$O.$HP_V
  if [ ! -s $OV.00.vcf ] ; then
    gridss --jar $HP_JDIR/gridss.jar -r $HP_RDIR/$HP_MT.fa -o $OV.vcf.gz $O.bam  -t 1 -w $OV
    cat $HP_SDIR/gridss.vcf > $OV.fix.vcf
    bcftools view -i 'FILTER="PASS"' $OV.vcf.gz | bcftools query  -f "%CHROM\t%POS\t.\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT:DP:AF\t[%GT:%REF:%AF]\n" | grep -v -P 'h\t'  >> $OV.fix.vcf
    cat $OV.fix.vcf | filterVcf.pl -sample $S -source $HP_V -depth $HP_DP  > $OV.00.vcf 
    rm -rf $OV $OV.vcf.gz.* $OV.fix.vcf
  fi
fi

########################################################################################################################################
# get haplogroup

if [ "$HP_O" == "Human" ] ; then
  if [ ! -s $OS.haplogroup ] ; then
    if [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "chrM" ] ; then java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup
    elif [ "$HP_MT" == "RSRS" ] ; then                          java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup --rsrs
    fi

    if [ -s $OS.haplogroup ] ; then java $HP_JOPT -jar $HP_JDIR/haplocheck.jar --out $OS.haplocheck $OS.vcf ; fi
  fi
fi

#########################################################################################################################################
# get new consensus; format reference

if  [ ! -s $OS.fa ]  ; then
  bcftools consensus -f $HP_RDIR/$HP_MT.fa $OS.max.vcf.gz -H A | perl -ane 'chomp; if($.==1) { print ">$ENV{S}\n" } else { s/N//g; print } END {print "\n"}' > $OS.fa
  rm -f $OS.max.vcf.gz $OS.max.vcf.gz.tbi
  samtools faidx $OS.fa
  rm -f $OS.dict
  java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $OS.fa --OUTPUT $OS.dict

  export MTLEN=`cut -f2 $OS.fa.fai`

  circFasta.sh   $N $OS $HP_E $OSC
  rotateFasta.sh $N $OS $HP_E $OSR
fi

########################################################################################################################################
# realign reads; check coverage

if  [ ! -s $OS.bam ] ; then
  cat $O.fq | \
    bwa mem $OSC - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view  -F 0x90C -h | \
    tee $OS.sam  | \
    samtools  view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix  | \
    count.pl -i 3 -j 4  | sort > $OS.score

  join $OS.score $ON.score -a 1 --nocheck-order | perl -ane 'next if(@F==3 and $F[2]>$F[1]);print join "\t",@F;print "\n"' |\
     intersectSam.pl $OS.sam - > $OS.sam2 ; mv $OS.sam2 $OS.sam

  cat $OS.sam |\
     circSam.pl -ref_len $OS.fa.fai -offset 0 |\
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P -o $OS.bam 
  samtools index $OS.bam

  cat $OS.sam | \
     circSam.pl -ref_len $OS.fa.fai -offset $HP_E | \
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P -o  $OSR.bam 
  samtools index $OSR.bam

  bedtools bamtobed -i $OS.bam -ed | perl -ane 'print if($F[-2]==0);' | bedtools merge -d -3 | bed2bed.pl -min 3 > $OS.merge.bed
  rm -f $OS.*sam $OS.score
  #!!! rm -f $O.fq $ON.score
fi

### !!! rm -f $O.bam*  $OSR.bam*
rm -f  $OSC.*

# exit if the number of iterations is set to 1
if [ $HP_I -lt 2 ] || [ $M == "mutserve" ] ; then
  rm -f $OS.bam*  $OSR.bam*
  #!!! rm -f $O.fq $ON.score $O.bam* $OR.bam*
  exit 0
fi

#########################################################################################################################################
# 2nd iteration
# count aligned reads; compute cvg; compute cvg stats; gets split alignments

#if [ ! -s $OS.count ]  ; then samtools idxstats $OS.bam  -@ $HP_P | idxstats2count.pl -sample $S -chrM $S > $OS.count ; fi
if [ ! -s $OS.cvg ]    ; then cat $OS.bam  | bedtools bamtobed -cigar | grep "^$S" | bedtools genomecov -i - -g $OS.fa.fai -d  | tee $OS.cvg  | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $OS.cvg.stat  ; fi
if [ ! -f $OS.sa.bed ] ; then samtools view -h $OS.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $OS.sa.bed ; fi

#########################################################################################################################################
# compute SNP/INDELs using mutect2/freebayes
if [ ! -s $OSS.vcf ] ; then
  if [ "$M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $OS.fa -I $OS.bam       -O $OSS.orig.vcf  $HP_GOPT --native-pair-hmm-threads $HP_P --callable-depth 6 --max-reads-per-alignment-start 0  -min-AF $MINAF # --mitochondria-mode
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $OS.fa -V $OSS.orig.vcf -O $OSS.vcf  --min-reads-per-strand 2

    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $OSR.fa -I $OSR.bam       -O $OSSR.orig.vcf  $HP_GOPT --native-pair-hmm-threads $HP_P --callable-depth 6 --max-reads-per-alignment-start 0 -L "$S:285-315"   -min-AF $MINAF
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $OSR.fa -V $OSSR.orig.vcf -O $OSSR.vcf  --min-reads-per-strand 2
    cat  $OSSR.vcf | perl -ane 'next if(/^#/) ; $F[1]=($F[1]-$ENV{HP_E})%$ENV{MTLEN} ; print join "\t",@F; print "\n" '  >> $OSS.vcf
    rm   $OSS.orig.* $OSSR.* $OSS.vcf.*

  elif [ "$M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction $MINAF $OS.bam -f $OS.fa  > $OSS.vcf

  elif [ "$M" == "varscan" ] ; then
    samtools mpileup -f $OS.fa $OS.bam -B -d $MAXDP | java -jar $HP_JDIR/VarScan.jar pileup2snp   -B  --variants  --min-var-freq $MINAF > $OSS.txt
    samtools mpileup -f $OS.fa $OS.bam -B -d $MAXDP | java -jar $HP_JDIR/VarScan.jar pileup2indel -B  --variants  --min-var-freq $MINAF | grep -v "^#" >> $OSS.txt

    cat $HP_SDIR/$M.vcf  > $OSS.vcf
    echo "##sample=$S" >> $OSS.vcf
    fa2Vcf.pl $OS.fa | grep -m 1 contig= >> $OSS.vcf
    fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OSS.vcf
    cat $OSS.txt | \
      perl -ane 'next if($.==1 or $F[2]=~/N/i or $F[3]=~/N/i); $DP=$F[4]+$F[5];$F[6]=~/(.+)%/;$AF=$1/100; print join "\t",($F[0],$F[1],".",$F[2],$F[-1],".",".",".","GT:DP:AF","0/1:$DP:$AF");print "\n";' | \
      perl -ane 'if($F[4]=~/\+(.+)/) {($F[4],$F[7])=("$F[3]$1","INDEL")} elsif($F[4]=~/\-(.+)/) {($F[3],$F[4],$F[7])=("$F[3]$1",$F[3],"INDEL")} print join "\t",@F; print "\n";' | uniqVcf.pl >> $OSS.vcf
  else
    echo "Unsuported SNV caller 2"
    exit 1
  fi
   rm  $OSR.* 
fi

if [ ! -s $OSS.00.vcf ] ; then
  # filter SNP
  cat $OSS.vcf | bcftools norm -m-any -f $OS.fa   | fix${M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header> $OSS.fix.vcf
  cat $OSS.fix.vcf | fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -rlen $HP_MTLEN -mfile $OS.max.vcf  | \
    filterVcf.pl -sample $S -source $M -header $HP_SDIR/$M.vcf  -depth $HP_DP | bedtools sort  -header > $OSS.00.vcf  

  annotateVcf.sh $OSS.00.vcf 
  intersectVcf.pl $OS.00.vcf $OS.max.vcf | cat - $OSS.00.vcf |  uniqVcf.pl | bedtools sort -header > $OSS.00.vcf.tmp ; mv $OSS.00.vcf.tmp $OSS.00.vcf
  #intersectVcf.pl $OS.00.vcf $OS.max.vcf | differenceVcf.pl - $OSS.00.vcf  | perl -ane 'if(/(.+):0\.\d+$/) { print "$1:1\n"} else { print }' | cat - $OSS.00.vcf | uniqVcf.pl | bedtools sort -header > $OSS.00.vcf.tmp # new(rna-seq)
fi

rm -f $OS.bam*  $OS.dict $OS.fa.fai
rm -f $O.fq $ON.score $O.bam* $OR.bam* #!!!

