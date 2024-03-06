#!/bin/bash

export S=$1

###########################


mkdir -p out/$1

if [ ! -s out/$1/$1.varscan.txt ] ; then
  #--output-vcf   --min-reads2 4  --variants 1
  if [ ! -s out/$1/$1.varscan.snp.txt ] ; then
    samtools mpileup -f $HP_RDIR/hs38DH.fa bams/$1.bam -r chrM -B -d 2000  | tee \
      >(java -jar $HP_JDIR/VarScan.v2.4.6.jar pileup2snp   -q 0 -B  --variants  > out/$1/$1.varscan.snp.txt) \
      >(java -jar $HP_JDIR/VarScan.v2.4.6.jar pileup2indel -q 0 -B  --variants  > out/$1/$1.varscan.indel.txt) > /dev/null
  fi

  cat out/$1/$1.varscan.snp.txt out/$1/$1.varscan.indel.txt | sort -k2,2n -u >  out/$1/$1.varscan.txt
  rm  out/$1/$1.varscan.snp.txt out/$1/$1.varscan.indel.txt
fi

cat out/$1/$1.varscan.txt  | uniq | \
  perl -ane 'next if($.==1 or $F[2]=~/N/i or $F[3]=~/N/i); $DP=$F[4]+$F[5];$F[6]=~/(.+)%/;$AF=$1/100; print join "\t",($F[0],$F[1],".",$F[2],$F[-1],".",".","SM=$ENV{S}","GT:DP:AF","0/1:$DP:$AF");print "\n";' | \
  perl -ane 'if($F[4]=~/\+(.+)/) {($F[4],$F[5])=("$F[3]$1","INDEL")} elsif($F[4]=~/\-(.+)/) {($F[3],$F[4],$F[5])=("$F[3]$1",$F[3],"INDEL")} print join "\t",@F; print "\n";' | uniqVcf.pl > out/$1/$1.varscan.vcf

#cat out/*/*vcf | concat2cat.pl | sort -k2,2n > out/varscan.00.vcf
