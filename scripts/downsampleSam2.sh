#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16GB

P=4

#if [ ! -f bams+/$1.bam ] ; then
#  fastq-dump $1  --split-spot  -Z |  \
#    minimap2 -ax sr -t 8 /home/dpuiu1/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.mmi /dev/stdin | samtools view -f 0x2 -h -b > bams+/$1.bam
#fi

#if [ ! -s bams+/$1.cvg ]  ; then
#  cat bams+/$1.bam | bedtools bamtobed -cigar | grep ^chrM | bedtools genomecov -i - -g ~/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.fa.fai  -d | tee bams+/$1.cvg  | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > bams+/$1.cvg.stat
#fi


#if [ ! -s bams/$1.bam ] ; then
#  samtools fastq  bams+/$1.bam |  \
#    minimap2 -ax sr -t 8  ~/Homo_sapiens/hs38DH/hs38DH.mmi /dev/stdin  |\
#    perl -ane '$F[6]=$F[2] if($F[6] eq "="); print  if(/^\@/ or ($F[2] eq "chrM" or $F[2] eq "chr1" and 629084<$F[3] and $F[3]<634672 or $F[2] eq "chr17" and 22521208<$F[3] and $F[3]<22521639) and ($F[6] eq "chrM" or $F[6] eq "chr1" and 629084<$F[7] and $F[7]<634672 or $F[6] eq "chr17" and 22521208<$F[7] and $F[7]<22521639));' | \
#    samtools view -f 0x2 -h -b /dev/stdin | samtools sort > bams/$1.bam

#  samtools index bams/$1.bam
#  samtools idxstats bams/$1.bam > bams/$1.idxstats
#fi

if [ ! -s bams+/$1.cvg.stat ] ; then
  cat bams+/$1.bam | bedtools bamtobed -cigar | grep ^chrM | bedtools genomecov -i - -g ~/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.fa.fai  -d | tee bams+/$1.cvg  | \
    cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > bams+/$1.cvg.stat
fi


#if [ ! -s bams+/$1.dedup.filtered.bam ] ; then
  samtools view -h bams+/$1.bam | \
    samblaster -i /dev/stdin -o /dev/stdout --removeDups --addMateTags | \
    samtools view -b > bams+/$1.dedup.bam

  #  grep -v '\@'  | perl -ane '@F=split /\t/; $F[3]=($F[3]-$1)%16569 if($F[5]=~/^(\d+)[SH]/); $F[7]=($F[7]-$1)%16569 if(/MC:Z:(\d+)[SH]/);print join "\t",@F;' | \
  #  tee bams+/$1.dedup.adjust.sam | \
  #  count.pl -i 3    >  bams+/$1.dedup.adjust.start.count

  #join.pl -i 3 bams+/$1.dedup.adjust.sam bams+/$1.dedup.adjust.start.count 	| \
  #  join.pl -i 7 - bams+/$1.dedup.adjust.start.count |\
  #  perl -lane '@F=split /\t/; $last=pop @F; @F=($last,@F);  $last=pop @F; @F=($last,@F);  print join "\t",@F;'| \
  #  sort -k1,1n -k2,2n | \
  #  uniq.pl -i 5 -c 15 | cut -f 3| \
  #  samtools view -N /dev/stdin bams+/$1.bam  -b > bams+/$1.dedup.filtered.bam

  bedtools bamtobed -i  bams+/$1.dedup.bam -bedpe | \
    perl -lane '@F[1,2,4,5]=@F[4,5,1,2] if($F[1]>$F[4]); print join "\t",@F;' | \
    sort -k2,2n -k5,5n |  \
    perl -ane '$u=0; $s=0; $n=0; $m=$c[$F[1]]; foreach ($F[1]..$F[2],$F[4]..$F[5]) { $s+=$c[$_];$n++; $m=$c[$_] if($m<$c[$_])}  if($s/$n<2000 and $m<2200) { print;  foreach ($F[1]..$F[2],$F[4]..$F[5])  { $c[$_]++}}'  | \
    cut -f7 | \
    samtools view -N /dev/stdin bams+/$1.bam  -b > bams+/$1.dedup.filtered.bam

  #map functiomn
  #@array[$start_index .. $end_index] = map { $_ * 2 } @array[$start_index .. $end_index];

  cat bams+/$1.dedup.filtered.bam | bedtools bamtobed -cigar | grep ^chrM | bedtools genomecov -i - -g ~/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.fa.fai  -d | tee bams+/$1.dedup.filtered.cvg  | \
    cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > bams+/$1.dedup.filtered.cvg.stat
#fi


mkdir -p bams/
if [ ! -s bams/$1.bam ] ; then
  samtools fastq  bams+/$1.dedup.filtered.bam |  \
    minimap2 -ax sr -t 8  ~/Homo_sapiens/hs38DH/hs38DH.mmi /dev/stdin  |\
    perl -ane '$F[6]=$F[2] if($F[6] eq "="); print  if(/^\@/ or ($F[2] eq "chrM" or $F[2] eq "chr1" and 629084<$F[3] and $F[3]<634672 or $F[2] eq "chr17" and 22521208<$F[3] and $F[3]<22521639) and ($F[6] eq "chrM" or $F[6] eq "chr1" and 629084<$F[7] and $F[7]<634672 or $F[6] eq "chr17" and 22521208<$F[7] and $F[7]<22521639));' | \
    samtools view -f 0x2 -h -b /dev/stdin | samtools sort > bams/$1.bam

  samtools index bams/$1.bam
  samtools idxstats bams/$1.bam > bams/$1.idxstats

  cat bams/$1.bam | bedtools bamtobed -cigar | grep ^chrM | bedtools genomecov -i - -g ~/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.fa.fai  -d | tee bams/$1.cvg  | \
    cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > bams/$1.cvg.stat
fi





