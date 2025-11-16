#!/bin/bash

S=$1
s=${S,,}

#identity:mean,max,stdev  (99,100,0.5)
##params="--length 10000,25000 --identity ... --junk_reads 0.001 --random_reads 0.001 --chimeras 0.001 --glitches 0"
#params="--error_model pacbio2021 --qscore_model pacbio2021 --identity 99,100,0.5"   # used

#ref="/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/RefSeq/$S.fa"
#badread simulate --reference $ref --quantity "1600x" $params 2>$S.log 1>$S.fq

#ref="/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/RefSeq/$s.fa"
#badread simulate --reference $ref --quantity "400x" $params 2>$S.log 1>>$S.fq

#samtools faidx $S.fq

minimap2 ~/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.fa $S.fq -a | samtools sort > ../bams/$S.bam
samtools index ../bams/$S.bam

#bwa mem ~/Homo_sapiens_mito/MitoHPC2/RefSeq/chrM.fa $S.fq | samtools sort > ../bams2/$S.bam
#samtools index ../bams2/$S.bam

