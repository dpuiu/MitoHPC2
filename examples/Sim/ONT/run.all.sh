#!/usr/bin/env bash
set -eux

export HP_SDIR=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts/
export HP_BDIR=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2//bin/
export HP_JDIR=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2//java/
export HP_RDIR=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2//RefSeq/
export HP_ODIR=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out/
export HP_RNAME=hs38DH
export HP_RMT=chrM
export HP_RNUMT=""
export HP_RCOUNT=
export HP_RURL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
export HP_O=Human
export HP_MT=chrM
export HP_MTR=chrMR
export HP_MTC=chrMC
export HP_MTLEN=16569
export HP_NUMT=NUMT
export HP_CN=1
export HP_E=
export HP_L=4000
export HP_IN=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/in.txt
export HP_ODIR=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out/
export HP_M=bcftools
export HP_I=2
export HP_T1=05
export HP_T2=10
export HP_T3=20
export HP_V=
export HP_DP=10
export HP_FRULE="perl -ane 'print unless(/Homopolymer/ and /:0\.[01234]\d*$/);' |  bcftools filter -e 'DP<10'"
export HP_FOPT=""
export HP_DOPT=""
export HP_GOPT=""
export HP_JOPT="-Xms3G -Xmx3G -XX:ParallelGCThreads=1"
export HP_MM="3G"
export HP_P=1
export PATH=/data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts/:/data/darking1/projects/mito/Heteroplasmy/MitoHPC2//bin/:$PATH

sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_A	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_A.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_A/sample_A
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_B	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_B.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_B/sample_B
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_C	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_C.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_C/sample_C
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_D	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_D.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_D/sample_D
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_E	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_E.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_E/sample_E
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_F	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_F.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_F/sample_F
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_G	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_G.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_G/sample_G
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_HV	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_HV.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_HV/sample_HV
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_H	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_H.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_H/sample_H
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_I	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_I.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_I/sample_I
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_J	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_J.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_J/sample_J
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_K	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_K.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_K/sample_K
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_L0	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_L0.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_L0/sample_L0
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_L1	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_L1.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_L1/sample_L1
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_L2	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_L2.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_L2/sample_L2
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_L3	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_L3.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_L3/sample_L3
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_L4	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_L4.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_L4/sample_L4
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_L5	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_L5.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_L5/sample_L5
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_M	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_M.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_M/sample_M
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_N	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_N.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_N/sample_N
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_P	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_P.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_P/sample_P
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_R	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_R.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_R/sample_R
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_S	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_S.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_S/sample_S
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_T	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_T.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_T/sample_T
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_U	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_U.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_U/sample_U
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_V	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_V.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_V/sample_V
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_W	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_W.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_W/sample_W
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_X	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_X.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_X/sample_X
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_Y	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_Y.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_Y/sample_Y
sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//filter.lr.sh sample_Z	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/bams//sample_Z.bam	/data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out//sample_Z/sample_Z

sbatch -J HP_1458225 --cpus-per-task=1 --nodes=1 --mem=3G --time=20:00 -d singleton /data/darking1/projects/mito/Heteroplasmy/MitoHPC2//scripts//getSummary.sh /data/darking1/projects/mito/Heteroplasmy/MitoHPC2/examples/Sim/ONT/out/
