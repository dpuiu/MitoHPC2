#!/usr/bin/env bash 
#set -e

#if [ -z $HP_HDIR ] ; then echo "Variable HP_HDIR not defined. Make sure you followed the SETUP ENVIRONMENT instructions" ;  fi

##############################################################################################################

# Program that setups the environmnet

# Variable HP_HDIR must be pre-set !!!

##############################################################################################################
#DIRECTORY PATHS

#export HP_HDIR=`readlink -f $HP_SDIR/..`	#HP home directory
export HP_SDIR=$HP_HDIR/scripts/
export HP_BDIR=$HP_HDIR/bin/			#bin directory
export HP_JDIR=$HP_HDIR/java/			#java directory
export HP_RDIR=$HP_HDIR/RefSeq/			#Human reference directory

###############################################################
#SOFTWARE PATH

export PATH=$HP_SDIR:$HP_BDIR:$PATH

################################################################
#ALIGNMNET REFERENCE

#hs38DH(default)
export HP_RNAME=hs38DH
export HP_RMT=chrM
export HP_RURL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
export HP_RFILE=$HP_RDIR/$HP_RNAME

################################################################
#GENOME REFERENCES

export HP_O=Human		 # organism: Human, Mouse...
export HP_MT=chrM                # chrM, rCRS or RSRS, FASTA file available under $HP_RDIR
export HP_MTLEN=16569
export HP_NUMT=NUMT

################################################################
#PARAMETERS

export HP_TYPE="ont"            # hifi/ont ("minimap2 -x map-$HP_TYPE")
export HP_L=4000                 # maximum number of alignments; increased from 2000 to 4000
export HP_M=bcftools 	         # SNV caller: bcftools

export HP_CN=1                   # do compute mtDNA copy number
export HP_I=2		         # number of SNV iterations : 0: compute read counts,mtDNA-CN; 1:1 iteration (mutect2,mutserve) ;  2:2 iterations (mutect2)
export HP_T1=03                  # heteroplasmy tholds
export HP_T2=05
export HP_T3=10
export HP_DP=10    		 # minimum coverage: Ex 25
export HP_V=                     # SV caller: gridss

export HP_FRULE="perl -ane 'print unless(/Homopolymer/ and /:0\.[01234]\d*$/);' |  bcftools filter -e 'DP<$HP_DP'"   # filter rule (or just "tee")

################################################################
#COMPUTATION 
export HP_P=1				       		            # number of processors
export HP_MM="3G"                                                   # maximum memory
export HP_JOPT="-Xms$HP_MM -Xmx$HP_MM -XX:ParallelGCThreads=$HP_P"  # JAVA options

#export HP_SH="bash" ;                                                                        export HP_SHS="$HP_SH"                     # bash
export HP_SH="sbatch -J HP_$$ --cpus-per-task=$HP_P --nodes=1 --mem=$HP_MM --time=20:00" ;  export HP_SHS="$HP_SH -d singleton"        # SLURM job scheduling
#export HP_SH="qsub -V -N HP_$$ -l mem_free=$HP_MM,h_vmem=$HP_MM -pe local $HP_P -cwd" ;     export HP_SHS="$HP_SH -hold_jid HP_$$"     # SGE job scheduling

################################################################
#INPUT/OUTPUT DIRS

PWD=`pwd -P`
#export HP_FDIR=$PWD/fastq/      # fastq input file directory ; .fq or .fq.gz file extension
export HP_ADIR=$PWD/bams/	# bams or crams input file directory
export HP_ODIR=$PWD/out/        # output dir
export HP_IN=$PWD/in.txt        # input file to be generated


################################################################
# PRINT ENV

printenv | egrep '^HP_|^PATH=' | sed 's|=|="|' | sed 's|$|"|' | sort > env.txt

################################################################
#IN FILE

if [ -d $HP_ADIR ] ; then
  if [ ! -s $HP_IN ] ; then
    find $HP_ADIR/  -name "*.bam" -o -name "*.cram" -readable | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN
  fi
fi

