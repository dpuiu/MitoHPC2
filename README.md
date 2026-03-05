# MitoHPC2: Mitochondrial High Performance Caller v2

A pipeline for detecting mitochondrial **homoplasmies** and **heteroplasmies** from sequencing data.  

Capabilities:
- Supports short and long reads
- Multiple SNV callers
- Human and mouse genomes
- Multiple reference versions
- Multiple heteroplasmy thresholds
- Calls maternal haplogroup
- Detects contamination
- Processes multiple samples in parallel and combines results
- Optimized for low CPU, memory, and runtime usage

---

## Citing

If you use this pipeline, please cite:

Battle et al., *NAR 2022*: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/)

---

## Reference

For additinal information please check 
(MitoHPC v1)[https://github.com/dpuiu/MitoHPC]

---

## Prerequisites

The pipeline requires the following tools:


- [bwa 0.7.17](https://github.com/lh3/bwa/releases), 
[minimap2 2.28](https://github.com/lh3/minimap2/releases), 
[htslib 1.21](https://github.com/samtools/htslib/releases), 
[samtools 1.21](https://github.com/samtools/samtools/releases), 
[bcftools 1.21](https://github.com/samtools/bcftools/releases), 
[bedtools 2.31.1](https://github.com/arq5x/bedtools2/releases), 
[fastp 0.24.0](http://opengene.org/fastp/fastp), 
[samblaster 0.1.26](https://github.com/GregoryFaust/samblaster/releases), 
[GATK Mutect2 4.6.0.0](https://github.com/broadinstitute/gatk/releases/), 
[mutserve 2.0.0-rc15](https://github.com/seppinho/mutserve/releases), 
[freebayes 1.3.6](https://github.com/freebayes/freebayes/releases), 
[VarScan 2.4.6](https://github.com/dkoboldt/varscan/releases), 
[Clair3 1.2.0](https://github.com/HKU-BAL/Clair3/releases), 
[ClairS-TO 0.4.2](https://github.com/HKU-BAL/ClairS-TO), 
[DeepVariant 1.10.0-beta](https://github.com/google/deepvariant), 
[DeepSomatic 1.9.0](https://github.com/google/deepsomatic/), 
[GRIDSS 2.13.2](https://github.com/PapenfussLab/gridss/releases), 
[Haplogrep 2.4.0](https://github.com/seppinho/haplogrep-cmd/releases), 
[Haplocheck 1.3.3](https://github.com/genepi/haplocheck/releases), 
[plink2 2.0](https://www.cog-genomics.org/plink/2.0/)
   
- [R](https://cran.r-project.org/src/base/), 
[Apptainer](https://github.com/apptainer/apptainer.git)

---

## Reference Genomes

- Default: hg38
- Setting available in  scripts/init.sh

----

## Input Samples

The pipeline expects pre-aligned, coordinate sorted reads:

- **File formats**: `.bam` or `.cram` 
- **Sequencing platforms**: Illumina paired-end, PacBio HiFi, ONT  
- **Index files**: `.bai` or `.crai`  
- **Indexstats files**: `.idxstats`: contain the number of reads aligned to each chromosome  

---

## Installation

### 1. Download the Pipeline

```bash
git clone https://github.com/dpuiu/MitoHPC2.git
```

### 2. Update the Pipeline (optional)

```bash
cd MitoHPC2/
git pull
# or
git checkout .
```

### 3. Set up the Environment (important)

```bash
cd MitoHPC2/scripts
export HP_SDIR=`pwd`                           # set script directory
echo "export HP_SDIR=`pwd`" >> ~/.bashrc       # add variable to .bashrc
. ./init.sh                                    # or init.{hs38DH,hg19,mm39}.sh for different references
```

### 4. Install System Prerequisites (optional)

```bash
sudo $HP_SDIR/install_sysprerequisites.sh      # installs perl, python, java, wget, etc.
```

### 5. Install Pipeline Prerequisites

```bash
$HP_SDIR/install_prerequisites.sh             # installs bwa, samtools, bedtools, etc.
# or force latest versions
$HP_SDIR/install_prerequisites.sh -f
```

### 6. Install Long-Read SNV Callers (optional)

Download pre-built Singularity images to `$HP_BDIR`:

```bash
curl -s ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/MitoHPC2/bin/
# example files:
# clair3_v1.2.0.sif
# clairs-to_v0.4.2.sif
# deepsomatic_1.9.0.sif
# deepvariant_1.10.0-beta.sif
```

Convert to sandbox for faster execution:

```bash
singularity build --sandbox ~/clairs-to_sandbox/   $HP_BDIR/clairs-to.sif
singularity build --sandbox ~/deepsomatic_sandbox/ $HP_BDIR/deepsomatic.sif

# Check size
du -hs ~/clairs-to_sandbox/ ~/deepsomatic_sandbox/
```

### 7. Check Installation

```bash
$HP_SDIR/checkInstall.sh
cat checkInstall.log
```

"Success message!" expected

---

## Custom Annotation

- The pipeline uses `MitoHPC2/RefSeq/*.{vcf,bed}.gz` for annotation.  
- Copy any custom VCF/BED files to `MitoHPC2/RefSeq/`.  
- Ensure files are bgzipped and indexed.

---

## Running the Pipeline on Illumina Data

### SETUP ENVIRONMENT

```
# go to your working directory 
# copy over the init.sh file

cp -i $HP_SDIR/init.sh .                      # copy init to working dir.

# view/edit init.sh file
nano init.sh                                  # check/edit local init file
    ...
    export HP_O=Human                         # organism
    export HP_MT=chrM
    export HP_MTLEN=16569
    export HP_RNAME=hs38DH                    # reference assembly used for alignment (Ex: hs38DH, hg19)   
    export HP_RMT=chrM                        # MT sequence name
    export HP_RNUMT="chr1:629084-634672 ..."  # NUMT regions to be considered     
    export HP_RURL=https://hgdownload..."     # URL location (if needed to be downloaded)
    export HP_IN=$PWD/in.txt                  # input TSV file, contains sample names, file paths (automatically generated by runnng ". ./init.sh")
    export HP_ADIR=$PWD/bams/                 # pre-existing alignment dir; contains .bam, .bai, [.idxstats] files, or
    export HP_ADIR=$PWD/crams/                #                                       .cram, .crai, [.idxstats] files
    export HP_ODIR=$PWD/out/                  # output dir  
    export HP_L=222000                        # Use at most 222K reads (random subsampling; ~2K cvg of MT for 150bp reads); if empty, use all the reads
    export HP_M=mutect2                       # SNV caller: mutect2, mutserve, or freebayes
    export HP_I=2                             # number of SNV calling iterations
                                              #  1: one, uses rCRS or RSRS as reference 
                                              #  2: two, the 2nd one used the sample consensus 
                                              #     computed in the 1st iteration; higher accuracy 
    export HP_T1=03                           # heteroplasmy thold 1: 3%
    export HP_T2=05                           # heteroplasmy thold 2: 5%
    export HP_T3=10                           # heteroplasmy thold 3: 10%
    export HP_SH=bash                         # job scheduling: bash, qsub,sbatch, ..

```

### Single SNV Caller

```bash
cp $HP_SDIR/init.sh .
nano ./init.sh             # edit SNV caller if needed
# e.g., export HP_M=mutect2
. ./init.sh
$HP_SDIR/run.sh | tee run.all.sh | bash
```

### Multiple SNV Callers

```bash
cp $HP_SDIR/init3.sh .
nano ./init3.sh            # set HP_M1, HP_M2, HP_M3
. ./init3.sh
$HP_SDIR/run3.sh | tee run.all.sh | bash
```

## PacBio HiFi Data

```bash
cp $HP_SDIR/init.hifi.sh .
nano ./init.hifi.sh        # set HP_M, HP_PLATFORM, HP_MODELTYPE
. ./init.hifi.sh
$HP_SDIR/run.lr.sh | tee run.all.sh | bash
```

## ONT Data

```bash
cp $HP_SDIR/init.ont.sh .
nano ./init.ont.sh         # set HP_M, HP_PLATFORM, HP_MODELTYPE
. ./init.ont.sh
$HP_SDIR/run.lr.sh | tee run.all.sh | bash
```

---

## Output

```bash
ls $HP_ODIR/.*
```

---

## Evaluation

Compare results to Illumina mutect2 (T=10):

```bash
$HP_SDIR/eval.sh Illumina/out/mutect2.10.concat.vcf $HP_M.10.concat.vcf | uniq.pl | column -t > $HP_M.10.eval
```

---

## Examples

### HPRC Samples

- [FTP download](ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/MitoHPC2/bin/examples/HPRC/)

### Simulated Samples

- 30 simulated samples (haplogroups A-Z)  
- 39 HPRC samples ([sequence info](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/summary.txt))  

#### HiFi10 & ONT10 Evaluations

- Multiple SNV caller evaluations: `clairs-to`, `deepsomatic`, `varscan`  
- Links provided for different heteroplasmy thresholds (short-read / long-read)

##### HiFi10

* clairs-to:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.10.10.eval),
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.05.10.eval),
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.10.05.eval),
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.05.05.eval),
[10.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.10.00.eval),
[05.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/clairs-to.05.00.eval)

* deepsomatic:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.10.10.eval),
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.05.10.eval),
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.10.05.eval),
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.05.05.eval),
[10.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.10.00.eval),
[05.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/deepsomatic.05.00.eval)

* varscan:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/varscan.10.10.eval),
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/varscan.05.10.eval),
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/varscan.10.05.eval),
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/varscan.05.05.eval),
[10.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/varscan.10.00.eval),
[05.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10/out/varscan.05.00.eval)

##### ONT10

* clairs-to:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.10.10.eval),
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.05.10.eval),
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.10.05.eval),
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.05.05.eval),
[10.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.10.00.eval),
[05.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/clairs-to.05.00.eval)

* deepsomatic:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.10.10.eval),
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.05.10.eval),
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.10.05.eval),
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.05.05.eval),
[10.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.10.00.eval),
[05.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/deepsomatic.05.00.eval)

* varscan:
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/varscan.10.10.eval),
[05.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/varscan.05.10.eval),
[10.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/varscan.10.05.eval),
[05.05](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/varscan.05.05.eval),
[10.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/varscan.10.00.eval),
[05.00](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10/out/varscan.05.00.eval)


### 39 HPRC bases Simulated samples 

#### HiFi10Sim

* 500x cvg of Illumina paired end reads (2*150bp) were simulated using the best HiFi10 raeds (avg Q30)
* commands:
  wgsim -N 22200 -1 150 -2 150 -d 300 -e 0 -r 0 -R 0
* There reads were processed using the MitoHPC2 pipeline ("mutect2" SNV caller)
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/HiFi10Sim/out/mutect2.10.10.eval)

#### ONT10Sim 

* 500x cvg of Illumina paired end reads (2*150bp) were simulated using the best ONT10 raeds (avg Q20)
* commands:
  wgsim -N 22200 -1 150 -2 150 -d 300 -e 0 -r 0 -R 0
* There reads were processed using the MitoHPC2 pipeline ("mutect2" SNV caller)
[10.10](https://github.com/dpuiu/MitoHPC2/blob/main/examples/HPRC/ONT10Sim/out/mutect2.10.10.eval)


Note: 1st number: short read heteroplasmy thold;  2nd: long read heteroplasmy thold


