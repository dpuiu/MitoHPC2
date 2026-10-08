# HPRC Data Download Scripts

These scripts are used to download **HPRC sequencing samples from AWS S3** for use with MitoHPC2.  

The scripts support: Illumina, PacBio (HiFi) , Oxford Nanopore (ONT)  

Both **indexed and unindexed** files can be downloaded.  

The processing requirements differ by sequencing technology:

* **Illumina:** Downloaded CRAM files must be **aligned/sorted** before mitochondrial reads can be extracted.
* **HiFi/ONT:** Downloaded BAM files do **not** need to be aligned before mitochondrial reads are extracted.

For **Illumina CRAM files**, the complete **hg38 reference genome** is required.

---

## Prerequisites

Before running the HPRC scripts, the MitoHPC2 environment and reference files must be installed, initialized and tested

## 0. Update and Initialize MitoHPC2

```bash
cd MitoHPC2
git pull

export HP_SDIR="$PWD/scripts"
. "$HP_SDIR/init.sh"

printenv | grep '^HP_'
```

## 1. Create the Conda environment

The Conda environment replaces the previous `install_prerequisites.sh` software installation procedure.

Run this **once**:

```bash
conda init
conda env create -f $HP_SDIR/scripts/install_prerequisites.yaml
```

The `mitohpc2` environment provides the main tools required by the scripts, including: s5cmd, samtools, bwa, minimap2 ...

## 2. Activate the Conda environment

```bash
conda activate mitohpc2
```

## 3. Install/configure reference files

Reference files are installed separately. Run:

```bash
$HP_SDIR/scripts/install_prerequisites.ref.sh
```

## 4. Initialize MitoHPC2 Long Read (if necessary)

Set the required `HP_` variables:

```bash
$HP_SDIR/init.hifi.sh # or
$HP_SDIR/init.ont.sh
```

## 5. Check the installation

Verify that the software and reference files are correctly configured:

```bash
$HP_SDIR/checkInstall.sh
$HP_SDIR/checkInstall.ref.sh
```
---

## 6. Illumina

Create a directory for the Illumina data:

```bash
mkdir Illumina
cd Illumina
```

### Download one Illumina sample

To see the command used for the first sample:

```bash
head -1 $HP_SDIR/../scriptsHPRC/download_aws_illumina.all.sh
```

For example:

```text
download_aws_illumina.sh HG00097 s3://human-pangenomics/submissions/59C50DDF-5FAF-4841-AC3E-6C02D636C57F--Y4_1000G_DATA/HG00097.final.cram
```

Run the first download:

```bash
head -1 $HP_SDIR/../scriptsHPRC/download_aws_illumina.all.sh | bash
```

### Download all Illumina samples

Once the test download works:

```bash
$HP_SDIR/../scriptsHPRC/download_aws_illumina.all.sh
```

Downloading and processing all HPRC samples requires substantial disk space. Check available storage before starting a large download.

---

## 7. HiFi

*HiFi download and processing instructions go here.*

## 8. ONT

*ONT download and processing instructions go here.*
