#snakemake -s download_aws.smk --cores 2 --config tsv=Illumina.tsv
#snakemake -s download_aws.smk --cores 2 --config sample=HG02970 file=s3://human-pangenomics/.../HG02970.cram

import pandas as pd

RMT = "chrM"
MTLEN = 16569
RNUMT = "chr1:629084-634672 chr17:22521208-22521639"

if "tsv" in config:
    data = pd.read_csv(config["tsv"], sep="\t")

    SAMPLES = data["sample"].tolist()
    FILES = dict(zip(data["sample"], data["file"]))

elif "sample" in config and "file" in config:
    SAMPLES = [config["sample"]]
    FILES = {config["sample"]: config["file"]}

else:
    raise ValueError("Specify either tsv=FILE or both sample=SAMPLE and file=FILE")



def file(wildcards):
    return FILES[wildcards.sample]

rule ALL:
    input:
        expand("{sample}.MT.bam", sample=SAMPLES),
        expand("{sample}.MT.bam.bai", sample=SAMPLES),
        expand("{sample}.MT.idxstats", sample=SAMPLES),
        expand("{sample}.MT.count", sample=SAMPLES)

rule DOWNLOAD:
    output:
        temp("{sample}.cram")
    params:
        remote=file
    shell:
        """
        s5cmd cp {params.remote} {output}
        """

rule INDEX:
    input:
        "{sample}.cram"
    output:
        temp("{sample}.cram.crai")
    threads: 2
    shell:
        """
        samtools index -@ {threads} {input}
        """

rule COUNT:
    input:
        cram="{sample}.cram",
        crai="{sample}.cram.crai"
    output:
        "{sample}.MT.count"
    threads: 2
    shell:
        """
        samtools idxstats -@ {threads} {input.cram} |
        ./idxstats2count.pl \
            -sample {wildcards.sample} \
            -chrM {RMT} \
            > {output}
        """

rule EXTRACT:
    input:
        cram="{sample}.cram",
        crai="{sample}.cram.crai"
    output:
        "{sample}.MT.bam"
    shell:
        """
        samtools view -h {input.cram} \
            {RMT}:1-{MTLEN} {RNUMT} \
            -F 0x90C |
            ./filterSam.pl \
            {RMT}:1-{MTLEN} {RNUMT} |
            samtools view -b > {output}
        """

rule INDEX_EXTRACTED:
    input:
        ancient("{sample}.MT.bam")
    output:
        "{sample}.MT.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule IDXSTATS:
    input:
        bam=ancient("{sample}.MT.bam"),
	bai=ancient("{sample}.MT.bam.bai")
    output:
        "{sample}.MT.idxstats"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """
