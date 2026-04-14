# ── Snakefile ────────────────────────────────────────────────────────────────
# Nanopore long-read pipeline: restrander → minimap2 → samtools index → IsoQuant
#
# Usage:
#   snakemake --configfile config.yaml  --rerun-incomplete
# ─────────────────────────────────────────────────────────────────────────────

import os
from glob import glob

# ── config ────────────────────────────────────────────────────────────────────
configfile: "config.yaml"

FASTA      = config["fasta"]
GTF        = config["gtf"]
REST_BIN   = config["restrander_bin"]
REST_CFG   = config["restrander_config"]
INPUT_DIR  = config["input_dir"]
OUT        = config["output_dir"]

T_MM2      = config["threads_minimap2"]
T_SORT     = config["threads_sort"]
T_ISO      = config["threads_isoquant"]

# ── sample discovery ──────────────────────────────────────────────────────────
SAMPLES, = glob_wildcards(os.path.join(INPUT_DIR, "{sample}.fastq.gz"))

if not SAMPLES:
    raise ValueError(f"No *.fastq.gz files found in {INPUT_DIR}")

#directory layout 
RESTRANDED = os.path.join(OUT, "restranded")
STATS      = os.path.join(OUT, "stats")
BAM        = os.path.join(OUT, "bam")
ISOQUANT   = os.path.join(OUT, "isoquant")
LOGS       = os.path.join(OUT, "logs")

# ── target rule
rule all:
    input:
        expand(
            os.path.join(ISOQUANT, "{sample}", "{sample}", "{sample}.transcript_models.gtf"),
            sample=SAMPLES
        )

#rule 1: restrander 
rule restrander:
    input:
        fq = os.path.join(INPUT_DIR, "{sample}.fastq.gz")
    output:
        fq   = os.path.join(RESTRANDED, "{sample}_restranded.fastq.gz"),
        stats = os.path.join(STATS,      "{sample}_stats.json")
    log:
        os.path.join(LOGS, "restrander", "{sample}.log")
    shell:
        """
        mkdir -p {RESTRANDED} {STATS} $(dirname {log})
        {REST_BIN} {input.fq} {output.fq} {REST_CFG} \
            > {output.stats} \
            2> {log}
        """

#rule 2: minimap2 + samtools sort 
rule minimap2:
    input:
        fq   = os.path.join(RESTRANDED, "{sample}_restranded.fastq.gz"),
        fasta = FASTA
    output:
        bam = os.path.join(BAM, "{sample}.minimap2.bam")
    threads: T_MM2 + T_SORT
    log:
        os.path.join(LOGS, "minimap2", "{sample}.log")
    shell:
        """
        mkdir -p {BAM} $(dirname {log})
        minimap2 \
            -ax splice -uf -k14 \
            --secondary=no \
            -t {T_MM2} \
            {input.fasta} {input.fq} \
        2> {log} \
        | /home/mohamad/miniconda3/envs/rnaseq_align/bin/samtools sort \
            -@ {T_SORT} \
            -T {BAM}/{wildcards.sample}_tmp \
            -o {output.bam}
        """

#rule 3: samtools index
rule samtools_index:
    input:
        bam = os.path.join(BAM, "{sample}.minimap2.bam")
    output:
        bai = os.path.join(BAM, "{sample}.minimap2.bam.bai")
    log:
        os.path.join(LOGS, "samtools_index", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        /home/mohamad/miniconda3/envs/rnaseq_align/bin/samtools index {input.bam} 2> {log}
        """

#rule 4: IsoQuant
rule isoquant:
    input:
        bam  = os.path.join(BAM, "{sample}.minimap2.bam"),
        bai  = os.path.join(BAM, "{sample}.minimap2.bam.bai"),
        fasta = FASTA,
        gtf   = GTF
    output:
        gtf = os.path.join(ISOQUANT, "{sample}", "{sample}", "{sample}.transcript_models.gtf")
    threads: T_ISO
    log:
        os.path.join(LOGS, "isoquant", "{sample}.log")
    shell:
        """
        mkdir -p {ISOQUANT}/{wildcards.sample} $(dirname {log})
        /home/mohamad/miniconda3/envs/isoquant/bin/python /home/mohamad/miniconda3/envs/isoquant/bin/isoquant.py \
            --reference  {input.fasta} \
            --genedb     {input.gtf} \
            --bam        {input.bam} \
            --data_type  nanopore \
            --output     {ISOQUANT}/{wildcards.sample} \
            --threads    {threads} \
            --complete_genedb \
            --sqanti_output \
            --prefix     {wildcards.sample} \
            > {log} 2>&1
        """
