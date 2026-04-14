# DISCO
A Snakemake-based end-to-end pipeline for processing Oxford Nanopore Technology (ONT) long-read RNA sequencing data, from raw reads through isoform discovery and downstream biological analysis.

This pipeline automates strand-aware alignment, quality control, transcript reconstruction, and downstream analysis of nanopore direct RNA or cDNA sequencing data. It is designed for multi-sample experiments and produces isoform-resolved expression quantification alongside comparative QC metrics, differential expression results, pathway enrichment, and long non-coding RNA biomarker discovery.
Pipeline stages
1. Strand correction — Restrander
Corrects strand orientation of raw FASTQ reads using a user-supplied configuration. Outputs restranded FASTQ files and per-sample JSON statistics.
2. Splice-aware alignment — Minimap2 + samtools
Aligns restranded reads to a reference genome using minimap2 in splice-aware mode (-ax splice -uf -k14). Alignments are coordinate-sorted and indexed with samtools.
3. Quality control — NanoComp
Runs NanoComp across all samples simultaneously to generate comparative read-level QC metrics, including read length distribution, base quality, and sequencing yield. Outputs an interactive HTML report and summary statistics.
4. Isoform quantification — IsoQuant
Reconstructs transcript models and quantifies gene and isoform expression per sample using IsoQuant with a complete gene database. Produces novel GTF annotations and SQANTI-compatible output for downstream filtering.
5. Downstream analysis

Differential expression — DESeq2 or edgeR on IsoQuant count matrices, visualised as volcano plots (log₂FC vs −log₁₀ adjusted p-value) to identify up- and down-regulated isoforms across conditions.
Pathway enrichment — Over-representation analysis and GSEA using clusterProfiler against GO, KEGG, and Reactome databases at isoform resolution.
lncRNA biomarker discovery — Coding-potential filtering of novel transcripts with FEELnc and CPAT, co-expression network analysis via WGCNA, and biomarker scoring against phenotype of interest.

Requirements

Snakemake ≥ 7.x
minimap2
samtools
NanoComp
IsoQuant (conda env: isoquant)
Restrander binary
R with DESeq2 / edgeR, clusterProfiler, WGCNA
Python ≥ 3.8
