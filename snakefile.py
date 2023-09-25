# Snakemake script for WGBS data analysis

import pandas as pd
import os
import argparse
import sys

configfile: "config.yaml"


def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_path(sample_df, wildcards, col):
    return sample_df.loc[[wildcards.sample], col].dropna().values[0]

samples_df = parse_samples(config["sample"])


OUTPUT_DIR = config["outdir"]
REFERENCE_DIR = config["ref"]
SAMPLES = samples_df['id'].tolist()
#print(samples_df.loc[["sample1","sample2"],"fq1"])

rule all:
    input:
        expand(f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.deduplicated.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/4.bismark_methylation/{{sample}}", sample=SAMPLES)


# Rule for trimming using Soapnuke
rule soapnuke:
    input:
        fastq_1 = lambda wildcards: get_path(samples_df, wildcards, 'fq1'),
        fastq_2 = lambda wildcards: get_path(samples_df, wildcards, 'fq2')
    output:
        trim_dr = directory(f"{OUTPUT_DIR}/1.trimmed/{{sample}}"),
        trimmed_fastq_1 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_1_trimmed.fq.gz",
        trimmed_fastq_2 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_2_trimmed.fq.gz"
        #trimmed_report_1 = f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_1_trimmed_fastqc.zip",
        #trimmed_report_2 = f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_2_trimmed_fastqc.zip"
    benchmark:
        "benchmarks/1.trimmed/{sample}.benchmark.txt"
    log:
        "logs/1.trimmed/{sample}.log"
    params:
        out_dir=directory(f"{OUTPUT_DIR}/1.trimmed/{{sample}}"),
        clean1=f"{{sample}}_1_trimmed.fq.gz",
        clean2=f"{{sample}}_2_trimmed.fq.gz"
    threads: config["threads"]["soapnuke"]
    shell:
        'SOAPnuke filter -l 30 -q 0.5 -1 {input.fastq_1} -2 {input.fastq_2} -C {params.clean1} -D {params.clean2} -o {params.out_dir} -T {threads}'


# Rule for alignment using Bismark
rule bismark_alignment:
    input:
        trimmed_fastq_1 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_1_trimmed.fq.gz",
        trimmed_fastq_2 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_2_trimmed.fq.gz"
    output:
        alin_dr = directory(f"{OUTPUT_DIR}/2.bismark_align/{{sample}}"),
        bam = f"{OUTPUT_DIR}/2.bismark_align/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.bam"
        #html = f"{OUTPUT_DIR}/3.bismark_align/{{sample}}/Testpaired_PE_report.html"
    benchmark:
        "benchmarks/2.bismark_align/{sample}.benchmark.txt"
    log:
        "logs/2.bismark_align/{sample}.log"
    threads: config["threads"]["bis_align"]
    shell:
        'bismark --bowtie2 --parallel {threads} --genome {REFERENCE_DIR} -o {output.alin_dr} -1 {input.trimmed_fastq_1} -2 {input.trimmed_fastq_2} 2> {log}'

rule bismark_duplicate:
    input:
        bam_file = f"{OUTPUT_DIR}/2.bismark_align/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.bam"
    output:
        dup_dr = directory(f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}"),
        duplicate = f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.deduplicated.bam"
    benchmark:
        "benchmarks/3.bismark_deduplicate/{sample}.benchmark.txt"
    log:
        "logs/3.bismark_deduplicate/{sample}.log"
    threads: config["threads"]["bis_dedul"]
    shell:
        'deduplicate_bismark -p --parallel {threads} --bam {input.bam_file} --output_dir {output.dup_dr} 2> {log}'

rule bismark_extract:
    input:
        bam_dedup_file = f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.deduplicated.bam"
    output:
        extr_dr = directory(f"{OUTPUT_DIR}/4.bismark_methylation/{{sample}}"),
        #txt = f"{OUTPUT_DIR}/5.bismark_methylation/{{sample}}/Testpaired_pe.deduplicated.nonCG_filtered.M-bias.txt"
    benchmark:
        "benchmarks/4.bismark_methylation/{sample}.benchmark.txt"
    log:
        "logs/4.bismark_methylation/{sample}.log"
    threads: config["threads"]["bis_call"]
    shell:
        'bismark_methylation_extractor --comprehensive --report -p --no_overlap --gzip --parallel {threads} --output {output.extr_dr} --genome_folder {REFERENCE_DIR} {input.bam_dedup_file} 2> {log}'

