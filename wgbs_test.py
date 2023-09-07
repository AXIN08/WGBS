# Snakemake script for WGBS data analysis

# Define input and output directories
INPUT_DIR = "/mnt/microbio_research/Raw_data/MagIC/MOMmy/WGBS/F23A430000203_HOMphtnH_20230725/soapnuke/clean"
OUTPUT_DIR = "/home/liuxin/data_pasteur/12_epigenome/Other_samples"
REFERENCE_DIR = "/home/liuxin/data_pasteur/12_epigenome/CHM13_ref"

# Define the samples
(DIR,SAMPLES) = glob_wildcards(f"{INPUT_DIR}/{{DIR}}/{{sample}}_1.fq.gz")  # Replace with your sample names
#SAMPLE_NAMES = [sample.sample for sample in SAMPLES]

# Define the workflow
rule all:
    input:
        #expand(f"{OUTPUT_DIR}/7.bismark_methylation/{{sample}}.CpG_context_test_dataset_bismark_bt2.txt.gz", sample=SAMPLES),
        expand(f"{INPUT_DIR}/{{sample}}/{{sample}}_1.fq.gz", sample=SAMPLES),
        expand(f"{INPUT_DIR}/{{sample}}/{{sample}}_2.fq.gz", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/1.fastqc/{{sample}}/{{sample}}_1_fastqc.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_1_val_1.fq.gz", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_2_val_2.fq.gz", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/3.bismark_align/{{sample}}/Testpaired_pe.bam", sample=SAMPLES),
        #expand(f"{OUTPUT_DIR}/3.bismark_align/{{sample}}", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/4.bismark_deduplicate/{{sample}}/Testpaired_pe.deduplicated.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/5.bismark_methylation/{{sample}}", sample=SAMPLES)
        #expand(f"{OUTPUT_DIR}/5.bismark_methylation/{{sample}}/CpG_context_test_dataset_bismark_bt2.txt.gz", sample=SAMPLES)


# Rule for quality control using FastQC
rule fastqc_quality_control:
    input:
        fastq_1 = f"{INPUT_DIR}/{{sample}}/{{sample}}_1.fq.gz",
        fastq_2 = f"{INPUT_DIR}/{{sample}}/{{sample}}_2.fq.gz"
    output:
        qc_dir = directory(f"{OUTPUT_DIR}/1.fastqc/{{sample}}"),
        html_report = f"{OUTPUT_DIR}/1.fastqc/{{sample}}/{{sample}}_1_fastqc.html",
        zip_report = f"{OUTPUT_DIR}/1.fastqc/{{sample}}/{{sample}}_1_fastqc.zip"
    benchmark:
        "2.benchmarks/{{sample}}.fastqc.benchmark.txt"
    shell:
        "fastqc -o {output.qc_dir} {input.fastq_1} {input.fastq_2}"

# Rule for trimming using Trim Galore
rule trim_galore:
    input:
        fastq_1 = f"{INPUT_DIR}/{{sample}}/{{sample}}_1.fq.gz",
        fastq_2 = f"{INPUT_DIR}/{{sample}}/{{sample}}_2.fq.gz"
    output:
        trim_dr = directory(f"{OUTPUT_DIR}/2.trimmed/{{sample}}"),
        trimmed_fastq_1 = f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_1_val_1.fq.gz",
        trimmed_fastq_2 = f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_2_val_2.fq.gz"
    benchmark:
        "2.benchmarks/{{sample}}.trim.benchmark.txt"
    shell:
        "trim_galore --phred33 -j 5 -q 30 --stringency 3 --fastqc --paired -o {output.trim_dr} {input.fastq_1} {input.fastq_2}"


# Rule for alignment using Bismark
rule bismark_alignment:
    input:
        trimmed_fastq_1 = f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_1_val_1.fq.gz",
        trimmed_fastq_2 = f"{OUTPUT_DIR}/2.trimmed/{{sample}}/{{sample}}_2_val_2.fq.gz"
    output:
        alin_dr = directory(f"{OUTPUT_DIR}/3.bismark_align/{{sample}}"),
        bam = f"{OUTPUT_DIR}/3.bismark_align/{{sample}}/Testpaired_pe.bam",
        #html = f"{OUTPUT_DIR}/3.bismark_align/{{sample}}/Testpaired_PE_report.html"
    benchmark:
        "2.benchmarks/{{sample}}.bismark.align.benchmark.txt"
    shell:
        "bismark --bowtie2 --genome {REFERENCE_DIR} -o {output.alin_dr} -B Testpaired -1 {input.trimmed_fastq_1} -2 {input.trimmed_fastq_2}"

rule bismark_duplicate:
    input:
        bam_file = f"{OUTPUT_DIR}/3.bismark_align/{{sample}}/Testpaired_pe.bam"
    output:
        dup_dr = directory(f"{OUTPUT_DIR}/4.bismark_deduplicate/{{sample}}"),
        duplicate = f"{OUTPUT_DIR}/4.bismark_deduplicate/{{sample}}/Testpaired_pe.deduplicated.bam"
    benchmark:
        "2.benchmarks/{{sample}}.bismark.deduplicate.benchmark.txt"
    shell:
        "deduplicate_bismark -p --bam {input.bam_file} --output_dir {output.dup_dr}"

rule bismark_extract:
    input:
        bam_dedup_file = f"{OUTPUT_DIR}/4.bismark_deduplicate/{{sample}}/Testpaired_pe.deduplicated.bam"
    output:
        extr_dr = directory(f"{OUTPUT_DIR}/5.bismark_methylation/{{sample}}"),
        #txt = f"{OUTPUT_DIR}/5.bismark_methylation/{{sample}}/Testpaired_pe.deduplicated.nonCG_filtered.M-bias.txt"
    benchmark:
        "2.benchmarks/{{sample}}.bismark.extract.benchmark.txt"
    shell:
        "bismark_methylation_extractor --comprehensive --report -p --no_overlap --gzip --mbias_only --parallel 30 --output {output.extr_dr} --genome_folder {REFERENCE_DIR} {input.bam_dedup_file}"




