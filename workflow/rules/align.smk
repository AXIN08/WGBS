rule bismark_alignment:
    input:
        trimmed_fastq_1 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_1_trimmed.fq.gz",
        trimmed_fastq_2 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_2_trimmed.fq.gz"
    output:
        alin_dir=directory(f"{OUTPUT_DIR}/2.bismark_align/{{sample}}"),
        bam_file = f"{OUTPUT_DIR}/2.bismark_align/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.bam"
    benchmark:
        "1.benchmarks/2.bismark_align/{sample}.benchmark.txt"
    log:
        "1.logs/2.bismark_align/{sample}.log"
    threads: config["threads"]["bis_align"]
    shell:
        'bismark --bowtie2 --multicore {threads} --genome {REFERENCE_DIR} -o {output.alin_dir} -1 {input.trimmed_fastq_1} -2 {input.trimmed_fastq_2} 2> {log}'
