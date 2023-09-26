rule bismark_alignment:
    input:
        trimmed_fastq_1 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_1_trimmed.fq.gz",
        trimmed_fastq_2 = f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_2_trimmed.fq.gz"
    output:
        alin_dr = directory(f"{OUTPUT_DIR}/2.bismark_align/{{sample}}"),
        bam = f"{OUTPUT_DIR}/2.bismark_align/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.bam"
        #html = f"{OUTPUT_DIR}/3.bismark_align/{{sample}}/Testpaired_PE_report.html"
    benchmark:
        "{OUTPUT_DIR}/benchmarks/2.bismark_align/{sample}.benchmark.txt"
    log:
        "{OUTPUT_DIR}/logs/2.bismark_align/{sample}.log"
    resources: mem_mb = 122880 
    threads: config["threads"]["bis_align"]
    shell:
        'bismark --bowtie2 --multicore {threads} --genome {REFERENCE_DIR} -o {output.alin_dr} -1 {input.trimmed_fastq_1} -2 {input.trimmed_fastq_2} 2> {log}'

