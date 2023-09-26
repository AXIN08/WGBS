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

