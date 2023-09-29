rule soapnuke:
    input:
        fastq_1 = lambda wildcards: get_path(samples_df, wildcards, 'fq1'),
        fastq_2 = lambda wildcards: get_path(samples_df, wildcards, 'fq2')
    output:
        temp(f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_1_trimmed.fq.gz"),
        temp(f"{OUTPUT_DIR}/1.trimmed/{{sample}}/{{sample}}_2_trimmed.fq.gz")
    benchmark:
        "1.benchmarks/1.trimmed/{sample}.benchmark.txt"
    log:
        "1.logs/1.trimmed/{sample}.log"
    params:
        out_dir=directory(f"{OUTPUT_DIR}/1.trimmed/{{sample}}"),
        clean1=f"{{sample}}_1_trimmed.fq.gz",
        clean2=f"{{sample}}_2_trimmed.fq.gz"
    threads: config["threads"]["soapnuke"]
    shell:
        'SOAPnuke filter -l 30 -q 0.5 -1 {input.fastq_1} -2 {input.fastq_2} -C {params.clean1} -D {params.clean2} -o {params.out_dir} -T {threads}'
