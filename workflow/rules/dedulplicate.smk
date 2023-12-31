rule bismark_duplicate:
    input:
        bam_file = f"{OUTPUT_DIR}/2.bismark_align/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.bam"
    output:
        dup_dr = directory(f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}"),
        duplicate = f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.deduplicated.bam"
    benchmark:
        "1.benchmarks/3.bismark_deduplicate/{sample}.benchmark.txt"
    log:
        "1.logs/3.bismark_deduplicate/{sample}.log"
    threads: config["threads"]["bis_dedul"]
    priority: 50
    shell:
        'deduplicate_bismark -p --parallel {threads} --bam {input.bam_file} --output_dir {output.dup_dr} 2> {log}'

