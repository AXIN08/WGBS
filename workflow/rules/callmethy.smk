rule bismark_extract:
    input:
        bam_dedup_file = f"{OUTPUT_DIR}/3.bismark_deduplicate/{{sample}}/{{sample}}_1_trimmed_bismark_bt2_pe.deduplicated.bam"
    output:
        extr_dr = directory(f"{OUTPUT_DIR}/4.bismark_methylation/{{sample}}")
        #txt = f"{OUTPUT_DIR}/5.bismark_methylation/{{sample}}/Testpaired_pe.deduplicated.nonCG_filtered.M-bias.txt"
    benchmark:
        "1.benchmarks/4.bismark_methylation/{sample}.benchmark.txt"
    log:
        "1.logs/4.bismark_methylation/{sample}.log"
    threads: config["threads"]["bis_call"]
    priority: 100
    shell:
        'bismark_methylation_extractor --cytosine_report -p --no_overlap --gzip --multicore {threads} --output {output.extr_dr} --genome_folder {REFERENCE_DIR} {input.bam_dedup_file} 2> {log} && cp {input.bam_dedup_file} /mnt/magic_nas/liuxin/WGBS_bamdata'

