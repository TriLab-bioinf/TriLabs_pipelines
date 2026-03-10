import os

# Run Bismark on trimmed reads
rule bismark:
    input:
        fq1 = lambda wildcards: "results/1-trim/{sample}.P.R1.fastq.gz",
        fq2 = lambda wildcards: "results/1-trim/{sample}.P.R2.fastq.gz",
        genome = f"{config['bismark']['path']}"
    output: 
        bam = "results/2-bismark/{sample}/bismark.bam",
        report = "results/2-bismark/{sample}/bismark_report.txt",
        tmp = directory("results/2-bismark/{sample}.tmp/"),
        dir = directory("results/2-bismark/{sample}/"),
    threads: 38
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    benchmark: "benchmarks/2-bismark/{sample}.bismark.tsv"
    params: 
        other_bismark_params = f"{config['bismark']['other_bismark_params']}"
    log: 
        logfile = "logs/2-bismark/{sample}.bismark.log"
    shell:
        """
        module load bismark/0.25.0 bowtie/2-2.5.3
            
        bismark -p 3 -n 1 --multi 4 \
            --temp_dir {output.tmp} \
            --output_dir {output.dir} \
            --rg_id {wildcards.sample} \
            --rg_sample {wildcards.sample} \
            {params.other_bismark_params} \
            {input.genome} -1 {input.fq1} -2 {input.fq2} > {log.logfile} 2>&1
        
        mv {output.dir}/{wildcards.sample}.P.R1_bismark_bt2_pe.bam {output.bam}
        mv {output.dir}/{wildcards.sample}.P.R1_bismark_bt2_PE_report.txt {output.report}
        """