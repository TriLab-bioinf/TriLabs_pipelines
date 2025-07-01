# Trim reads with trimmomatic
rule trimming:
    input:  
        fq1 = lambda wildcards: glob.glob(f"data/reads/{get_fastq(wildcards.sample, wildcards.type,'fastq_1')}")
    output:
        fq1P = "results/1-trim/{sample}_{type}.trimmed.fastq.gz" 
    params: 
        jar="/usr/local/apps/trimmomatic/0.39/trimmomatic-0.39.jar",
        read ="SE",
        adapter=f"ILLUMINACLIP:{adapters}:2:30:12",
        leading="LEADING:0",
        trailing="TRAILING:0",
        slidingwindow="SLIDINGWINDOW:4:20",
        minlen="MINLEN:20",
        headcrop="HEADCROP:0"
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 32000
    threads: 16
    log:    logfile = "results/1-trim/{sample}_{type}.log"
    benchmark:
        "benchmarks/1-trim/{sample}_{type}.tsv"
    shell:
        """
        module load trimmomatic

        trimmomatic {params.read} \
                -threads {threads} \
                {input.fq1} \
                {output.fq1P} \
                {params.adapter} \
                {params.leading} \
                {params.headcrop} \
                {params.trailing} \
                {params.slidingwindow} \
                {params.minlen} > {log.logfile} 2>&1
        """
