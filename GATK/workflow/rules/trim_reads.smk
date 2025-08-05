# Trim reads with trimmomatic
rule trimming:
    input:  
        fq1 = lambda w: get_fastqs(w.sample)[0],
        fq2 = lambda w: get_fastqs(w.sample)[1]
    output:
        fq1P = "results/1-trim/{sample}_R1.paired.fastq.gz",
        fq2P = "results/1-trim/{sample}_R2.paired.fastq.gz",
        fq1U = "results/1-trim/{sample}_R1.unpaired.fastq.gz",
        fq2U = "results/1-trim/{sample}_R2.unpaired.fastq.gz",
    params: 
        jar="/usr/local/apps/trimmomatic/0.39/trimmomatic-0.39.jar",
        read ="PE",
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
    log:    logfile = "logs/1-trim/{sample}.log"
    benchmark:
        "benchmarks/1-trim/{sample}.tsv"
    shell:
        """
        module load trimmomatic

        trimmomatic {params.read} \
                -threads {threads} \
                {input.fq1} \
                {input.fq2} \
                {output.fq1P} \
                {output.fq1U} \
                {output.fq2P} \
                {output.fq2U} \
                {params.adapter} \
                {params.leading} \
                {params.headcrop} \
                {params.trailing} \
                {params.slidingwindow} \
                {params.minlen} > {log.logfile} 2>&1
        """
