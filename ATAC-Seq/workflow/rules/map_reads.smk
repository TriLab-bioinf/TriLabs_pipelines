# Map reads to reference with bowtie2

# Set bowtiedb path and check if bowtiedb already exists, otherwise create one
if (bowtiedb_path == "None") and (not os.path.exists("data/bowtiedb/genome")):
    db_path: str = "data/bowtiedb/"

    # Make bowtiedb
    rule make_bowtie_index:
        input: gen = f"{genome}"
        output: db = "data/bowtiedb/genome",
        resources: 
            mem_mb = 32000,
            partition = "quick",
            runtime = 14 * 60
        threads: 8
        params: 
        shell:
            """
            module load bowtie

            bowtie2-build --threads {threads} {input.gen} {output.db} 
                
            touch {output.db}
            """

elif (bowtiedb_path == "None") and (os.path.exists("data/bowtiedb/genome")):
    # bowtiedb was already created => skip make_bowtie_index
    db_path: str = "data/bowtiedb/"
else:
    # Use an existent bowtiedb located somewhere else
    db_path: str = bowtiedb_path

# Mapping reads to reference
rule map_reads:
    input:
        fq1 = "results/1-trim/{sample}_R1.paired.fastq.gz",
        fq2 = "results/1-trim/{sample}_R2.paired.fastq.gz",
        db = lambda wildcards: rules.make_bowtie_index.output.db if ((bowtiedb_path == "None") and (not os.path.exists("data/bowtiedb/genome"))) else []
    output: 
        bam = "results/2-map_reads/{sample}.bam",
        bai = "results/2-map_reads/{sample}.bam.bai"
    threads: 16
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    benchmark:
        "benchmarks/2-map_reads/{sample}.bowtie2.tsv"
    params: 
        prefix = "{sample}",
        bowtiedb_dir = f"{db_path}",
        db_base = "data/bowtiedb/genome"
    log: "logs/2-map_reads/{sample}.bowtie2.log"        
    shell:
        """
        module load bowtie samtools

        bowtie2 -x {params.db_base} \
        --threads {threads} \
        --phred33 \
        --rg-id {params.prefix} --rg SM:{params.prefix} --rg LB:library --rg PL:ILLUMINA \
        -1 {input.fq1}  -2 {input.fq2} \
        --sensitive | samtools sort -@ 8 -O BAM -o {output.bam} -
        
        samtools index -@ 8 {output.bam}
        """