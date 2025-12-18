# Map reads to reference with bwa2

# Set bwadb path and check if bwadb already exists, otherwise create one
if (bwadb_path == "None") and (not os.path.exists("data/bwadb/genome.fa")):
    db_path: str = "data/bwadb/"

    # Make bwadb
    rule make_bwa_index:
        input: gen = f"{genome}"
        output: db = "data/bwadb/genome.fa",
        resources: 
            partition = "norm",
            runtime = 14 * 60,
            mem_mb = 1024 * 100,
            disk_mb= 1024 * 50
        threads: 8
        params: 
        shell:
            """
            module load bwa-mem2/2.2.1

            bwa-mem2 index {input.gen} -p {output.db} 
                
            touch {output.db}
            """

elif (bwadb_path == "None") and (os.path.exists("data/bwadb/genome.fa")):
    # bwadb was already created => skip make_bwa_index
    db_path: str = "data/bwadb/"
else:
    # Use an existent bwadb located somewhere else
    db_path: str = bwadb_path

# Mapping reads to reference

# Create dictionaries from metadata for @RG line
# metadata is imported from ./config/samplesheet.csv file
fq_1: dict = {s:fq1 for s, fq1 in zip(metadata['sample_ID'], metadata['fastq_1'])}
fq_2: dict = {s:fq2 for s, fq2 in zip(metadata['sample_ID'], metadata['fastq_2'])}
samples: list = list(metadata['sample_ID'])
ID: dict = {s:f"{fq1}{lan}" for s, fq1, lan in zip(metadata['sample_ID'], metadata['flowcell_ID'], metadata['lane'])}
SM = samples
PL: dict = {s:tech for s, tech in zip(metadata['sample_ID'], metadata['technology'])}
LB: dict = {s:lib for s, lib in zip(metadata['sample_ID'], metadata['library'])}

rule map_reads:
    input:
        fq1 = lambda wildcards: [f"results/1-trim/{wildcards.sample}.P.R1.fastq.gz" if (config["trimming"] == True) else f"data/reads/{fq_1[wildcards.sample]}"],
        fq2 = lambda wildcards: (
                [
                    f"results/1-trim/{wildcards.sample}.P.R2.fastq.gz" if (config["trimming"] == True) else f"data/reads/{fq_2[wildcards.sample]}"
                ]
                if config["paired"] == True
                else []
            ),
        db = lambda wildcards: rules.make_bwa_index.output.db if ((bwadb_path == "None") and (not os.path.exists("data/bwadb/genome.fa"))) else []
    output: 
        bam = "results/2-map_reads/{sample}.bam",
        bai = "results/2-map_reads/{sample}.bam.csi"
    threads: 16
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    benchmark:
        "benchmarks/2-map_reads/{sample}.bwa.tsv"
    params: 
        prefix = "{sample}",
        bwadb_dir = f"{db_path}",
        db_base = "genome.fa",
        rg_line = lambda wildcards: f"\"ID:{ID[wildcards.sample]}\t" + 
                                    f"SM:{wildcards.sample}\t" + 
                                    f"PL:{PL[wildcards.sample]}\t" + 
                                    f"LB:{LB[wildcards.sample]}\t" + 
                                    f"PU:{ID[wildcards.sample]}\""
    log: 
        logfile = "logs/2-map_reads/{sample}.bwa.log"        
    shell:
        """
        module load bwa/0.7.17 samtools/1.21

        bwa mem {params.bwadb_dir}/{params.db_base} \
            {input.fq1} {input.fq2} -t {threads} | \
        samtools view -hb - | \
        samtools addreplacerg -r {params.rg_line} -O BAM - | \
        samtools sort -@ {threads} \
            -O BAM -T {wildcards.sample}.tmp \
            --write-index \
            -o {output.bam} - > {log.logfile} 2>&1
        """