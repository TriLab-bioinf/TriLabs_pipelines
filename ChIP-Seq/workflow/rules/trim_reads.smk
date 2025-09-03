# Trim reads with trimmomatic
rule trimming:
    input:  
        fq1=lambda wildcards: [
            my_files for my_files in glob.glob(f"{reads_dir}/{get_fastq(wildcards.sample, wildcards.type,'fastq_1')}")
        ],
        fq2=lambda wildcards: (
            [
                my_files
                for my_files in glob.glob(f"{reads_dir}/{get_fastq(wildcards.sample, wildcards.type,'fastq_2')}")
            ]
            if config["paired"] == True
            else []
        ),
    output:
        fqP1="results/1-trim/{sample}_{type}.P.R1.fastq.gz",
        fqP2="results/1-trim/{sample}_{type}.P.R2.fastq.gz",
        fqUP1="results/1-trim/{sample}_{type}.UP.R1.fastq.gz",
        fqUP2="results/1-trim/{sample}_{type}.UP.R2.fastq.gz",
        failed="results/1-trim/{sample}_{type}.failed.fastq.gz",
        html="results/1-trim/{sample}_{type}.fastp.html",
        json="results/1-trim/{sample}_{type}.fastp.json",
    params:
        par=f"""--trim_poly_x --poly_x_min_len 10 --qualified_quality_phred 20 --length_required 20 --adapter_fasta {adapters} --cut_right --cut_right_mean_quality 20 --cut_right_window_size 5 --n_base_limit 2 """,
        in2=lambda wildcard, input: f"--in2 {str(input.fq2)}" if config["paired"] == True else "",
        out2=lambda wildcard, output: f"--out2 {str(output.fqP2)}" if config["paired"] == True else "",
        unpaired1=lambda wildcard, output: f"--unpaired1  {str(output.fqUP1)}" if config["paired"] == True else "",
        unpaired2=lambda wildcard, output: f"--unpaired2 {str(output.fqUP2)}" if config["paired"] == True else ""
    threads: 8
    resources:
        partition="quick",
        runtime=4 * 60,
        mem_mb=32000,
    log:
        logfile="logs/1-trim/{sample}_{type}.log",
    benchmark:
        "benchmarks/1-trim/{sample}_{type}.tsv"
    shell:
        """
        module load fastp
        
        fastp -i {input.fq1} \
        --out1 {output.fqP1} \
        --failed_out {output.failed} \
        --html {output.html} \
        --json {output.json} {params.in2} {params.out2} {params.unpaired1} {params.unpaired2} {params.par} > {log.logfile} 2>&1

        touch {output.fqP2} {output.fqUP1} {output.fqUP2}
        """