rule merge_vcfs:
    input:
        vcfs = expand("results/8-all-calls-by-chr/all_samples.{chrom}.vcf.gz", chrom = get_chr_ids(whole_chr_intervals)),
    output:
        vcf_files = "results/9-all-calls/all_vcf_files.list",
        merged_vcf = "results/9-all-calls/all.vcf.gz"
    log:
        "logs/9-all-calls/merge_vcfs.log"
    params:
        java_opts = "",  # optional
    resources:
        mem_mb = 8 * 1024,
        runtime = 4 * 60,
        partition = "normal"
    benchmark:
        "benchmarks/all-calls/merge-vcfs.tsv"
    shell:
        """
        echo {input.vcfs} | tr ' ' '\n' > {output.vcf_files}
        
        module load GATK/4.6.0.0
        
        gatk --java-options "-Xmx4g" MergeVcfs \
            -I {output.vcf_files} \
            -O {output.merged_vcf} > {log} 2>&1
        """


rule select_variants:
    input:
        vcf = "results/9-all-calls/all.vcf.gz"
    output:
        snps = "results/10-subset-snps/all_SNP.vcf.gz",
        indels = "results/11-subset-indels/all_INDELS.vcf.gz"
    log:
        snps = "logs/subset/snps.log",
        indels = "logs/subset/indels.log"
    params:
        extra = ""
    resources:
        mem_mb = 1024 * 32,
        runtime = 2 * 60,
        partition = "quick"
    benchmark:
        "benchmarks/subset-variants/subset.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk SelectVariants \
            -V {input.vcf} \
            -select-type SNP \
            -O {output.snps}  > {log.snps} 2>&1

        gatk SelectVariants \
            -V {input.vcf} \
            -select-type INDEL \
            -O {output.indels}  > {log.indels} 2>&1
        """
