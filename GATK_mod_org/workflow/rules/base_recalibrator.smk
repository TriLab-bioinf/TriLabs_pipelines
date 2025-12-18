rule gatk_base_recalibrator:
    input:
        bam = "results/3-dedup/{sample}.dedup.bam",
        ref = f"{genome}",
        refdict = f"{genome_dict}",
        known = f"{dbsnp}"  # optional known sites - single or a list
    output:
        recal_table = "results/4-recal/{sample}.grp",
    log:
        "logs/4-recal/{sample}.log",
    params:
        extra = "",  # optional
        java_opts = "",  # optional
    resources:
        partition = "normal",
        runtime = 24 * 60,
        mem_mb = 1024 * 16,
    shell:
        """
        module load GATK/4.6.0.0

        gatk BaseRecalibrator \
			-R {input.ref} \
			-I {input.bam} \
			--known-sites {input.known} \
			-O {output.recal_table} > {log} 2>&1
        
        """

# Apply recalibration with GATK ApplyBQSR
rule gatk_apply_bqsr:
    input:
        bam = "results/3-dedup/{sample}.dedup.bam",
        ref = f"{genome}",
        refdict = genome_dict,
        recal_table = "results/4-recal/{sample}.grp",
    output:
        bam = "results/5-recal/{sample}.recal.bam",
    log:
        "logs/5-recal/{sample}.log",
    params:
        extra = "",  # optional
        java_opts = "",  # optional
        embed_ref = True,  # embed the reference in cram output
    threads: 8
    resources:
        partition = "normal",
        runtime = 24 * 60,
        mem_mb = 1024 * 16
    shell:
        """
        module load GATK/4.6.0.0

        gatk ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.bam} > {log} 2>&1
        """
