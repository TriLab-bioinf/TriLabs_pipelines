rule VariantRecalibrator:
    input:
        vcf = "results/6-merged/all.wgs.combined.vcf.gz",
        ref = f"{genome}"
    output:
        SNP_recal = "results/7-filtered_vcf/all.wgs.SNP.recal",
        SNP_tranches = "results/7-filtered_vcf/all.wgs.SNP.tranches",
        SNP_Rscript = "results/7-filtered_vcf/all.wgs.SNP.plots.R",
        indel_recal = "results/7-filtered_vcf/all.wgs.indel.recal",
        indel_tranches = "results/7-filtered_vcf/all.wgs.indel.tranches",
        indel_Rscript = "results/7-filtered_vcf/all.wgs.indel.plots.R"
    threads: 16
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    log: 
        logfile = "logs/7-filtered_vcf/all.wgs.VariantRecalibrator.log"
    shell:
        """
        module load GATK/4.6.0.0

        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=8" VariantRecalibrator \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 \
          -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
          -tranche 95.0 -tranche 94.0 \
          -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
          -R {input.ref} \
          -V {input.vcf} \
          --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
          /fdb/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz  \
          --resource:omni,known=false,training=true,truth=false,prior=12.0 \
          /fdb/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
          --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
          /fdb/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
          -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
          -mode SNP -O {output.SNP_recal} --tranches-file {output.SNP_tranches} \
          --rscript-file {output.SNP_Rscript} > {log.logfile} 2>&1

        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=8" VariantRecalibrator \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 \
          -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
          -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
          -tranche 92.0 -tranche 91.0 -tranche 90.0 \
          -R {input.ref} \
          -V {input.vcf} \
          --resource:mills,known=false,training=true,truth=true,prior=12.0 \
          /fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
          --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
          /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz \
          -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
          -mode INDEL -O {output.indel_recal} --tranches-file {output.indel_tranches} \
          --rscript-file {output.indel_Rscript} >> {log.logfile} 2>&1

        """

rule apply_vqsr:
    input:
        SNP_recal = "results/7-filtered_vcf/all.wgs.SNP.recal",
        SNP_tranches = "results/7-filtered_vcf/all.wgs.SNP.tranches",   
        indel_recal = "results/7-filtered_vcf/all.wgs.indel.recal",
        indel_tranches = "results/7-filtered_vcf/all.wgs.indel.tranches",
        vcf = "results/6-merged/all.wgs.combined.vcf.gz"
    output:
        SNP_vcf = "results/7-filtered_vcf/all.wgs.SNP.recalibrated_99.9.vcf.gz",
        indel_vcf = "results/7-filtered_vcf/all.wgs.indel.SNP.recalibrated_99.9.vcf.gz"
    threads: 16
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    log: 
        logfile = "logs/7-filtered_vcf/all.wgs.apply_vqsr.log"
    shell:
        """
        module load GATK/4.6.0.0

        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=8" ApplyVQSR \
          -V {input.vcf} \
          --recal-file {input.SNP_recal} \
          -mode SNP \
          --tranches-file {input.SNP_tranches} \
          --truth-sensitivity-filter-level 99.9 \
          --create-output-variant-index true \
          -O {output.SNP_vcf} > {log.logfile} 2>&1

        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
          -V {output.SNP_vcf} \
          -mode INDEL \
          --recal-file {input.indel_recal} \
          --tranches-file {input.indel_tranches} \
          --truth-sensitivity-filter-level 99.9 \
          --create-output-variant-index true \
          -O {output.indel_vcf} >> {log.logfile} 2>&1

        """