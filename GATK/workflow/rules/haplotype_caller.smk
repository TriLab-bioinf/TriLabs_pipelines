
rule haplotype_caller:
  input:
    bam = "results/4-bqsr/{sample}.recal.bam",
    ref = f"{genome}"
  output:
    vcf = "results/5-vcf/{sample}.raw.g.vcf.gz"
  threads: 16
  resources:
    partition = "norm",
    runtime = 14 * 60,
    mem_mb = 1024 * 64,
    disk_mb= 1024 * 20
  log: 
    logfile = "logs/5-vcf/{sample}.haplotype_caller.log"
  shell:
    """
    module load GATK/4.6.0.0

    mkdir -p results/5-vcf

    gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms16G -Xmx16G -XX:ParallelGCThreads=8" HaplotypeCaller \
      -R {input.ref} \
      -I {input.bam} \
      -O {output.vcf} \
      -ERC GVCF > {log.logfile} 2>&1
    """