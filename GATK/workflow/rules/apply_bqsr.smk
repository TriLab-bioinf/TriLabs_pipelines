rule base_recalibration:
  input:
    bam = "results/3-dedup/{sample}.dedup.bam",
    ref = f"{genome}",
    known_sites = config["known_sites"]
  output:
    table = "results/4-bqsr/{sample}.recal_data.table"
  threads: 8
  resources:
    partition = "norm",
    runtime = 14 * 60,
    mem_mb = 1024 * 64,
    disk_mb= 1024 * 20
  params:
    known_sites = lambda wildcards, input: ' '.join(f"--known-sites {site}" for site in input.known_sites)
  shell:
    """
    module load GATK/4.6.0.0

    gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=4" BaseRecalibrator \
      -I {input.bam} \
      -R {input.ref} \
      {params.known_sites} \
      -O {output.table}
    """

rule apply_bqsr:
  input:
    bam = "results/3-dedup/{sample}.dedup.bam",
    table = "results/4-bqsr/{sample}.recal_data.table",
    ref = f"{genome}"
  output:
    "results/4-bqsr/{sample}.recal.bam"
  threads: 8
  resources:
    partition = "norm",
    runtime = 14 * 60,
    mem_mb = 1024 * 64,
    disk_mb= 1024 * 20
  shell:
    """
    module load GATK/4.6.0.0

    gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=4" ApplyBQSR \
      -R {input.ref} \
      -I {input.bam} \
      --bqsr-recal-file {input.table} \
      -O {output}
    """
