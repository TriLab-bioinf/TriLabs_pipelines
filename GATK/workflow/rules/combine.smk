rule combine_calls:
  input:
      ref=f"{genome}",
      gvcfs=expand("results/5-vcf/{sample}.raw.g.vcf.gz", sample=samples)
  output:
      gvcf="results/6-merged/all.wgs.g.vcf.gz"
  threads: 16
  resources:
      partition = "norm",
      runtime = 14 * 60,
      mem_mb = 1024 * 64,
      disk_mb= 1024 * 20
  log:
      "logs/6-merged/combinegvcfs.log"
  params:
      gvcf_args = lambda wildcards, input: " ".join("-V {}".format(f) for f in input.gvcfs)
  shell:
    """
    module load GATK/4.6.0.0

    gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms8G -Xmx8G -XX:ParallelGCThreads=8" CombineGVCFs \
      {params.gvcf_args} \
      -R {input.ref} \
      -O {output.gvcf}

    """

rule genotype_variants:
  input:
      ref=f"{genome}",
      gvcf="results/6-merged/all.wgs.g.vcf.gz",
  output:
      vcf="results/6-merged/all.wgs.combined.vcf.gz"
  threads: 16
  resources:
      partition = "norm",
      runtime = 14 * 60,
      mem_mb = 1024 * 64,
      disk_mb= 1024 * 20
  log:
      "logs/6-merged/genotypegvcfs.log"
  shell:
    """
    module load GATK/4.6.0.0

    gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms8G -Xmx8G -XX:ParallelGCThreads=8" GenotypeGVCFs \
      -V {input.gvcf} \
      -R {input.ref} \
      -O {output.vcf}

    """

