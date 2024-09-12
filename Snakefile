
# Snakefile to run GATK variant calling pipeline

configfile:"proj_config.yaml"

SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

localrules:

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}.paired_{dir}.fq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/bwa_mem/{sample}.bam", sample = SAMPLES)

rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule trimmomatic:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
	    rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
	    rev = "data/trimming/{sample}.paired_R2.fq.gz",
	    fwd_unpaired = "data/trimming/{sample}.unpaired_R1.fq.gz",
	    rev_unpaired = "data/trimming/{sample}.unpaired_R2.fq.gz"
    conda:
        "envs/trimmomatic.yaml"
    params:
        trimmer = ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"],
	    adapters = config["adapters"]
    shell:
        """
        trimmomatic \
	    PE \
	    -phred33 \
	    {input.fwd} \
	    {input.rev} \
	    {output.fwd} \
	    {output.fwd_unpaired} \
	    {output.rev} \
	    {output.rev_unpaired} \
	    {params.trimmer}
        """

rule bwa_mem:
    input:
        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
        rev = "data/trimming/{sample}.paired_R2.fq.gz"
    output:
        bam_file = "data/bwa_mem/{sample}.bam"
    conda:
        "envs/bwa_mem.yml"
    params:
        bwa_index = config["bwa_index"],
	    rg = config["rg"],
	    threads = 12
    shell:
        """
	    bwa-mem2 mem -v 1 -M -R {params.rg} {params.bwa_index} {input.fwd} {input.rev} | samtools sort -@ {params.threads} -o {output.bam_file} -
        """

#rule haplotype_caller:
#    input:
#        fasta = "/home/exacloud/gscratch/prime-seq/cachedGenomes/128/128_Mmul_10.fasta",
#        bam_file = "data/bwa_mem/{sample}.bam"
#    output:
#        gvcf = "data/haplotype/{sample}.g.vcf.gz"
#    params:
#        tmp_dir = "/path/to/tmp",
#        mem = "56g"
#    conda:
#        "envs/gatk.yml"
#    shell:
#        """
#        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" HaplotypeCaller \
#            -R {input.fasta} \
#            -I {input.bam} \
#            -O {output.gvcf} \
#            -ERC GVCF \
#            -A DepthPerSampleHC
#        """
