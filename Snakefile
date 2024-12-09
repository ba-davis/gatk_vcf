
# Snakefile to run GATK variant calling pipeline

configfile:"proj_config.yaml"

SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

localrules: samtools_flagstat, flagstat_summary, cleanup_unsorted_bam, sample_name_map, make_intervals

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}.paired_{dir}.fq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/bwa_mem/{sample}.bam", sample = SAMPLES),
        expand("data/bwa_mem/{sample}.sorted.bam", sample = SAMPLES),
        expand("data/bwa_mem/{sample}.flagstat.txt", sample = SAMPLES),
        "data/bwa_mem/flagstat_summary.txt",
        expand("data/haplotype/{sample}.g.vcf.gz", sample = SAMPLES),
        expand("data/reblock/{sample}.rb.g.vcf.gz", sample = SAMPLES),
        "data/genomicsdb/sample_name_map.txt",
        "data/genomicsdb/intervals.list",
        "data/genomicsdb/genomics_db",
        "data/genotype/genotyped.vcf.gz",
        "data/genotype/genotyped.filtered.vcf.gz",
        "data/variant_filtration/genotyped.filtered.vfilt.vcf.gz",
        "data/select_variants/genotyped.filtered.vfilt.selectVariants.vcf.gz"

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
        "envs/bwa_mem.yaml"
    params:
        bwa_index = config["bwa_index"],
	    rg = config["rg"],
	    threads = 12
    shell:
        """
	    bwa-mem2 mem -v 1 -M -R {params.rg} {params.bwa_index} {input.fwd} {input.rev} | samtools sort -@ {params.threads} -o {output.bam_file} -
        """

rule sort_and_index_bam:
    input:
        bam_file = "data/bwa_mem/{sample}.bam"
    output:
        sorted_bam_file = "data/bwa_mem/{sample}.sorted.bam",
        sorted_bam_index = "data/bwa_mem/{sample}.sorted.bam.bai"
    conda:
        "envs/bwa_mem.yaml"
    params:
        threads = 12
    shell:
        """
        samtools sort -@ {params.threads} -o {output.sorted_bam_file} {input.bam_file} &&
        samtools index {output.sorted_bam_file}
        """

rule cleanup_unsorted_bam:
    input:
        unsorted_bam = "data/bwa_mem/{sample}.bam"
    output:
        temp("data/bwa_mem/{sample}.unsorted.bam.deleted")  # Use a temporary output to ensure this rule is run after sorting
    shell:
        """
        rm {input.unsorted_bam}
        """

rule samtools_flagstat:
    input:
        bam_file = "data/bwa_mem/{sample}.sorted.bam"
    output:
        flagstat = "data/bwa_mem/{sample}.flagstat.txt"
    conda:
        "envs/bwa_mem.yaml"
    params:
        threads = 4
    shell:
        """
        samtools flagstat -@ {params.threads} {input.bam_file} > {output.flagstat}
        """

rule flagstat_summary:
    input:
        expand("data/bwa_mem/{sample}.flagstat.txt", sample=SAMPLES)
    output:
        summary = "data/bwa_mem/flagstat_summary.txt"
    conda:
        "envs/pandas.yaml"
    shell:
        """
        python scripts/parse_flagstat.py {output.summary} {input}
        """

rule haplotype_caller:
    input:
        fasta = config["fasta"],
        bam_file = "data/bwa_mem/{sample}.sorted.bam"
    output:
        gvcf = "data/haplotype/{sample}.g.vcf.gz"
    params:
        tmp_dir = "data/",
        mem = "56g"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" HaplotypeCaller \
            -R {input.fasta} \
            -I {input.bam_file} \
            -O {output.gvcf} \
            -ERC GVCF \
            -A DepthPerSampleHC
        """

rule reblock_gvcf:
    input:
        fasta = config["fasta"],
        gvcf = "data/haplotype/{sample}.g.vcf.gz"
    output:
        rb_gvcf = "data/reblock/{sample}.rb.g.vcf.gz"
    params:
        tmp_dir = "data/",
        mem = "56g"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" ReblockGVCF \
            -R {input.fasta} \
            -V {input.gvcf} \
            -O {output.rb_gvcf}
        """

rule sample_name_map:
    input:
        rb_gvcfs = expand("data/reblock/{sample}.rb.g.vcf.gz", sample=SAMPLES)
    output:
        sample_map = "data/genomicsdb/sample_name_map.txt"
    run:
        with open(output.sample_map, "w") as out_file:
            for gvcf in input.rb_gvcfs:
                sample_name = gvcf.split("/")[-1].replace(".rb.g.vcf.gz", "")
                out_file.write(f"{sample_name}\t{gvcf}\n")

rule make_intervals:
    input:
        sample_map = "data/genomicsdb/sample_name_map.txt",
        fasta_dict = config["fasta_dict"]
    output:
        intervals = "data/genomicsdb/intervals.list"
    shell:
        """
        (grep "^@SQ" {input.fasta_dict} | cut -f2 | sed 's/SN://') > {output.intervals}
        """

rule genomics_db_import:
    input:
        sample_map = "data/genomicsdb/sample_name_map.txt",
        fasta = config["fasta"],
        intervals = "data/genomicsdb/intervals.list"
    output:
        db_path = directory("data/genomicsdb/genomics_db")
    params:
        tmp_dir = "data/",
        mem = "56g",
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" GenomicsDBImport \
            --genomicsdb-workspace-path {output.db_path} \
            --sample-name-map {input.sample_map} \
            --bypass-feature-reader \
            --genomicsdb-shared-posixfs-optimizations \
            --reader-threads 4 \
            --batch-size 50 \
            -R {input.fasta} \
            -L {input.intervals}
        """

rule genotype_gvcfs:
    input:
        db_path = "data/genomicsdb/genomics_db",
        fasta = config["fasta"],
        intervals = "data/genomicsdb/intervals.list"
    output:
        vcf = "data/genotype/genotyped.vcf.gz"
    params:
        tmp_dir = "data/",
        mem = "72g"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" GenotypeGVCFs \
            -R {input.fasta} \
            --variant gendb://{input.db_path} \
            -O {output.vcf} \
            -L {input.intervals} \
            --annotate-with-num-discovered-alleles \
            -stand-call-conf 20 \
            --max-alternate-alleles 5 \
            --genomicsdb-max-alternate-alleles 8 \
            #--only-output-calls-starting-in-intervals \
            --genomicsdb-shared-posixfs-optimizations
        """

rule sample_specific_genotype_filtration:
    input:
        fasta = config["fasta"],
        vcf = "data/genotype/genotyped.vcf.gz",
        sample_map = "sample_type_map.txt"
    output:
        filtered_vcf = "data/genotype/genotyped.filtered.vcf.gz"
    params:
        tmp_dir = "data/",
        mem = "120g",
        jar = "scripts/DISCVRSeq-1.3.78.jar"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        java -Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m \
        -jar {params.jar} SampleSpecificGenotypeFiltration \
        -R {input.fasta} \
        -V {input.vcf} \
        -O {output.filtered_vcf} \
        --set-filtered-genotype-to-no-call \
        --genotype-filter-name DP-LT10 \
        --genotype-filter-expression "WGS:DP<10" \
        --genotype-filter-name DP-GT100 \
        --genotype-filter-expression "WGS:DP>100" \
        --genotype-filter-name GQ-LT20a \
        --genotype-filter-expression "WGS:g.hasGQ() && !g.hasExtendedAttribute('RGQ') && GQ<20" \
        --genotype-filter-name DP-LT10 \
        --genotype-filter-expression "WXS:DP<10" \
        --genotype-filter-name GQ-LT20a \
        --genotype-filter-expression "WXS:g.hasGQ() && !g.hasExtendedAttribute('RGQ') && GQ<20" \
        --sample-map {input.sample_map}
        """

# annotation rule? VariantAnnotator

rule variant_filtration:
    input:
        fasta = config["fasta"],
        vcf = "data/genotype/genotyped.filtered.vcf.gz",
        mask = config["repeat_mask_bed"]
    output:
        vcf = "data/variant_filtration/genotyped.filtered.vfilt.vcf.gz"
    params:
        tmp_dir = "data/",
        mem = "120g",
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" VariantFiltration \
            -R {input.fasta} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-name QualityFilter \
            --filter "vc.hasAttribute('QD') && QD < 2.0" \
            --filter-name FisherStrand \
            --filter "vc.hasAttribute('FS') && ((!vc.isIndel() && FS > 60.0) || (vc.isIndel() && FS > 200.0))" \
            --filter-name StrandOddsRatio \
            --filter "vc.hasAttribute('SOR') && ((!vc.isIndel() && SOR > 3.0) || (vc.isIndel() && SOR > 10.0))" \
            --filter-name MappingQuality \
            --filter "vc.hasAttribute('MQ') && (!vc.isIndel() && MQ < 40.0)" \
            --filter-name MQRankSum \
            --filter "vc.hasAttribute('MQRankSum') && (!vc.isIndel() && MQRankSum < -12.5)" \
            --filter-name ReadPosRankSum \
            --filter "vc.hasAttribute('ReadPosRankSum') && ((!vc.isIndel() && ReadPosRankSum < -8.0) || (vc.isIndel() && ReadPosRankSum < -20.0))" \
            --filter-name NoneCalled \
            --filter "vc.getCalledChrCount() == 0" \
            --mask-name RepeatMask \
            --mask {input.mask}
        """

rule select_variants:
    input:
        fasta = config["fasta"],
        vcf = "data/variant_filtration/genotyped.filtered.vfilt.vcf.gz"
    output:
        vcf = "data/select_variants/genotyped.filtered.vfilt.selectVariants.vcf.gz"
    params:
        tmp_dir = "data/",
        mem = "120g",
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk --java-options "-Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -Xms{params.mem} -Xss2m" SelectVariants \
            -R {input.fasta} \
            -V {input.vcf} \
            -O {output.vcf} \
            --select "vc.getCalledChrCount() > 0" \
            --allow-nonoverlapping-command-line-samples \
            --set-filtered-gt-to-nocall
        """
