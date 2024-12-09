# gatk_vcf
Snakemake pipeline for GATK variant calling.

Steps include:

- **FastQC**
  - sequencing quality check
- **Trimmomatic**
  - trim low quality bases and adapters
- **BWA mem2**
  - align trimmed reads to reference genome
- **GATK HaplotypeCaller**
  - identify variants from alignment files. Outputs a gVCF file per sample.
- **GATK ReblockGVCF**
  - reduces the size of gVCF files
- **GATK GenomicsDBImport**
  - merge multiple gVCF files into a GenomicsDB datastore for joint genotyping. This needs to be re-created every time you add samples that you want to joint genotype.
- **GATK GenotypeGVCFs**
  - performs joint genotyping on the combined samples. Converts the merged gVCF data into actual genotypes for each variant and sample. Limiting 'max-alternate-alleles' becomes important for performance when joint-genotyping a large number of WGS samples, less important for smaller jobs.
- **GATK SampleSpecificGenotypeFiltration**
  - filters out low-quality genotypes for individual samples within a multi-sample VCF.
- **GATK VariantFiltration**
  - applies filters to variant calls based on criteria such as quality scores and read depth.
- **GATK SelectVariants**
  - selects a subset of variants (only SNPs, indels, variants in certain regions, etc) from a VCF file.

**Example execution**   
snakemake --use-conda --jobs 100 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -A {cluster.lab} --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"

**Notes**
sample_name_map.txt file (and intervals file of entire genome) is produced automatically for GenomicsDBImport step.
SampleSpecificGenotypeFiltration step requires "data/sample_type_map.txt" file (to be made manually before execution).
Includes DISCVRSeq-1.3.78.jar in scripts directory, required for SampleSpecificGenotypeFiltration.
