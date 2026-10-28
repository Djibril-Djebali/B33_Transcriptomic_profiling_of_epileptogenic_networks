# Transcriptomic profiling of epileptogenic networks

## Authors

Djibril Djebali (@Djibril-Djebali), https://centuri-livingsystems.org/d-djebali

## Workflow Overview

This workflow performs a bulkRNA-seq analysis to study gene expression and gene variant from the sequencing output.

The different analysis steps:
 - [01_Clean](04_Workflow/01_Clean.smk) the quality of the raw reads are assessed using [FastQC v0.13.0](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC v1.35](https://github.com/multiqc/multiqc), and then they are trimmed with [Trimmomatic v0.36](https://github.com/usadellab/trimmomatic) to remove bad quality sequences and adapters are cut with [Cutadapt v5.2](https://cutadapt.readthedocs.io/en/stable/).
 - [02_Alignment](04_Workflow/02_Alignment.smk) the genome index is built and used to map the raw reads with [STAR v2.7.11b](https://github.com/alexdobin/STAR) and [Samtools v1.21](https://github.com/samtools/samtools) to convert resulting sam files in bam files.
 - [03_Counts](04_Workflow/03_Counts.smk) Quantification of gene-level read counts from aligned reads and generating count matrices with FeatureCounts from the package [Subreads v2.0.8](https://subread.sourceforge.net/featureCounts.html) and filtering of lowly expressed genes based on CPM thresholds to improve statistical robustness.

## Usage

### Step 1: Install workflow

You can either download the workflow and extract the latest release, or clone it. If you intend to modify and further extend this workflow or want to work under version control, you can fork this repository.

```bash
git clone https://github.com/Djibril-Djebali/B33_Transcriptomic_profiling_of_epileptogenic_networks.git
```
We would be pleased if you use this workflow and participate in its improvement. If you use it in a paper, don't forget to give credits to the author by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files and repositories:

- 00_RawData needs pair-end fastq file of each run to analyse.
- [config.yaml](07_Config/config.tsv) indicating the workflow parameters to use.
- [samples.tsv](07_Config/samples.tsv) to indicate the fastq IDs and their sample names.
- [metadata.tsv](07_Config/metadata.tsv) to indicate the fastq IDs and the metadata associated.

### Step 3: Execute workflow

- You need [Apptainer v1.4.4](https://github.com/apptainer/apptainer/blob/main/INSTALL.md) installed on your computer or cluster.
- Retrieve annovar.latest.tar.gz file in the following site: [Annovar download](https://annovar.openbioinformatics.org/en/latest/user-guide/download/), untar it un put the annovar folder obtained in 02_Container/.
- Pull the Snakemake v9.27.0 from the docker repositories and convert it into an Apptainer image with the following command from the workflow folder.

```bash
mkdir -p ./.apptainer
apptainer pull \
	--disable-cache \ 
	./.apptainer/snakemake_v9.27.0.sif \
	docker://snakemake/snakemake:v9.27.0
```

- Then execute the workflow with the following command from the workflow folder.

```bash
apptainer exec ./.apptainer/snakemake_v9.27.0.sif \
		snakemake \
			-s Snakefile \
			--use-conda \
			--sdm apptainer \
			--cores all 
```

Otherwise you can use [GNU make v4.3](https://github.com/mirror/make) and pull the Apptainer image of snakemake and launch the workflow by using `make build` and `make run` command from the workflow folder.


### Step 4: Investigate results 

After successful execution, you can access the results in: