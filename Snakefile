from snakemake.utils import validate, min_version
import pandas as pd

# ============================================================ 
# Set minimum snakemake version
# ============================================================ 
min_version("9.1.0")

# ============================================================ 
# Load config file
# ============================================================ 
configfile: "07_Config/config.yaml"

# ============================================================ 
# Set variables
# ============================================================ 
ROOTDIR   = os.getcwd()
DATA      = ROOTDIR +  "/00_RawData/"
REFERENCE = ROOTDIR +  "/01_Reference/"
CONTAINER = ROOTDIR +  "/02_Container/"
SCRIPT    = ROOTDIR +  "/03_Script/"
WORKFLOW  = ROOTDIR +  "/04_Workflow/"
OUTPUT    = ROOTDIR +  "/05_Output/"
REPORT    = ROOTDIR +  "/06_Report/"
CONFIG    = ROOTDIR +  "/07_Config/"

# ============================================================ 
# Access variable in the config or other files
# ============================================================ 
# Filenames of Samples
SAMPLES = pd.read_table(CONFIG + config["samples"])
SAMPLES = SAMPLES["ID"].tolist()

# Nomenclature for the paired-end reads
TYPES = config["types"]

# Tools usage variable
TRIMMOMATIC = config["trimmomatic"]["use"]
CUTADAPT    = config["cutadapt"]["use"]
BAM_INDEX   = config["bam_index"]["use"]

# ============================================================ 
# Target output 
# ============================================================ 
# Rules expand variables
FASTQ_STATES = (["multiqc_raw", "multiqc_cleaned"] if (TRIMMOMATIC or CUTADAPT) else ["multiqc_raw"])
ALIGN_EXT    = ([".bam", ".bai"] if BAM_INDEX else [".bam"])
MATRIX       = ["count_matrix", "cpm_matrix", "filtered_cpm_matrix", "filtered_count_matrix", "removed_cpm_matrix", "removed_count_matrix"]

rule all:
	input:
		##### 00_Clean: Read Quality Control and Preprocessing #####
		multiqc_clean  = expand(REPORT + "00_Clean/{fastq_state}.html", fastq_state = FASTQ_STATES),
		##### 01_Alignement: Mapping reads into the reference genome #####
		results_star   = expand(OUTPUT + "01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out{align_ext}", sample = SAMPLES, align_ext = ALIGN_EXT),
		mutiqc_star    = REPORT + "01_Alignment/multiqc_alignment.html",
		##### 02_Counts: Quantification of expression from the mapping results #####
		sample_counts  = expand(OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.tsv", sample = SAMPLES),
		matrix_counts  = expand(OUTPUT + "02_Counts/02_Matrix/{matrix}.tsv", matrix = MATRIX),
		multiqc_counts = REPORT + "02_Counts/multiqc_featurecounts.html",


# ============================================================ 
# Loading rules 
# ============================================================ 
include: WORKFLOW + "00_Clean.smk"
include: WORKFLOW + "01_Alignement.smk"
include: WORKFLOW + "02_Counts.smk"
