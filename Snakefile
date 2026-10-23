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

# ============================================================ 
# Rules 
# ============================================================ 
# Rules expand variables
FASTQ_STATES = (["multiqc_raw", "multiqc_cleaned"] if (TRIMMOMATIC or CUTADAPT) else ["multiqc_raw"])

rule all:
	input:
		##### 00_Clean: Read Quality Control and Preprocessing #####
		multiqc = expand(REPORT + "00_Clean/{fastq_state}.html", fastq_state = FASTQ_STATES),
		##### 01_Alignement: Mapping reads into the reference genome #####
		bam = expand()

# ============================================================ 
# Rules 
# ============================================================ 
include: WORKFLOW + "00_Clean.smk"
include: WORKFLOW + "01_Alignement.smk"
