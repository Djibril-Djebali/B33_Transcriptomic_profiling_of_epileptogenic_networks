# ============================================================
#  00_Clean : Read Quality Control and Preprocessing
# ============================================================
#
#  This workflow processes raw paired-end FASTQ reads through
#  a configurable quality control and cleaning pipeline.
#
#  ──────────────────────── PIPELINE ─────────────────────────
#
#  [1] FastQC (Raw)
#      └─ Quality assessment of raw reads.
#         HTML report + ZIP per sample per read type.
#
#  [2] MultiQC (Raw)
#      └─ Aggregates all raw FastQC reports into a single
#         interactive HTML summary report.
#
#  [3] Cutadapt  ──  optional (config["cutadapt"]["use"])
#      └─ Removes primers from paired-end reads.
#         Parameters:
#           · threads         : Number of CPU threads
#           · primer_forward  : Forward primer sequence (-g)
#           · primer_reverse  : Reverse primer sequence (-G)
#           · primer_overlap  : Minimum overlap between primer and read (-O)
#           · primer_err      : Maximum error rate allowed in primer matching (-e)
#
#  [4] Trimmomatic  ──  optional (config["trimmomatic"]["use"])
#      └─ Trims low-quality bases from paired-end reads (PE mode).
#         Parameters:
#           · threads               : Number of CPU threads
#           · leading               : Min quality to keep a leading base (LEADING)
#           · trailing              : Min quality to keep a trailing base (TRAILING)
#           · slidingwindow_size    : Sliding window size (SLIDINGWINDOW)
#           · slidingwindow_quality : Min average quality in the sliding window
#           · minlen                : Minimum read length to retain (MINLEN)
#
#  [5] FastQC (Cleaned)
#      └─ Quality assessment of reads after cleaning.
#         Automatically adapts to whichever step(s) were run.
#
#  [6] MultiQC (Cleaned)
#      └─ Aggregates all cleaned FastQC reports into a single
#         interactive HTML summary report.
#
#  ────────────── CONDITIONAL FOLDER STRUCTURE ───────────────
# 
#  Cutadapt only      →  02_Cutadapt/    →  03_CleanedFastQC/
#  Trimmomatic only   →  02_Trimmomatic/ →  03_CleanedFastQC/
#  Both               →  02_Cutadapt/    →  03_Trimmomatic/ →  04_CleanedFastQC/
#
#  ──────────────────────── OUTPUTS ──────────────────────────
#
#  Reports  →  REPORT + "/00_Clean/multiqc_raw.html"
#              REPORT + "/00_Clean/multiqc_cleaned.html"
#  Reads    →  FASTQ_FOLDER + "/{sample}_{type}" + FASTQ_SUFFIX + ".fastq.gz"
#              (path exported for use in downstream rules)
#
# ============================================================
# FastQC to check the reads quality
# ============================================================

rule fastqc_raw:
    input:
        fastq = DATA + "{sample}_{type}.fastq.gz",
    output:
        fastqc_report_html = OUTPUT + "00_Clean/01_RawFastQC/{sample}_{type}_fastqc.html",
        fastqc_report_zip  = OUTPUT + "00_Clean/01_RawFastQC/{sample}_{type}_fastqc.zip",
    conda:
        CONTAINER + "00_Clean.yaml"
    threads:
        config["fastqc"]["threads"]
    message:
        "Creating FastQC reports for {wildcards.sample} raw reads."
    shell:
        """
        mkdir -p $(dirname {output.fastqc_report_zip})
        fastqc \
            --threads {threads} \
            --outdir  $(dirname {output.fastqc_report_zip}) \
            {input.fastq} 
        """

# ============================================================
# Summarize FastQC results in a single report
# ============================================================

rule multiqc_raw:
    input:
        fastqc_reports_zip = expand(OUTPUT + "00_Clean/01_RawFastQC/{sample}_{type}_fastqc.zip", sample = SAMPLES, type = TYPES)
    output:
        multiqc_report = REPORT + "00_Clean/multiqc_raw.html"
    conda:
        CONTAINER + "00_Clean.yaml"
    message:
        "Aggregating FastQC reports of raw reads with MultiQC."
    shell:
        """
        multiqc \
            --no-data-dir \
            --outdir      $(dirname {output.multiqc_report}) \
            --filename    $(basename {output.multiqc_report}) \
            {input.fastqc_reports_zip}
		"""

# ============================================================
# Cutadapt to remove primers 
# ============================================================
def build_primer_flags(flag, value):
    # Build cutadapt primer flag(s) from config value.
    # Handles: single string, list of strings, empty string, or None.
    # Returns an empty string if no primer is set (flag is silently skipped).
    if not value:
        return ""
    if isinstance(value, list):
        return " ".join(f"{flag} {v}" for v in value if v)
    return f"{flag} {value}"

if CUTADAPT:

    rule cutadapt:
        input:
            fastq_1 = DATA + "{sample}_" + TYPES[0] + ".fastq.gz",
            fastq_2 = DATA + "{sample}_" + TYPES[1] + ".fastq.gz"
        output:
            fastq_1_cutadapt = OUTPUT + "00_Clean/02_Cutadapt/{sample}_" + TYPES[0] + "_Cutadapt.fastq.gz",
            fastq_2_cutadapt = OUTPUT + "00_Clean/02_Cutadapt/{sample}_" + TYPES[1] + "_Cutadapt.fastq.gz",
        params:
            primer_5_fwd   = lambda w: build_primer_flags("-g", config["cutadapt"].get("primer_5_forward")),
            primer_5_rev   = lambda w: build_primer_flags("-G", config["cutadapt"].get("primer_5_reverse")),
            primer_3_fwd   = lambda w: build_primer_flags("-a", config["cutadapt"].get("primer_3_forward")),
            primer_3_rev   = lambda w: build_primer_flags("-A", config["cutadapt"].get("primer_3_reverse")),
            primer_overlap = config["cutadapt"]["primer_overlap"],
            primer_err     = config["cutadapt"]["primer_err"],
        threads:
            config["cutadapt"]["thread"]
        conda:
            CONTAINER + "00_Clean.yaml"
        message:
            "Removing primers of {wildcards.sample} with Cutadapt"
        shell:
            """
            cutadapt \
                {params.primer_5_fwd} \
                {params.primer_5_rev} \
                {params.primer_3_fwd} \
                {params.primer_3_rev} \
                -j {threads} \
                -O {params.primer_overlap} \
                -e {params.primer_err} \
                -o {output.fastq_1_cutadapt} \
                -p {output.fastq_2_cutadapt} \
                --match-read-wildcards \
                {input.fastq_1} \
                {input.fastq_2} \
            """

# ============================================================
# Trimmomatic to trim bad quality reads 
# ============================================================
SEPARATOR = "_"

if TRIMMOMATIC:

    TRIMMED_FOLDER  = "02_Trimmomatic/"
    SAMPLE_FOLDER   = DATA
    CUTADAPT_SUFFIX = ""

    if CUTADAPT:

        TRIMMED_FOLDER  = "03_Trimmomatic/"
        CUTADAPT_SUFFIX = "_Cutadapt"
        SAMPLE_FOLDER   = OUTPUT + "00_Clean/02_Cutadapt/"
        SEPARATOR       = ""

    rule trimmomatic:
        input:
            fastq_1 = SAMPLE_FOLDER + "{sample}_" + TYPES[0] + CUTADAPT_SUFFIX + ".fastq.gz",
            fastq_2 = SAMPLE_FOLDER + "{sample}_" + TYPES[1] + CUTADAPT_SUFFIX + ".fastq.gz"
        output:
            fastq_1_trimmed  = OUTPUT + "00_Clean/" + TRIMMED_FOLDER + "{sample}_" + TYPES[0] + CUTADAPT_SUFFIX + SEPARATOR + "Trimmed.fastq.gz",
            fastq_2_trimmed  = OUTPUT + "00_Clean/" + TRIMMED_FOLDER + "{sample}_" + TYPES[1] + CUTADAPT_SUFFIX + SEPARATOR + "Trimmed.fastq.gz",
            fastq_1_unpaired = OUTPUT + "00_Clean/" + TRIMMED_FOLDER + "{sample}_" + TYPES[0] + CUTADAPT_SUFFIX + SEPARATOR + "Unpaired.fastq.gz",
            fastq_2_unpaired = OUTPUT + "00_Clean/" + TRIMMED_FOLDER + "{sample}_" + TYPES[1] + CUTADAPT_SUFFIX + SEPARATOR + "Unpaired.fastq.gz",
        params:
            leading               = config["trimmomatic"]["leading"],
            trailing              = config["trimmomatic"]["trailing"],
            slidingwindow_size    = config["trimmomatic"]["slidingwindow_size"],
            slidingwindow_quality = config["trimmomatic"]["slidingwindow_quality"],
            minlen                = config["trimmomatic"]["minlen"]
        threads:
            config["trimmomatic"]["thread"],
        conda:
            CONTAINER + "00_Clean.yaml"
        message:
            "Trimming {wildcards.sample} raw reads with Trimmomatic"
        shell:
            """
            trimmomatic \
                PE \
                -threads {threads} \
                {input.fastq_1} \
                {input.fastq_2} \
                {output.fastq_1_trimmed} \
                {output.fastq_1_unpaired} \
                {output.fastq_2_trimmed} \
                {output.fastq_2_unpaired} \
                LEADING:{params.leading} \
                TRAILING:{params.trailing} \
                SLIDINGWINDOW:{params.slidingwindow_size}:{params.slidingwindow_quality} \
                MINLEN:{params.minlen}
            """

# ============================================================
# FastQC to check the trimmed reads quality
# ============================================================
FASTQ_FOLDER  = DATA
FASTQ_SUFFIX  = ""
FASTQC_FOLDER = "03_CleanedFastQC"

if CUTADAPT:
   
    FASTQ_FOLDER  = OUTPUT + "00_Clean/02_Cutadapt/"
    FASTQ_SUFFIX  = "_Cutadapt"

if TRIMMOMATIC:
    
    FASTQ_FOLDER  = OUTPUT + "00_Clean/" + TRIMMED_FOLDER
    FASTQ_SUFFIX  = FASTQ_SUFFIX + SEPARATOR + "Trimmed"


if CUTADAPT and TRIMMOMATIC:

    FASTQC_FOLDER = "04_CleanedFastQC"

rule fastqc_cleaned:
    input:
        fastq = FASTQ_FOLDER + "{sample}_{type}" + FASTQ_SUFFIX + ".fastq.gz",
    output:
        fastqc_report_html = OUTPUT + "00_Clean/" + FASTQC_FOLDER + "/{sample}_{type}" + FASTQ_SUFFIX + "_fastqc.html",
        fastqc_report_zip  = OUTPUT + "00_Clean/" + FASTQC_FOLDER + "/{sample}_{type}" + FASTQ_SUFFIX + "_fastqc.zip",
    conda:
        CONTAINER + "00_Clean.yaml"
    threads:
        config["fastqc"]["threads"]
    message:
        "Creating FastQC reports for {wildcards.sample} trimmed reads."
    shell:
        """
        mkdir -p $(dirname {output.fastqc_report_zip})
        fastqc \
            --threads {threads} \
            --outdir  $(dirname {output.fastqc_report_zip}) \
            {input.fastq}
        """

# ============================================================
# MultiQC to check the reads quality after cleaning
# ============================================================

rule multiqc_cleaned:
    input:
        fastqc_report_zip = expand(OUTPUT + "00_Clean/{fastqc_folder}/{sample}_{type}" + FASTQ_SUFFIX + "_fastqc.zip", sample = SAMPLES, type = TYPES, fastqc_folder = FASTQC_FOLDER)
    output:
        multiqc_report = REPORT + "00_Clean/multiqc_cleaned.html"
    conda:
        CONTAINER + "00_Clean.yaml"
    message:
        "Aggregating FastQC reports of trimmed reads with MultiQC."
    shell:
        """		
        multiqc \
            --no-data-dir \
            --outdir   $(dirname {output.multiqc_report}) \
            --filename $(basename {output.multiqc_report}) \
            {input.fastqc_report_zip}
		"""

# ============================================================
# Input FastQ to use for downstream analysis
# ============================================================
# This the path to use for the following rule if you need the cleaned reads
# fastq_1 = FASTQ_FOLDER + "{sample}_" + TYPES[0] + FASTQ_SUFFIX + ".fastq.gz",
# fastq_2 = FASTQ_FOLDER + "{sample}_" + TYPES[1] + FASTQ_SUFFIX + ".fastq.gz"