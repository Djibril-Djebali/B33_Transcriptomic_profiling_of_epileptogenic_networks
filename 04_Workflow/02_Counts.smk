# ============================================================
#  02_Counts: Quantification of expression from the mapping 
#  results
# ============================================================
#
# ============================================================
#  Quantifying expression with FeatureCounts
# ============================================================

# rule featurecounts:
#     input:
#         bam = OUTPUT + "01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.bam",
#         gtf  = REFERENCE + config["reference"]["ref_gtf"]
#     output:
#         counts  = OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.tsv",
#         summary = OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.summary"
#     params:
#         type           = config["featurecounts"]["type"],
#         group          = config["featurecounts"]["group"],
#         strand         = config["featurecounts"]["strand"],
#         B              = ("-B" if config["featurecounts"]["B"] else ""),
#         C              = ("-C" if config["featurecounts"]["C"] else ""),
#         Q              = ("-Q" if config["featurecounts"]["Q"] else ""),
#         multimapped    = ("-M" if config["featurecounts"]["multimapped"] else ""),
#         fraction       = ("--fraction" if config["featurecounts"]["fraction"] else ""),
#         extraAttribute = ("--extraAttribute " + config["featurecounts"]["extraAttribute"] if config["featurecounts"]["extraAttribute"] else "")
#     threads:
#         config["featurecounts"]["threads"]
#     log:
#         OUTPUT + "02_Counts/01_Samples/log/{sample}.txt"
#     conda:
#         CONTAINER + "02_Counts.yaml"
#     message:
#         "Quantifying gene expression from BAMs for {wildcards.sample}"
#     shell:
#         """
#         featureCounts \
#             -p \
#             -B \
#             -C \
#             --countReadPairs \
#             -T {threads} \
#             -a {input.gtf} \
#             -o {output.counts} \
#             -t {params.type} \
#             -g {params.group} \
#             -s {params.strand} \
#             {params.B} \
#             {params.C} \
#             {params.Q} \
#             {params.multimapped} \
#             {params.fraction} \
#             {params.extraAttribute} \
#             {input.bam} \
#             > {log} 2>&1

#         mv {output.counts}.summary {output.summary}

#         BAMDIR=$(dirname {input.bam})/
#         sed -e "s|$BAMDIR||g" \
#             -e 's/_Aligned.sortedByCoord.out.bam//g' \
#             {output.counts} > {output.counts}.temp && mv {output.counts}.temp {output.counts}
#         """

# ============================================================
#  Aggregating FeatureCounts results into a counts matrix
# ============================================================

rule count_matrix:
    input:
        counts = expand(OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.tsv", sample = SAMPLES),
    output:
        counts = OUTPUT + "02_Counts/02_Matrix/count_matrix.tsv"
    conda:
        CONTAINER + "02_Counts.yaml"
    message:
        "Aggregating FeatureCounts results into a counts matrix"
    script:
        SCRIPT + "02_Counts/01_Count_matrix.R"

# ============================================================
#  Converting counts matrix into CPM
# ============================================================

rule cpm_matrix:
    input:
        counts = OUTPUT + "02_Counts/02_Matrix/count_matrix.tsv"
    output:
        counts = OUTPUT + "02_Counts/02_Matrix/cpm_matrix.tsv"
    conda:
        CONTAINER + "02_Counts.yaml"
    message:
        "Converting counts matrix into a CPM matrix"
    script:
        SCRIPT + "02_Counts/02_CPM_matrix.R"

# ============================================================
#  Filtering count matrix with a cpm threshold
# ============================================================

rule cpm_filter:
    input:
        metadata     = CONFIG + config["metadata"],
        cpm_matrix   = OUTPUT + "02_Counts/02_Matrix/cpm_matrix.tsv",
        count_matrix = OUTPUT + "02_Counts/02_Matrix/count_matrix.tsv"
    output:
        cpm_matrix      = OUTPUT + "02_Counts/02_Matrix/filtered_cpm_matrix.tsv",
        count_matrix    = OUTPUT + "02_Counts/02_Matrix/filtered_count_matrix.tsv",
        rm_cpm_matrix   = OUTPUT + "02_Counts/02_Matrix/removed_cpm_matrix.tsv",
        rm_count_matrix = OUTPUT + "02_Counts/02_Matrix/removed_count_matrix.tsv"
    params:
        threshold = config["cpm"]["threshold"]
    conda:
        CONTAINER + "02_Counts.yaml"
    message:
        "Filtering count matrix with a cpm threshold"
    script:
        SCRIPT + "02_Counts/03_CPM_filter.R"

# ============================================================
#  MultiQC report for FeatureCounts output
# ============================================================

rule multiqc_counts:
    input:
        counts_summary = expand(OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.summary", sample = SAMPLES)
    output:
        multiqc_report = REPORT + "02_Counts/multiqc_featurecounts.html"
    conda:
        CONTAINER + "00_Clean.yaml"
    message:
        "Aggregating FeatureCounts summary in MultiQC."
    shell:
        """
        multiqc \
            --no-data-dir \
            --outdir   $(dirname {output.multiqc_report}) \
            --filename $(basename {output.multiqc_report}) \
            {input.counts_summary}
        """