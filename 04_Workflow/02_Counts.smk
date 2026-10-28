# ============================================================
#  WORKFLOW : Expression Quantification Pipeline
# ============================================================
#
#  This workflow quantifies gene expression from aligned reads
#  and produces filtered count matrices ready for downstream
#  differential expression analysis.
#
#  ─────────────────── PIPELINE OVERVIEW  ────────────────────
#
#  [1] FeatureCounts  ──  per sample
#      └─ Counts reads overlapping genomic features (genes)
#         from sorted BAM files using a GTF annotation.
#         Parameters:
#           · type          : Feature type to count over (-t) e.g. exon
#           · group         : Attribute to group features by (-g) e.g. gene_id
#           · strand        : Strandedness (-s): 0 = unstranded, 1 = forward, 2 = reverse
#           · B             : Both mates must map (-B)
#           · C             : Don't count chimeric fragments (-C)
#           · Q             : Minimum mapping quality (-Q)
#           · multimapped   : Count multimapped reads (-M)
#           · fraction      : Fractional counting for multimappers (--fraction, requires -M)
#           · extraAttribute: Additional GTF attribute to report e.g. gene_name
#         Note: -p and --countReadPairs are always enabled
#         (paired-end mode). BAM path and suffix are stripped
#         from the output header for cleaner sample names and
#         DESeq2 compatibility.
#
#  [2] Count Matrix
#      └─ Aggregates all per-sample FeatureCounts outputs into
#         a single gene × sample raw count matrix.
#         Script: 02_Counts/01_Count_matrix.R
#
#  [3] CPM Matrix
#      └─ Normalizes the raw count matrix into CPM
#         (Counts Per Million) values.
#         Script: 02_Counts/02_CPM_matrix.R
#
#  [4] CPM Filter
#      └─ Removes lowly expressed genes based on a CPM
#         threshold. Produces both filtered and removed
#         matrices for counts and CPM.
#         Parameters:
#           · threshold : Minimum CPM value to retain a gene
#         Script: 02_Counts/03_CPM_filter.R
#
#  [5] MultiQC (Counts)
#      └─ Aggregates all FeatureCounts summary files into a
#         single interactive HTML report.
#
#  ──────────────────────── OUTPUTS ──────────────────────────
#
#  Per sample  →  02_Counts/01_Samples/{sample}_expression_counts.tsv
#                 02_Counts/01_Samples/{sample}_expression_counts.summary
#
#  Matrices    →  02_Counts/02_Matrix/count_matrix.tsv
#                 02_Counts/02_Matrix/cpm_matrix.tsv
#                 02_Counts/02_Matrix/filtered_count_matrix.tsv
#                 02_Counts/02_Matrix/filtered_cpm_matrix.tsv
#                 02_Counts/02_Matrix/removed_count_matrix.tsv
#                 02_Counts/02_Matrix/removed_cpm_matrix.tsv
#
#  Report      →  REPORT + /02_Counts/multiqc_featurecounts.html
#
# ============================================================
#  Quantifying expression with FeatureCounts
# ============================================================

rule featurecounts:
    input:
        bam = OUTPUT + "01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.bam",
        gtf  = REFERENCE + config["reference"]["ref_gtf"]
    output:
        counts  = OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.tsv",
        summary = OUTPUT + "02_Counts/01_Samples/{sample}_expression_counts.summary"
    params:
        type           = config["featurecounts"]["type"],
        group          = config["featurecounts"]["group"],
        strand         = config["featurecounts"]["strand"],
        B              = ("-B" if config["featurecounts"]["B"] else ""),
        C              = ("-C" if config["featurecounts"]["C"] else ""),
        Q              = ("-Q" if config["featurecounts"]["Q"] else ""),
        multimapped    = ("-M" if config["featurecounts"]["multimapped"] else ""),
        fraction       = ("--fraction" if config["featurecounts"]["fraction"] else ""),
        extraAttribute = ("--extraAttribute " + config["featurecounts"]["extraAttribute"] if config["featurecounts"]["extraAttribute"] else "")
    threads:
        config["featurecounts"]["threads"]
    log:
        OUTPUT + "02_Counts/01_Samples/log/{sample}.txt"
    conda:
        CONTAINER + "02_Counts.yaml"
    message:
        "Quantifying gene expression from BAMs for {wildcards.sample}"
    shell:
        """
        featureCounts \
            -p \
            -B \
            -C \
            --countReadPairs \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -t {params.type} \
            -g {params.group} \
            -s {params.strand} \
            {params.B} \
            {params.C} \
            {params.Q} \
            {params.multimapped} \
            {params.fraction} \
            {params.extraAttribute} \
            {input.bam} \
            > {log} 2>&1

        mv {output.counts}.summary {output.summary}

        BAMDIR=$(dirname {input.bam})/
        sed -e "s|$BAMDIR||g" \
            -e 's/_Aligned.sortedByCoord.out.bam//g' \
            {output.counts} > {output.counts}.temp && mv {output.counts}.temp {output.counts}
        """

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