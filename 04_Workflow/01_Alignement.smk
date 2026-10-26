# ============================================================
#  01_Alignment : Mapping reads into the reference genome
# ============================================================
#
#  This workflow downloads and prepares a reference genome,
#  builds a STAR index, aligns cleaned paired-end reads, and
#  produces an aggregated alignment quality report.
#
#  ──────────────────────── PIPELINE ─────────────────────────
#
#  [1] Download Reference
#      └─ Fetches the reference genome FASTA and GTF annotation
#         files from remote URLs defined in config.
#           · link_ref_fasta : URL for the genome FASTA file
#           · link_ref_gtf   : URL for the genome GTF file
#
#  [2] Gunzip Reference
#      └─ Decompresses the downloaded .gz FASTA and GTF files
#         into the reference directory for downstream use.
#
#  [3] STAR Index
#      └─ Builds a STAR genome index from the reference FASTA
#         and GTF files.
#         Parameters:
#           · threads      : Number of CPU threads
#           · sjdbOverhang : Overhang length for splice junctions (typically read length − 1)
#           · ram_limit    : Maximum RAM allowed for index generation
#
#  [4] STAR Alignment
#      └─ Aligns cleaned paired-end reads to the indexed genome.
#         Outputs coordinate-sorted BAM files.
#         Parameters:
#           · threads               : Number of CPU threads
#           · outFilterMultimapNmax : Max number of loci a read can map to (multimappers filter)
#           · outSAMstrandField     : Strand field in SAM output
#           · outSAMattributes      : SAM attributes to include
#           · readFilesCommand      : Command to decompress input reads (e.g. zcat)
#           · ram_limit             : Maximum RAM for BAM sorting
#           · twopass               : STAR two-pass mapping mode
#
#  [5] BAM Indexation
#      └─ Index BAM files produced by the STAR alignement
#
#  [6] MultiQC (Alignment)
#      └─ Aggregates all STAR final log files into a single
#         interactive HTML summary report.
#
#  ──────────────────────── FOLDER STRUCTURE ─────────────────
#
#  01_Index/  →  STAR genome index (SA, SAindex, …)
#  02_BAM/    →  Coordinate-sorted BAM per sample
#  03_Log/    →  Per-sample STAR logs (used by MultiQC)
#
#  ──────────────────────── OUTPUTS ──────────────────────────
#
#  Index   →  OUTPUT + "/01_Alignment/01_Index/"  (SA, SAindex, etc...)
#  BAMs     →  OUTPUT + "/01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.bam"
#  BAIs     →  OUTPUT + "/01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.baI"
#  Report  →  REPORT + "/01_Alignment/multiqc_alignement.html"
#
# ============================================================
# Download and gunzip of reference genome
# ============================================================

rule dl_reference:
    output:
        fasta = REFERENCE + config["reference"]["ref_fasta"] + ".gz",
        gtf   = REFERENCE + config["reference"]["ref_gtf"] + ".gz"
    params:
        link_ref_fasta = config["reference"]["link_ref_fasta"],
        link_ref_gtf   = config["reference"]["link_ref_gtf"]
    conda:
        CONTAINER + "01_Alignment.yaml"
    message:
        "Downloading reference genome's .fasta and .gtf files."
    shell:
        """
        wget -O {output.fasta} {params.link_ref_fasta}
        wget -O {output.gtf} {params.link_ref_gtf}
        """


rule gunzip_reference:
    input:
        fasta = REFERENCE + config["reference"]["ref_fasta"] + ".gz",
        gtf   = REFERENCE + config["reference"]["ref_gtf"] + ".gz"
    output:
        fasta = REFERENCE + config["reference"]["ref_fasta"],
        gtf   = REFERENCE + config["reference"]["ref_gtf"]
    conda:
        CONTAINER + "01_Alignment.yaml"
    message:
        "Unzipping reference genome's .fasta and .gtf files."
    shell:
        """
        gunzip -c {input.fasta} > {output.fasta}
        gunzip -c {input.gtf} > {output.gtf}
        """


# ============================================================
# STAR indexation
# ============================================================

rule star_index:
    input:
        fasta = REFERENCE + config["reference"]["ref_fasta"],
        gtf   = REFERENCE + config["reference"]["ref_gtf"]
    output:
        genome_sa  = OUTPUT + "01_Alignment/01_Index/SA",
        genome_sai = OUTPUT + "01_Alignment/01_Index/SAindex"
    params:
        sjdbOverhang = config["star_index"]["sjdbOverhang"],
        ram_limit    = config["star_index"]["ram_limit"],
        genome_name  = config["reference"]["name"]
    threads:
        config["star_index"]["threads"]
    conda:
        CONTAINER + "01_Alignment.yaml"
    message:
        "Indexing {params.genome_name} Genome."
    shell:
        """
        mkdir -p $(dirname {output.genome_sa})
        STAR \
            --runThreadN             {threads} \
            --genomeDir              $(dirname {output.genome_sa}) \
            --genomeFastaFiles       {input.fasta} \
            --sjdbGTFfile            {input.gtf} \
            --sjdbOverhang           {params.sjdbOverhang} \
            --limitGenomeGenerateRAM {params.ram_limit} \
            --runMode                genomeGenerate
        """

# ============================================================
# STAR alignement
# ============================================================

rule star:
    input:
        fastq_1    = FASTQ_FOLDER + "{sample}_" + TYPES[0] + FASTQ_SUFFIX + ".fastq.gz",
        fastq_2    = FASTQ_FOLDER + "{sample}_" + TYPES[1] + FASTQ_SUFFIX + ".fastq.gz",
        genome_sa  = OUTPUT + "01_Alignment/01_Index/SA",
        genome_sai = OUTPUT + "01_Alignment/01_Index/SAindex"
    output:
        OUTPUT + "01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        index_dir             = OUTPUT + "01_Alignment/01_Index/",
        outprefix             = OUTPUT + "01_Alignment/03_Log/{sample}/",
        outFilterMultimapNmax = config["star_alignment"]["outFilterMultimapNmax"],
        outSAMstrandField     = config["star_alignment"]["outSAMstrandField"],
        outSAMattributes      = config["star_alignment"]["outSAMattributes"],
        readFilesCommand      = config["star_alignment"]["readFilesCommand"],
        ram_limit             = config["star_alignment"]["ram_limit"],
        twopass               = config["star_alignment"]["twopass"],
        genome_name           = config["reference"]["name"]
    threads:
        config["star_alignment"]["threads"]
    conda:
        CONTAINER + "01_Alignment.yaml"
    resources:
        star_slot = 1
    message:
        "Aligning {wildcards.sample} on {params.genome_name}."
    shell: 
        """
        mkdir -p {params.outprefix}
        mkdir -p $(dirname {output})
        STAR \
            --runThreadN            {threads} \
            --genomeDir             {params.index_dir} \
            --readFilesIn           {input.fastq_1} {input.fastq_2} \
            --readFilesCommand      {params.readFilesCommand} \
            --outFilterMultimapNmax {params.outFilterMultimapNmax} \
            --outFilterType         BySJout \
            --outFileNamePrefix     {params.outprefix}{wildcards.sample}_ \
            --outSAMstrandField     {params.outSAMstrandField} \
            --outSAMtype            BAM SortedByCoordinate \
            --outSAMattributes      {params.outSAMattributes} \
            --limitBAMsortRAM       {params.ram_limit} \
            --twopassMode           {params.twopass} 
        mv {params.outprefix}{wildcards.sample}_Aligned.sortedByCoord.out.bam {output}
        """

# ============================================================
# Indexing BAM files to obtain BAI files 
# ============================================================

if BAM_INDEX:

    rule bam_index:
        input:
            OUTPUT + "01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.bam"
        output:
            OUTPUT + "01_Alignment/02_BAM/{sample}_Aligned.sortedByCoord.out.bai"
        conda:
            CONTAINER + "01_Alignment.yaml"
        message:
            "Indexing {wildcards.sample}'s BAM file"
        shell:
            """
            samtools index {input} {output}
            """


# ============================================================
# MultiQC report for STAR alignement
# ============================================================

rule multiqc_alignment:
    input:
        star_logs = expand( OUTPUT + "01_Alignment/03_Log/{sample}/{sample}_Log.final.out", sample = SAMPLES)
    output:
        multiqc_report = REPORT + "01_Alignment/multiqc_alignment.html"
    conda:
        CONTAINER + "00_Clean.yaml"
    message:
        "Aggregating STAR alignment reports with MultiQC."
    shell:
        """
        multiqc \
            --no-data-dir \
            --outdir   $(dirname {output.multiqc_report}) \
            --filename $(basename {output.multiqc_report}) \
            {input.star_logs}
        """