# ============================================================
#  01_Alignment : Mapping reads into the reference genome
# ============================================================
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
		fasta = REF + config["reference"]["ref_fasta"],
		gtf   = REF + "formated_" + config["reference"]["ref_gtf"]
	output:
		index_done = OUTPUTDIR + "02_Alignement/01_Index/.done"
	params:
		index_dir    = OUTPUTDIR + "02_STAR/01_Index/",
		sjdbOverhang = config["star"]["sjdbOverhang"],
		cpus         = config["star"]["cpus_index"],
		ram          = config["star"]["ram_index"]
	conda:
		CONTAINER + "01_Alignment.yaml"
	shell:
		"""
		mkdir -p {params.index_dir}
		STAR \
			--runThreadN {params.cpus} \
			--runMode genomeGenerate \
			--genomeDir {params.index_dir} \
			--genomeFastaFiles {input.fasta} \
			--sjdbGTFfile {input.gtf} \
			--sjdbOverhang {params.sjdbOverhang} \
			--limitGenomeGenerateRAM {params.ram}
		touch {output.index_done}
		""" 

# ============================================================
# STAR alignement
# ============================================================