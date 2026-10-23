PHONY=help run
MAKEFILE=Makefile

SHELL := /bin/bash

help:
	@echo "Available targets"
	@grep -P '^[0-9A-Za-z]+\:' $(MAKEFILE) | sed 's/:.*//' | sed 's/^/- /'


build:
	@mkdir -p ./.apptainer
	@apptainer pull \
		--disable-cache \
		./.apptainer/snakemake_v9.27.0.sif \
		docker://snakemake/snakemake:v9.27.0

run:
	@apptainer exec ./.apptainer/snakemake_v9.27.0.sif \
		snakemake \
			-s Snakefile \
			--use-conda \
			--sdm apptainer \
			--cores 24 

dry:
	@apptainer exec ./.apptainer/snakemake_v9.27.0.sif \
		snakemake \
			--quiet rules \
			-s Snakefile \
			--dry-run \
			--rerun-triggers mtime 

clean_envs:
	@apptainer exec ./.apptainer/snakemake_v9.27.0.sif \
    	            snakemake \
        	                -s Snakefile \
            	            --dry-run \
                	        --conda-cleanup-envs \
                    	    --container-cleanup-images \
                        	--rerun-triggers mtime