# ============================================================
#  03_CPM_filter: Filtering count matrix with a cpm threshold
# ============================================================

# Loading library
suppressPackageStartupMessages(library(dplyr))

# Loading snakemake variables
input.metadata         <- snakemake@input[["metadata"]]
input.count_matrix     <- snakemake@input[["count_matrix"]]
input.cpm_matrix       <- snakemake@input[["cpm_matrix"]]
output.count_matrix    <- snakemake@output[["count_matrix"]]
output.cpm_matrix      <- snakemake@output[["cpm_matrix"]]
output.rm_count_matrix <- snakemake@output[["rm_count_matrix"]]
output.rm_cpm_matrix   <- snakemake@output[["rm_cpm_matrix"]]
params.threshold       <- snakemake@params[["threshold"]]

# ============================================================
# Loading files
# ============================================================

metadata     <- read.delim(input.metadata,     header = TRUE, sep = "\t")
count_matrix <- read.delim(input.count_matrix, header = TRUE, sep = "\t")
cpm_matrix   <- read.delim(input.cpm_matrix,   header = TRUE, sep = "\t")

# ============================================================
#  Filtering genes based on the CPM threshold
# ============================================================

# Every conditions of the metadata
conditions <- unique(metadata$Condition)

# This function check if a gene is expressed above the CPM 
# threshold in all samples of at least one of the conditions
filter_gene <- function(row_counts, 
                        metadata,
                        conditions,
                        threshold
                        ){
  
  sapply(conditions, function(x){
    
    # Retrieving every samples counts of one condition
    sample_counts <- row_counts[metadata[metadata$Condition == x, ]$ID]
    # Checking if every samples pass the CPM threshold
    all(sample_counts > params.threshold)
    
    }) %>% 
    any() %>%
    # Checking if any condition pass the filtering
    return()
}

# Check expression and a list of boolean if the gene pass the 
# filtering (TRUE) or not (FALSE)
kept_genes <- apply(cpm_matrix[,-1],
                    MARGIN = 1,
                    FUN    = filter_gene,
                    conditions = conditions,
                    metadata   = metadata,
                    threshold  = params.threshold
                    )

# Filtering count and cpm matrix, and keeping the genes removed
# separately
kept_genes_count_matrix    <- count_matrix[kept_genes, ]
kept_genes_cpm_matrix      <- cpm_matrix[kept_genes, ]
removed_genes_count_matrix <- count_matrix[!kept_genes, ]
removed_genes_cpm_matrix   <- cpm_matrix[!kept_genes, ]

# ============================================================
# Output
# ============================================================
write.table(kept_genes_count_matrix,    file = output.count_matrix,    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(kept_genes_cpm_matrix,      file = output.cpm_matrix,      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(removed_genes_count_matrix, file = output.rm_count_matrix, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(removed_genes_cpm_matrix,   file = output.rm_cpm_matrix,   sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




