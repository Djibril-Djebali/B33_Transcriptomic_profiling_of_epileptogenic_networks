# ============================================================
#  02_CPM_matrix : Convertion count into cpm
# ============================================================

# Loading library
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(edgeR))

# Loading snakemake variables
input.counts  <- snakemake@input[["counts"]]
output.counts <- snakemake@output[["counts"]]

# ============================================================
# Loading files
# ============================================================

# Loading count matrix
matrix_counts <- read.delim(input.counts, header = TRUE, sep = "\t", comment.char = "#", quote = "")

# ============================================================
# Converting counts into cpm
# ============================================================

# Keeping gene id
gene_id <- matrix_counts$Geneid

# Convert to class matrix
matrix_counts <- as.matrix(matrix_counts[, -1])
class(matrix_counts) <- "numeric"
# Convert to CPM
matrix_cpm <- cpm(matrix_counts)

# Transform in dataframe
matrix_cpm <- as.data.frame(matrix_cpm)
matrix_cpm <- cbind(Gene_id = gene_id, matrix_cpm)

# Convert to numerical value the cpm
matrix_cpm[] <- lapply(matrix_cpm, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

# ============================================================
# Output
# ============================================================
write.table(matrix_cpm, file = output.counts, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
