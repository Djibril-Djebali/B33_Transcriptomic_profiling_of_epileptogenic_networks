# ============================================================
#  01_Count_matrix : Aggregating FeatureCounts results into a 
#  counts matrix
# ============================================================

# Loading snakemake variables
input.counts  <- snakemake@input[["counts"]]
output.counts <- snakemake@output[["counts"]]

# ============================================================
# Creating and outputing the count matrix
# ============================================================

# Retrieve all counts from each samples and binding them together
for (i in 1:length(input.counts)) {
  sample_counts <- read.delim(input.counts[[i]], header = FALSE, sep = "\t", comment.char = "#", quote = "")
  if(i == 1){
    matrix_counts <- sample_counts[, c(1, 7)]
  }else{
    matrix_counts <- cbind(matrix_counts, sample_counts[, 7])
  }
}

# ============================================================
# Output
# ============================================================
write.table(matrix_counts, file = output.counts, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)