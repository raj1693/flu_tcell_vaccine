# Load necessary libraries
library(parallel)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(readxl)

# Input and output paths
input_file <- "/Users/raj/Desktop/rajesh_t_cell_prediction/iTOL/H3_HA.xlsx"
output_dir <- "/Users/raj/Desktop/rajesh_t_cell_prediction/iTOL/H3_HA"
consensus_file <- file.path(output_dir, "consensus_sequences.fasta")

# Ensure directories exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read Excel file
data <- read_excel(input_file)

# Extract sequences and clades
sequences <- data$Sequence
clades <- data$`Global H3 Clade`
sample_names <- data$`Strain Name`

# Validate input data
if (anyNA(sequences) || anyNA(clades) || anyNA(sample_names)) {
  stop("Error: Missing values detected in input data. Please check your Excel file.")
}

# Get unique clades
unique_clades <- unique(clades)
if (length(unique_clades) == 0) {
  stop("Error: No unique clades found in input data. Please check your Excel file.")
}

# Step 1: Write separate FASTA files for each clade
for (clade in unique_clades) {
  clade_indices <- which(clades == clade)
  clade_sequences <- sequences[clade_indices]
  clade_names <- sample_names[clade_indices]
  
  fasta_file <- file.path(output_dir, paste0("clade_", clade, ".fasta"))
  write.fasta(sequences = as.list(clade_sequences), names = clade_names, file.out = fasta_file)
}

