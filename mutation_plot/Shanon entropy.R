# Load necessary libraries
library(parallel)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(readxl)
library(entropy)

# Input and output paths
input_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/NP_H1.xlsx"  # Replace with your Excel file path
output_dir <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot" # Directory to store output files
aligned_file <- file.path(output_dir, "aligned_sequences.fasta")
dir.create(output_dir, showWarnings = FALSE)

# Read Excel file
data <- read_excel(input_file)

# Extract sequences
sequences <- data$Sequence
sample_names <- data$`Strain Name`

# Step 1: Write sequences to a FASTA file
unaligned_fasta <- file.path(output_dir, "unaligned_sequences.fasta")
write.fasta(sequences = as.list(sequences), 
            names = sample_names, 
            file.out = unaligned_fasta)
cat("Unaligned sequences written to:", unaligned_fasta, "\n")

# Step 2: Align sequences using MAFFT
cmd <- paste("mafft --auto --thread 4", shQuote(unaligned_fasta), ">", shQuote(aligned_file))
cat("Running MAFFT alignment command:", cmd, "\n")

tryCatch({
  system(cmd)
  if (!file.exists(aligned_file)) {
    stop("Alignment failed: Aligned file not created.")
  }
  cat("Aligned sequences written to:", aligned_file, "\n")
}, error = function(e) {
  stop("Error during alignment:", e$message)
})

# Step 3: Load aligned sequences
aligned_sequences <- readAAStringSet(aligned_file)

# Convert aligned sequences to a matrix (rows: sequences, columns: positions)
sequence_matrix <- do.call(rbind, strsplit(as.character(aligned_sequences), split = ""))

# Step 4: Calculate Shannon entropy for each position
num_positions <- ncol(sequence_matrix)

# Function to calculate Shannon entropy for a single position
calculate_entropy <- function(column) {
  probs <- table(column) / length(column)  # Probability of each amino acid
  entropy::entropy(probs, method = "ML")  # Shannon entropy
}

# Calculate entropy for each position
position_entropy <- apply(sequence_matrix, 2, calculate_entropy)

# Create a data frame with entropy values
entropy_df <- data.frame(
  Position = 1:num_positions,
  Shannon_Entropy = position_entropy
)

# Step 5: Save entropy results to a CSV file
entropy_file <- file.path(output_dir, "shannon_entropy.csv")
write.csv(entropy_df, file = entropy_file, row.names = FALSE)

cat("Shannon entropy results saved to:", entropy_file, "\n")

# Step 6: Plot Shannon entropy
ggplot(entropy_df, aes(x = Position, y = Shannon_Entropy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(
    y = "Shannon Entropy"
  ) +
  theme_minimal()
