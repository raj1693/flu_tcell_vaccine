# Load necessary libraries
library(ape)         # Phylogenetic tree construction
library(phangorn)    # Distance matrix & clustering
library(Biostrings)  # Handling protein sequences

# Define input and output files (Update paths)
fasta_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/NP_H3/combined.fasta"
output_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/NP_H3/NP_H3_consensus.fasta"

# Read protein sequences
protein_sequences <- readAAStringSet(fasta_file)

# Save sequences to a temporary file for alignment
temp_fasta <- tempfile(fileext = ".fasta")
writeXStringSet(protein_sequences, temp_fasta)

# Define output file for aligned sequences
aligned_fasta <- tempfile(fileext = ".fasta")

# Run MAFFT for sequence alignment
mafft_cmd <- paste("mafft --auto", temp_fasta, ">", aligned_fasta)
system(mafft_cmd)

# Read aligned sequences
aligned_sequences <- readAAStringSet(aligned_fasta)

# Convert sequences to character format
seq_list <- as.character(aligned_sequences)

# Convert to matrix (one row per sequence, one column per residue)
seq_matrix <- do.call(rbind, strsplit(seq_list, ""))

# Ensure sequences are properly aligned
seq_lengths <- unique(nchar(seq_list))
if (length(seq_lengths) > 1) {
  stop("Error: Sequences are not properly aligned. Check MAFFT output.")
}

# Convert to phyDat format for phylogenetic analysis
seq_matrix[seq_matrix %in% c("X", "?", "*")] <- "-"  # Replace unknowns with gaps
phyDat_sequences <- phyDat(seq_matrix, type = "AA")

# Compute distance matrix using JTT model
dist_matrix <- dist.ml(phyDat_sequences, model = "JTT")

# Build a phylogenetic tree using UPGMA (always ultrametric)
upgma_tree <- upgma(dist_matrix)
hc_tree <- as.hclust(upgma_tree)

# Cut tree into k clusters (Adjust k as needed)
k <- 3  
clusters <- cutree(hc_tree, k = k)

# Find the largest cluster (most sequences)
largest_cluster <- names(which.max(table(clusters)))

# Extract sequences from the largest cluster
cluster_seqs <- seq_matrix[clusters == largest_cluster, , drop = FALSE]

# Compute consensus sequence for the largest cluster
if (nrow(cluster_seqs) == 1) {
  consensus_seq <- paste(cluster_seqs, collapse = "")
} else {
  consensus_seq <- apply(cluster_seqs, 2, function(column) {
    column <- column[column != "-"]  # Ignore gaps
    if (length(column) == 0) return("-")  # Keep gap if all positions are gaps
    return(names(which.max(table(column))))  # Most frequent residue
  })
  consensus_seq <- paste(consensus_seq, collapse = "")
}

# Remove excessive gaps to fix length
consensus_seq <- gsub("-+", "-", consensus_seq)  # Collapse consecutive gaps
consensus_seq <- gsub("^-|-$", "", consensus_seq)  # Trim leading/trailing gaps

# Convert to AAString
final_microconsensus <- AAString(consensus_seq)

# Write final microconsensus to a FASTA file
sink(output_file)
cat(">microconsensus_phylo\n", as.character(final_microconsensus), "\n", sep = "")
sink()

message(paste("Microconsensus sequence saved to:", output_file))
