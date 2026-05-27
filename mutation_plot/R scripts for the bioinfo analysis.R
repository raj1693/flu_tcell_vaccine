library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

# Load sequences from a FASTA file
fasta_file <- "NP_Primary_consensus_H1N1" # Replace with your FASTA file path
sequences <- read.dna(fasta_file, format = "fasta", as.character = TRUE)
sample_names <- rownames(sequences)

# Convert sequences to amino acids
aa_sequences <- apply(sequences, 1, paste, collapse = "")
aa_matrix <- t(sapply(aa_sequences, strsplit, split = ""))

# Create a distance matrix based on sequence differences
dist_matrix <- dist(aa_matrix, method = "manhattan")

# Build a phylogenetic tree
phylo_tree <- nj(dist_matrix)

# Prepare data for the mutation heatmap
aa_colors <- c("A" = "red", "R" = "blue", "N" = "green", "D" = "orange", 
               "C" = "purple", "Q" = "pink", "E" = "yellow", "G" = "cyan",
               "H" = "brown", "I" = "black", "L" = "grey", "K" = "magenta",
               "M" = "gold", "F" = "navy", "P" = "darkgreen", "S" = "darkred",
               "T" = "lightblue", "W" = "darkblue", "Y" = "darkorange",
               "V" = "darkgrey", "-" = "white") # Add more if needed

heatmap_data <- data.frame(
  Sample = rep(rownames(aa_matrix), each = ncol(aa_matrix)),
  Position = rep(1:ncol(aa_matrix), times = nrow(aa_matrix)),
  AminoAcid = as.vector(t(aa_matrix))
)

heatmap_data$Color <- aa_colors[heatmap_data$AminoAcid]

# Plot the phylogenetic tree
tree_plot <- ggtree(phylo_tree, layout = "rectangular") + 
  theme_tree2()

# Plot the mutation heatmap
heatmap_plot <- ggplot(heatmap_data, aes(x = Position, y = Sample, fill = AminoAcid)) +
  geom_tile() +
  scale_fill_manual(values = aa_colors, na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Position", y = NULL)

# Combine the plots
library(patchwork)
tree_plot + heatmap_plot + plot_layout(widths = c(1, 3))
