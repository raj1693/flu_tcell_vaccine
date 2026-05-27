library(ggplot2)
library(dplyr)

# Load aligned protein sequences from the provided FASTA file
fasta_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/M1_H1/consensus_sequences_aligned.fasta" # Path to the uploaded file

# Function to parse the FASTA file
read_fasta <- function(file) {
  lines <- readLines(file)
  seqs <- list()
  current_seq <- NULL
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_seq <- sub("^>", "", line)
      seqs[[current_seq]] <- ""
    } else {
      seqs[[current_seq]] <- paste0(seqs[[current_seq]], line)
    }
  }
  return(seqs)
}

# Parse the sequences
sequences <- read_fasta(fasta_file)

# Extract sequence names and aligned data
Clade_names <- names(sequences)
sequence_data <- t(sapply(sequences, function(seq) strsplit(seq, "")[[1]]))

# Ensure sequences are aligned (all rows should have the same number of columns)
if (!all(nchar(unlist(sequences)) == nchar(sequences[[1]]))) {
  stop("The input sequences are not aligned.")
}

# Reference sequence (the first sequence)
reference_sequence <- sequence_data[1, ]

# Remove the reference sequence from the plot data
Clade_names <- Clade_names[-1]
sequence_data <- sequence_data[-1, , drop = FALSE]

# Compare all sequences to the reference sequence
difference_matrix <- t(sequence_data) != reference_sequence

# Prepare data for the mutation heatmap
heatmap_data <- data.frame(
  Clade = rep(Clade_names, each = ncol(sequence_data)),
  Position = rep(1:ncol(sequence_data), times = nrow(sequence_data)),
  AminoAcid = as.vector(t(sequence_data)),
  IsDifferent = as.vector(difference_matrix)
)

# Define colors for amino acids (ensure all observed amino acids are covered)
unique_amino_acids <- unique(c(as.vector(sequence_data), reference_sequence))

# Dynamically assign colors (adjust color list if needed)
default_colors <- c(
  "red", "blue", "green", "orange", "purple", "pink", "yellow", "cyan",
  "brown", "black", "grey", "magenta", "gold", "navy", "darkgreen",
  "darkred", "lightblue", "darkblue", "darkorange", "darkgrey", "white"
)
amino_acid_keys <- c(
  "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
  "S", "T", "W", "Y", "V", "-", "*"
)

# Ensure keys match all observed amino acids (avoid extra colors)
amino_acid_keys <- unique(amino_acid_keys)  # Remove extra keys if any

# Create the color mapping using only relevant colors
aminoAcidColors <- setNames(default_colors[1:length(amino_acid_keys)], amino_acid_keys)

# Assign FillColor ensuring alignment
heatmap_data$FillColor <- ifelse(
  heatmap_data$IsDifferent,
  aminoAcidColors[heatmap_data$AminoAcid],
  "lightgrey"
)

# Add reference sequence as a top row
reference_data <- data.frame(
  Clade = "Consensus",  # Change the name here
  Position = 1:ncol(sequence_data),
  AminoAcid = reference_sequence,
  IsDifferent = FALSE,  # Reference sequence has no differences
  FillColor = aminoAcidColors[reference_sequence]
)

# Add reference sequence as a top row and combine
heatmap_data_with_reference <- rbind(reference_data, heatmap_data)

# Reorder the levels of the "Clade" factor to ensure "Consensus" is at the top
heatmap_data_with_reference$Clade <- factor(heatmap_data_with_reference$Clade, levels = c("Consensus", unique(Clade_names)))

# Plot the mutation heatmap
# Plot the mutation heatmap
heatmap_plot_with_reference <- ggplot(heatmap_data_with_reference, aes(x = Position, y = Clade, fill = FillColor)) +
  geom_tile(color = "white") +
  scale_fill_identity(
    name = "Amino Acid",
    guide = guide_legend(
      nrow = 11,  # Adjust the number of rows as needed
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.4
    ),
    labels = names(aminoAcidColors),
    breaks = aminoAcidColors
  ) +
  theme_classic() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.text = element_text(size = 8)
  ) +
  labs(x = "Amino acid position", y = "Clade")

# Print the updated heatmap with the reference row
print(heatmap_plot_with_reference)