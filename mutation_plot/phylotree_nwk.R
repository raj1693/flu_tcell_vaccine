# Load required libraries
library(ape)
library(phangorn)
library(Biostrings)

# Step 1: Read the protein sequences from a FASTA file
fasta_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/consensus_sequences_aligned.fasta  # Replace with your actual file path

# Load the protein sequences into an AAStringSet object
sequences <- readAAStringSet(fasta_file)

# Step 2: Calculate the distance matrix for protein sequences
dist_matrix <- dist.ml(alignment)

# Step 3: Calculate a starting tree using the neighbor-joining method
# Create a distance matrix
dist_matrix <- dist.ml(alignment)
nj_tree <- nj(dist_matrix)

# Step 4: Fit a maximum likelihood tree
# Define a model of evolution (e.g., JTT)
fit <- pml(nj_tree, data = alignment)

# Optimize the tree using maximum likelihood
fit_ml <- optim.pml(fit, model = "JTT")

# Step 5: Write the ML tree to a Newick file
output_file <- "/path/to/output/my_ml_tree.nwk"  # Replace with your desired output path
write.tree(fit_ml$tree, file = output_file)

# Confirm the process
cat("Maximum likelihood tree written to:", output_file, "\n")


