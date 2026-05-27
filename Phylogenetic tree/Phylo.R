# Load necessary libraries
library(ggtree)       # For tree visualization
library(treeio)       # For tree file handling
library(ape)          # To read Newick tree
library(readxl)       # For reading Excel files
library(ggplot2)      # For visualization
library(dplyr)        # For data manipulation
library(ggtreeExtra)  # For concentric heatmaps around the tree
library(janitor)      # For cleaning column names

# Step 1: Load the Tree File
tree_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/Phylogenetic tree/NP_H1_combined.tre"  
tree <- read.tree(tree_file)  # Load Newick format tree

# Step 2: Load Metadata and Clean Names
metadata_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/Phylogenetic tree/NP_H1_phylogen.xlsx"
metadata <- read_excel(metadata_file) %>%
  clean_names() %>%  # Standardize column names
  rename(
    label = strain_name,  # Correct the name for strain
    Year_of_Isolation = collection_date,
    Clade = global_swine_h1_clade
  ) %>%
  distinct(label, .keep_all = TRUE) %>%
  mutate(
    Clade = as.factor(Clade),  # Ensure Clade is treated as categorical
    Year_of_Isolation = as.numeric(Year_of_Isolation)  # Convert Year_of_Isolation to numeric
  )

# Step 3: Merge Metadata with Tree Tips
tree_data <- full_join(metadata, data.frame(label = tree$tip.label), by = "label")

# Step 4: Prepare Bootstrap Data for Internal Nodes
node_data <- data.frame(
  node = (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)),  # Internal nodes
  bootstrap = as.numeric(tree$node.label)            # Bootstrap values (if available)
)

# Step 5: Plot the Circular Phylogenetic Tree
p <- ggtree(tree, layout = "circular") %<+% tree_data

# Step 5.1: Add Concentric Heatmaps
p <- p + 
  geom_fruit(
    geom = geom_tile,
    aes(fill = Year_of_Isolation),
    offset = 0.03,
    width = 0.05
  ) +
  scale_fill_viridis_c(name = "Year of Isolation", na.value = "grey") +
  geom_fruit(
    geom = geom_tile,
    aes(fill = Clade),
    offset = 0.07,
    width = 0.05
  ) +
  scale_fill_brewer(palette = "Set2", name = "Clade")

# Step 5.2: Add Bootstrap Values at Internal Nodes (if available)
if(!all(is.na(node_data$bootstrap))) {
  p <- p + 
    geom_point(
      data = node_data,  # Use fortified tree for positions
      aes(x = x, y = y, size = bootstrap),
      shape = 21, color = "black", fill = "white", stroke = 0.8, alpha = 0.6
    ) +
    scale_size_continuous(name = "Bootstrap Values", range = c(2, 6))
}

# Step 5.3: Add Strain Names as Tip Labels
p <- p + 
  geom_tiplab(aes(label = label), size = 2.5, align = TRUE, offset = 0.1)

# Step 5.4: Finalize the Theme and Add Title
p <- p + 
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Circular Phylogenetic Tree with Concentric Heatmaps and Bootstrap Values")

# Step 6: Save the Final Tree Plot
output_file <- "/Users/raj/Desktop/Codes_bioinfro/R scripts/Phylogenetic tree/NP_H1_combined_tree_with_metadata.png"
ggsave(output_file, p, width = 14, height = 14)

# Print the tree plot
print(p)
