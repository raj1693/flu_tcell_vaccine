# Load required libraries
library(ggplot2)
library(dplyr)
library(ggtree)
library(ape)
library(cowplot)

# Path to your Newick file
nwk_file <- ""  # Replace with your actual path

# Read the phylogenetic tree
phylo_tree <- read.tree(nwk_file)

# Create a ggtree plot
tree_plot <- ggtree(phylo_tree, layout = "rectangular") +
  theme_tree2() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Align clade names with heatmap rows
tree_plot <- tree_plot %<+% data.frame(Clade = Clade_names)

# Plot the mutation heatmap
heatmap_plot <- ggplot(heatmap_data_with_reference, aes(x = Position, y = Clade, fill = FillColor)) +
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
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_blank(),  # Hide y-axis text for the heatmap
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
  labs(x = "Position", y = "")

# Combine the tree and heatmap using cowplot
combined_plot <- plot_grid(
  tree_plot, heatmap_plot,
  align = "h", nrow = 1, rel_widths = c(1, 2)
)

# Print the combined plot
print(combined_plot)
