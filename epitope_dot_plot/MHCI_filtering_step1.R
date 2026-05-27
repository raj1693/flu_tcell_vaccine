library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(scales)

# Read the Excel data
data <- read_excel("/Users/raj/Desktop/Codes_bioinfro/R scripts/epitope_dot_plot/NP_M1_MHC_1.xlsx")

# Assuming your data is in a data frame named 'data'

# Group by the desired columns and count occurrences
data_counts <- data %>%
  group_by(`MHC allele`, `Amino acid position`, facet_column) %>%
  summarize(count = n())

# Filter for rows where count is 3
filtered_data <- data_counts %>%
  filter(count == 3)

# Reorder the facet column
filtered_data$facet_column <- factor(filtered_data$facet_column, levels = c("NP H1N1", "M1 H1N1", "NP H3N2", "M1 H3N2"))

# Determine the maximum 'Amino acid position' for dynamic scaling
max_position <- max(filtered_data$`Amino acid position`)

ggplot(filtered_data, aes(x = `Amino acid position`, y = `MHC allele`)) +
  geom_point() +
  facet_wrap(~ facet_column, ncol = 2) +
  labs(x = "Amino acid position", y = "MHC allele") +
  theme_minimal() 