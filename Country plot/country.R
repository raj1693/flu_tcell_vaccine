# Install necessary packages if missing
if (!require("readxl")) install.packages("readxl", dependencies=TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies=TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies=TRUE)
if (!require("tidyr")) install.packages("tidyr", dependencies=TRUE)
if (!require("sf")) install.packages("sf", dependencies=TRUE)
if (!require("rnaturalearth")) install.packages("rnaturalearth", dependencies=TRUE)
if (!require("rnaturalearthdata")) install.packages("rnaturalearthdata", dependencies=TRUE)

# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Step 1: Load your Excel file
# Replace the path with your actual file location
data <- read_excel("/Users/raj/Desktop/Codes_bioinfro/R scripts/Country plot/NP_M1_H1_H3_plot_cleaned.xlsx")

# Step 2: Remove rows with missing data
data <- data[complete.cases(data), ]

# Step 3: Standardize country names for better matching
data <- data %>%
  mutate(Country = tolower(Country)) %>% # Convert all country names to lowercase
  mutate(Country = case_when(
    Country == "czech republic" ~ "czechia",  # Standardize to "czechia"
    Country == "the unites states of america" ~ "united states of america", # Fix typo in "unites"
    TRUE ~ Country
  ))

# Preview the cleaned data
print(head(data))

# Step 4: Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Standardize country names in the world map for better matching
world <- world %>%
  mutate(name = tolower(name))

# Step 5: Reshape your data to long format for faceted plotting
data_long <- data %>%
  pivot_longer(cols = -Country, names_to = "Metric", values_to = "Value")

# Step 6: Merge map data with your dataset
map_data <- world %>%
  left_join(data_long, by = c("name" = "Country"))

# Filter out rows with missing values before plotting
map_data_filtered <- map_data[!is.na(map_data$Value), ]

# Step 7: Create the heatmap
heatmap <- ggplot(map_data_filtered, aes(fill = Value)) +
  geom_sf() + 
  scale_fill_gradientn(
    colours = c("lightyellow", "orange", "red"),
    name = "No. of Seqs",
    trans = "log10",  # Log transformation for better visualization
    na.value = "transparent"
  ) + 
  facet_wrap(~ Metric, ncol = 2) + 
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each facet
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

# Step 8: Display the heatmap
print(heatmap)
