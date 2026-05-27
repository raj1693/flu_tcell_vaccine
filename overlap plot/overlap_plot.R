library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggpubr)

# Load data from Excel
file_path <- "/Users/raj/Desktop/Codes_bioinfo/R scripts/overlap plot/Overlap plot.xlsx"
sheets <- excel_sheets(file_path)
data_list <- lapply(sheets, function(sheet) {
  read_excel(file_path, sheet = sheet) %>%
    dplyr::mutate(Start = as.numeric(Start) , 
                  End = as.numeric(End),
                  T_Cell = case_when(MHC == "II" ~ "CD4" , 
                                       MHC == "I" ~ "CD8",
                                       TRUE ~ MHC),
                  Host_Dummy = case_when(Host == "HUMAN" ~ "Human",
                                         Host == "MICE" ~ "Mice" , 
                                         TRUE ~ Host),
                  Group = paste0(T_Cell, " " ,Host_Dummy)) %>%
    dplyr::select(Allele , Start , End , Host , MHC, Group) %>%
    group_by(Allele) %>%
    mutate(Allele_Modified = str_pad(Allele, width = nchar(Allele) + row_number() - 1, side = "left")) %>%
    ungroup() %>%
    dplyr::mutate(Group = factor(Group , levels = c("CD8 Human", "CD4 Human", "CD8 Mice", "CD4 Mice"))) %>%
    dplyr::arrange(Start) %>%
    dplyr::mutate(Allele_Modified = factor(Allele_Modified , levels = .$Allele_Modified))
    }
  )

names(data_list) <- sheets

# Define colors for groups
color_map <- c("CD8 Human" = "red", "CD4 Human" = "orange", "CD8 Mice" = "blue", "CD4 Mice" = "yellow")

# Create individual plots as floating bar charts with more transparency for mice epitopes
plot_list <- lapply(sheets, function(sheet) {
  df <- data_list[[sheet]]
  
  # Set lower alpha for mice groups
  df$alpha_val <- case_when(
    df$Group == "CD8 Mice" ~ 0.3,
    df$Group == "CD4 Mice" ~ 0.3,
    TRUE ~ 1.0
  )
  
  df = df %>% 
    dplyr::mutate(edge_size = case_when(Host == "HUMAN" ~ 1,
                                        Host == "MICE" ~ 0.1 , 
                                        TRUE ~ 0.1))
  
  ggplot(df, aes(xmin = Start, xmax = End, 
                 y = factor(Allele_Modified, levels = unique(Allele_Modified)), 
                 fill = Group)) +
    geom_rect(aes(ymin = as.numeric(factor(Allele_Modified, levels = unique(Allele_Modified))) - 0.4, 
                  ymax = as.numeric(factor(Allele_Modified, levels = unique(Allele_Modified))) + 0.4,
                  alpha = alpha_val,
                  size = edge_size), 
              color = "black") +
    scale_fill_manual(values = color_map) +
    scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    scale_size_continuous(range = c(0.1 , 0.8), guide = "none") +
    labs(title = sheet, x = "Amino Acid Position", y = "Allele") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA), 
          panel.grid.major = element_line(color = "lightgrey"), 
          panel.grid.minor = element_line(color = "lightgrey"))
})


# Arrange all plots into a 2x2 grid
final_plot <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

# Save the final plot
ggsave("overlap_plot.png", final_plot, width = 12, height = 8)

# Print the final plot
print(final_plot)

