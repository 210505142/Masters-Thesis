######################### individual varients comparing top dep in pBABE

#test 1
library(tidyverse)
library(ggrepel)
library(ggpubr)


library(tidyverse)
library(ggrepel)
library(ggpubr)

# Define variants and labels
variants <- c("R113", "G184", "Pbabe", "Dual")
variant_labels <- c("R113" = "R113W", "G184" = "G184R", "Pbabe" = "pBABE", "Dual" = "R113W/G184R")

# Helper: assign colors per variant + subunit + fraction
fill_colors <- c(
  "R113_MRPL_Sol" = "#A65628", "R113_MRPL_Insol" = "#8B5A00",
  "R113_MRPS_Sol" = "#FDB863", "R113_MRPS_Insol" = "#FFA500",
  
  "G184_MRPL_Sol" = "#66C2A5", "G184_MRPL_Insol" = "#008B8B",
  "G184_MRPS_Sol" = "#B2DFDB", "G184_MRPS_Insol" = "#8EE5EE",
  
  "Pbabe_MRPL_Sol" = "#6495ED", "Pbabe_MRPL_Insol" = "#104E8B",
  "Pbabe_MRPS_Sol" = "#ADD8E6", "Pbabe_MRPS_Insol" = "#1E90FF",
  
  "Dual_MRPL_Sol" = "#9370DB", "Dual_MRPL_Insol" = "#68228B",
  "Dual_MRPS_Sol" = "#D8BFD8", "Dual_MRPS_Insol" = "#DA70D6"
)

# Classify ribosome type
classify_mito_ribo <- function(gene) {
  if (grepl("^MRPL", gene)) return("MRPL")
  if (grepl("^MRPS", gene)) return("MRPS")
  return(NA)
}

# Load data for one variant/fraction
load_variant_data <- function(variant, fraction) {
  file_path <- paste0("Merged_Final_", fraction, "_", variant, "vs", fraction, "_WT.txt")
  df <- read.delim(file_path)
  df$Variant  <- variant
  df$RiboType <- sapply(df$Genes, classify_mito_ribo)
  df <- df %>% filter(!is.na(RiboType))
  df$Group <- paste0(variant, "_", df$RiboType, "_", fraction)
  df$Fraction <- fraction
  return(df)
}

# Load all data
insol_all <- map_dfr(variants, ~ load_variant_data(.x, "Insol"))
sol_all   <- map_dfr(variants, ~ load_variant_data(.x, "Sol"))

# Determine axis limits
all_data <- bind_rows(insol_all, sol_all)
x_limits <- range(all_data$logFC, na.rm = TRUE)
y_limits <- c(0, max(-log10(all_data$P.Value), na.rm = TRUE))

# Identify top 10 MRPLs from pBABE
top_pbabe_mrpl <- bind_rows(
  insol_all %>% filter(Variant == "Pbabe", RiboType == "MRPL"),
  sol_all   %>% filter(Variant == "Pbabe", RiboType == "MRPL")
) %>%
  arrange(P.Value) %>%
  distinct(Genes, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Genes)

# Filter to those top genes across all data
insol_top <- insol_all %>% filter(Genes %in% top_pbabe_mrpl)
sol_top   <- sol_all   %>% filter(Genes %in% top_pbabe_mrpl)

# Updated plot function per variant


#######
make_mrpl_comparison_plot <- function(variant) {
  df <- bind_rows(
    insol_top   %>% filter(Variant == variant),
    sol_top     %>% filter(Variant == variant)
  )
  
  ggplot(df, aes(x = logFC, y = -log10(P.Value), fill = Group)) +
    geom_point(shape = 21, size = 3.5, color = "black", alpha = 0.9) +
    geom_text_repel(aes(label = Genes), size = 3.3, max.overlaps = 100) +
    
    # <-- BIGGER LEGEND HERE -->
    scale_fill_manual(
      values = fill_colors,
      name   = "Subunit + Fraction",
      guide  = guide_legend(
        title.position = "top",
        title.hjust    = 0.5,
        label.position = "right",
        direction      = "vertical",
        ncol           = 1,
        keywidth       = unit(1.2, "lines"),
        keyheight      = unit(1.2, "lines"),
        title.theme    = element_text(size = 16, face = "bold"),
        label.theme    = element_text(size = 14)
      )
    ) +
    
    facet_wrap(~ Fraction) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Top 10 MRPLs from pBABE applied to", variant_labels[variant]),
      x     = "log2 Fold Change",
      y     = "-log10(P-value)"
    ) +
    theme(
      legend.position   = "right",
      legend.background = element_rect(fill = "white", color = NA),
      # backup settings in case guide_legend() doesn't catch everything:
      legend.title      = element_text(size = 18, face = "bold"),
      legend.text       = element_text(size = 16),
      legend.key.size   = unit(1.2, "lines"),
      panel.border      = element_rect(color = "black", fill = NA)
    )
}
# Generate and combine all plots
variant_plots <- lapply(variants, make_mrpl_comparison_plot)

# Arrange all into one PDF
combined_plot <- ggarrange(
  plotlist = variant_plots,
  ncol = 2, nrow = 2,
  labels = NULL
)

# Save
ggsave("Top10_pBABE_MRPL_colored_by_fraction_across_variantstest.pdf", combined_plot, width = 16, height = 11)


#############################
##### Sol insol split

# Modified plot function for specific fraction
make_mrpl_fraction_plot <- function(variant, fraction) {
  df <- if (fraction == "Sol") {
    sol_top %>% filter(Variant == variant)
  } else {
    insol_top %>% filter(Variant == variant)
  }
  
  ggplot(df, aes(x = logFC, y = -log10(P.Value), fill = Group)) +
    geom_point(shape = 21, size = 3.5, color = "black", alpha = 0.9) +
    geom_text_repel(aes(label = Genes), size = 4, max.overlaps = 100) +
    scale_fill_manual(values = fill_colors) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme_minimal(base_size = 16) +
    labs(
      title = paste("Top 10 MRPLs from pBABE –", variant_labels[variant], fraction, "Fraction"),
      x = "log2 Fold Change",
      y = "-log10(P-value)"
    ) +
    theme(
      plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x    = element_text(size = 18),
      axis.title.y    = element_text(size = 18),
      axis.text.x     = element_text(size = 16),
      axis.text.y     = element_text(size = 16),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.title    = element_text(size = 18, face = "bold"),
      legend.text     = element_text(size = 16),
      legend.key.size = unit(1.2, "lines"),
      panel.border    = element_rect(color = "black", fill = NA)
    )
}


# Generate plots
sol_plots   <- lapply(variants, make_mrpl_fraction_plot, fraction = "Sol")
insol_plots <- lapply(variants, make_mrpl_fraction_plot, fraction = "Insol")

# Combine and save separately
sol_combined <- ggarrange(plotlist = sol_plots, ncol = 2, nrow = 2)
insol_combined <- ggarrange(plotlist = insol_plots, ncol = 2, nrow = 2)

# Save PDFs
ggsave("Top10_pBABE_MRPL_Soluble_only.pdf", sol_combined, width = 16, height = 11)
ggsave("Top10_pBABE_MRPL_Insoluble_only.pdf", insol_combined, width = 16, height = 11)
############################
################### axis title bigger
