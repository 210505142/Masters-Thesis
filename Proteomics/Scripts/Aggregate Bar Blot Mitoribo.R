# ----------------------------------------------------------------------------------
# ADDITIONAL AGGREGATION SCORE PLOT
# ----------------------------------------------------------------------------------
##### Add absolute value for inputs in agrregate function
# Add aggregation score column

###################################################
#### updated to use absoulte value
# Modified aggregation score calculation with absolute value
merged_all <- merged_all %>%
  mutate(Aggregation_Score = abs(logFC_insol - logFC_sol))  # Added abs()

# Identify top 5 aggregation-prone proteins per variant (now using absolute values)
top_aggregators <- merged_all %>%
  group_by(Variant) %>%
  arrange(desc(Aggregation_Score)) %>%  # Now sorting absolute values
  slice_head(n = 5) %>%
  ungroup()

# Updated aggregation score plot
agg_plot <- ggplot(merged_all, aes(x = reorder(Genes, Aggregation_Score), 
                                   y = Aggregation_Score, 
                                   fill = Group)) +
  geom_col(show.legend = FALSE) +
  geom_text_repel(
    data = top_aggregators,
    aes(label = Genes),
    size = 3,
    max.overlaps = 50,
    direction = "y",
    nudge_y = 0.5
  ) +
  facet_wrap(~ Variant, scales = "free_x", labeller = labeller(Variant = variant_labels)) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Absolute Aggregation Score by Variant",
    subtitle = "Score = |log2FC(insoluble) - log2FC(soluble)|",  # Updated notation
    x = "Proteins (ranked by absolute aggregation score)",
    y = "Absolute Aggregation Score"  # Updated axis label
  )
#######################
###################
## Full pipeline test with legend sucess
####

# Load required libraries
library(tidyverse)
library(ggrepel)

# Define variants and labels
variants <- c("R113", "G184", "Pbabe", "Dual")
variant_labels <- c("R113" = "R113W", "G184" = "G184R", "Pbabe" = "pBABE", "Dual" = "R113W/G184R")

# Custom fill colors for variant + subunit
fill_colors <- c(
  "R113_MRPL" = "#8B5A00",  "R113_MRPS" = "#FFA500",
  "G184_MRPL" = "#008B8B",  "G184_MRPS" = "#8EE5EE",
  "Pbabe_MRPL" = "#104E8B", "Pbabe_MRPS" = "#1E90FF",
  "Dual_MRPL"  = "#68228B", "Dual_MRPS"  = "#DA70D6"
)

# Legend display labels
legend_labels <- c(
  "R113_MRPL" = "R113W (Large)",  "R113_MRPS" = "R113W (Small)",
  "G184_MRPL" = "G184R (Large)",  "G184_MRPS" = "G184R (Small)",
  "Pbabe_MRPL" = "pBABE (Large)", "Pbabe_MRPS" = "pBABE (Small)",
  "Dual_MRPL"  = "R113W/G184R (Large)", "Dual_MRPS" = "R113W/G184R (Small)"
)

# Function to classify mitochondrial ribosome type
classify_mito_ribo <- function(gene) {
  if (grepl("^MRPL", gene)) return("MRPL")
  if (grepl("^MRPS", gene)) return("MRPS")
  return(NA)
}

# Load data for a given variant and fraction
load_variant_data <- function(variant, fraction) {
  file_path <- paste0("Merged_Final_", fraction, "_", variant, "vs", fraction, "_WT.txt")
  df <- read.delim(file_path)
  df$Variant  <- variant
  df$RiboType <- sapply(df$Genes, classify_mito_ribo)
  df <- df %>% filter(!is.na(RiboType))
  df$Group <- paste0(variant, "_", df$RiboType)
  return(df)
}

# Load and combine soluble and insoluble data across variants
merged_all <- map_dfr(variants, function(v) {
  insol <- load_variant_data(v, "Insol") %>%
    select(Genes, logFC_insol = logFC, RiboType, Group)
  sol <- load_variant_data(v, "Sol") %>%
    select(Genes, logFC_sol = logFC)
  combined <- inner_join(insol, sol, by = "Genes")
  combined$Variant <- v
  combined
})

# Calculate absolute aggregation score
merged_all <- merged_all %>%
  mutate(Aggregation_Score = abs(logFC_insol - logFC_sol))

# Identify top 5 aggregators per variant
top_aggregators <- merged_all %>%
  group_by(Variant) %>%
  arrange(desc(Aggregation_Score)) %>%
  slice_head(n = 5) %>%
  ungroup()

# Create aggregation bar plot
agg_plot <- ggplot(merged_all, aes(x = reorder(Genes, Aggregation_Score), 
                                   y = Aggregation_Score, 
                                   fill = Group)) +
  geom_col(show.legend = TRUE) +
  geom_text_repel(
    data = top_aggregators,
    aes(label = Genes),
    size = 3,
    max.overlaps = 50,
    direction = "y",
    nudge_y = 0.5
  ) +
  facet_wrap(~ Variant, scales = "free_x", labeller = labeller(Variant = variant_labels)) +
  scale_fill_manual(values = fill_colors, labels = legend_labels, name = "Variant + Subunit") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Absolute Aggregation Score by Variant",
    subtitle = "Score = |log2FC(insoluble) - log2FC(soluble)|",
    x = "Proteins (ranked by absolute aggregation score)",
    y = "Absolute Aggregation Score"
  )

# Save to PDF
ggsave("aggregation_score_plot.pdf", agg_plot, width = 14, height = 8)
