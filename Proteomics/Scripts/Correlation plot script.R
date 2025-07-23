##### updated correlation with custom legend to show aggregation potential in each varient using pearon coeff


########################################################################################
##### Trial — Updated Correlation Plot with Centered Axis and Absolute Outlier Comparison
########################################################################################

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(cowplot)

# Define variants and labels
variants <- c("R113", "G184", "Pbabe", "Dual")
variant_labels <- c("R113" = "R113W", "G184" = "G184R", "Pbabe" = "pBABE", "Dual" = "R113W/G184R")

# Define colors
fill_colors <- c(
  "R113_MRPL" = "#8B5A00",  "R113_MRPS" = "#FFA500",
  "G184_MRPL" = "#008B8B",  "G184_MRPS" = "#8EE5EE",
  "Pbabe_MRPL" = "#104E8B", "Pbabe_MRPS" = "#1E90FF",
  "Dual_MRPL"  = "#68228B", "Dual_MRPS"  = "#DA70D6"
)

# Legend labels
legend_labels <- c(
  "R113_MRPL" = "R113W (Large)",  "R113_MRPS" = "R113W (Small)",
  "G184_MRPL" = "G184R (Large)",  "G184_MRPS" = "G184R (Small)",
  "Pbabe_MRPL" = "pBABE (Large)", "Pbabe_MRPS" = "pBABE (Small)",
  "Dual_MRPL"  = "R113W/G184R (Large)", "Dual_MRPS" = "R113W/G184R (Small)"
)

# Classify MRPL/MRPS
classify_mito_ribo <- function(gene) {
  if (grepl("^MRPL", gene)) return("MRPL")
  if (grepl("^MRPS", gene)) return("MRPS")
  return(NA)
}

# Load variant-specific data
load_variant_data <- function(variant, fraction) {
  file_path <- paste0("Merged_Final_", fraction, "_", variant, "vs", fraction, "_WT.txt")
  df <- read.delim(file_path)
  df$Variant <- variant
  df$RiboType <- sapply(df$Genes, classify_mito_ribo)
  df <- df %>% filter(!is.na(RiboType))
  df$Group <- paste0(variant, "_", df$RiboType)
  return(df)
}

# Merge soluble and insoluble data
merge_sol_insol <- function(variant) {
  insol <- load_variant_data(variant, "Insol") %>%
    select(Genes, logFC_insol = logFC, RiboType, Group)
  sol <- load_variant_data(variant, "Sol") %>%
    select(Genes, logFC_sol = logFC)
  merged <- inner_join(insol, sol, by = "Genes")
  merged$Variant <- variant
  merged$Delta <- abs(merged$logFC_insol - merged$logFC_sol)
  return(merged)
}

# Combine all variants
merged_all <- map_dfr(variants, merge_sol_insol)

# Absolute max for centered axis
max_abs_fc <- max(abs(c(merged_all$logFC_sol, merged_all$logFC_insol)), na.rm = TRUE)
fc_range <- c(-max_abs_fc, max_abs_fc)

# Identify top outliers per variant
label_outliers <- merged_all %>%
  group_by(Variant) %>%
  arrange(desc(Delta)) %>%
  slice_head(n = 3) %>%
  ungroup()

# Correlation plot function
make_corr_plot <- function(df, label_df, variant_code, axis_limits) {
  corr_coef <- round(cor(abs(df$logFC_sol), abs(df$logFC_insol), use = "complete.obs"), 2)
  variant_title <- variant_labels[variant_code]
  
  ggplot(df, aes(x = logFC_sol, y = logFC_insol, fill = Group)) +
    geom_point(shape = 21, size = 3, color = "black", alpha = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_text_repel(data = label_df, aes(label = Genes), size = 3.5, max.overlaps = 100) +
    scale_fill_manual(values = fill_colors, labels = legend_labels) +
    scale_x_continuous(limits = axis_limits) +
    scale_y_continuous(limits = axis_limits) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Soluble vs Insoluble –", variant_title),
      subtitle = paste("Pearson r =", corr_coef),
      x = "log2FC (Soluble)",
      y = "log2FC (Insoluble)"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      panel.grid.minor = element_blank()
    )
}

# Generate plots for each variant
plots <- lapply(variants, function(v) {
  df <- merged_all %>% filter(Variant == v)
  labels <- label_outliers %>% filter(Variant == v)
  make_corr_plot(df, labels, v, axis_limits = fc_range)
})
###
# Extract clean legend with bigger text
legend_plot <- ggplot(
  data.frame(Group = names(fill_colors), x = 1, y = 1),
  aes(x = x, y = y, fill = Group)
) +
  geom_point(shape = 21, size = 4, color = "black") +
  scale_fill_manual(
    values = fill_colors,
    labels = legend_labels,
    name   = "Variant + Subunit",
    guide  = guide_legend(
      title.position  = "top",          # title above keys
      title.hjust      = 0.5,           # centered
      label.position  = "right",
      label.hjust     = 0,
      ncol            = 1,
      byrow           = TRUE,
      keywidth        = unit(1.2, "lines"),
      keyheight       = unit(1.2, "lines"),
      title.theme     = element_text(size = 16, face = "bold"),  # legend title
      label.theme     = element_text(size = 14)                  # legend labels
    )
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin     = margin(5, 5, 5, 5),
    # fallback in case guide_legend() doesn’t catch it
    legend.title      = element_text(size = 18, face = "bold"),
    legend.text       = element_text(size = 16),
    legend.key.size   = unit(1.2, "lines")
  )

legend <- get_legend(legend_plot)

# Combine with the main grid as before
main_grid  <- plot_grid(plotlist = plots, ncol = 2)
final_plot <- plot_grid(main_grid, legend, ncol = 2, rel_widths = c(4.5, 1.2))

ggsave(
  "correlation_plot_with_custom_labels_and_legendbigtext.pdf",
  final_plot, width = 14, height = 10
)

