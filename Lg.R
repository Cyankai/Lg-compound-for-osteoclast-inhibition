# ==============================================================================
# Script Name: analysis_Lg_osteoclast.R
# Description: Integrated Analysis of RNA-seq and H3K9me3 CUT&Tag Profiling 
#              Investigating the Inhibitory Mechanism of Lg on Osteoclastogenesis
# Input Data:  'res_df' (RNAseq_DESeq2_results), 'heat_matrix' (RNAseq_normalized
#              _counts_matrix), 'Samples_peakAnno_matrix_data' (CUT&Tag Peak Anno
#              tation), and Samples_analyze (Fitted DBA Object from DiffBind).
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(DiffBind)
})

setwd("/Users/kenny/Desktop/NRN/Cuttag_results")
 
1. Double-Track Dotplot: Osteoclast Markers Expression vs Promoter H3K9me3 

# Canonical osteoclast marker set
oc_markers <- c(
  "Nfatc1", "Fos", "Jun",               # Core differentiation transcription factors
  "Acp5", "Ctsk", "Mmp9", "Car2",       # Structural enzymes & resorption markers
  "Dcstamp", "Ocstamp", "Itgb3",        # Fusion & adhesion molecules
  "Tnfrsf11a", "Oscar"                  # Lineage receptors
)

# Extract promoter CUT&Tag signals (representative peak with max fold change)
peak_oc <- Samples_peakAnno_matrix_data_promoter %>%
  filter(SYMBOL %in% oc_markers) %>%
  arrange(desc(abs(Fold)), p.value, abs(distanceToTSS)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(
    gene = SYMBOL,
    chip_log2FC = Fold,
    chip_pval = p.value
  )

# Extract RNA-seq expression
rna_oc <- res_df %>%
  filter(gene %in% oc_markers) %>%
  dplyr::select(
    gene,
    rna_log2FC = log2FoldChange,
    rna_pval = pvalue,
    rna_padj = padj
  )

# Reshape into long-format dataframe for paired visualization
dot_df <- full_join(rna_oc, peak_oc, by = "gene") %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = c("Assay", ".value"),
    names_pattern = "(rna|chip)_(.*)"
  ) %>%
  mutate(
    Assay = ifelse(Assay == "rna", "RNA-seq (Expression)", "H3K9me3 (Promoter Peak)"),
    neg_log10_pval = -log10(pmax(pval, 1e-10)),
    is_sig = ifelse(pval < 0.05, "p < 0.05", "NS (Not Significant)")
  )

# Order genes by descending RNA downregulation
gene_order <- rna_oc %>% arrange(rna_log2FC) %>% pull(gene)
dot_df$gene <- factor(dot_df$gene, levels = gene_order)

# Export plot data table
write.csv(dot_df, "plot_data_tables/SourceData_Dotplot_OC_Markers.csv", row.names = FALSE)

# Render Dotplot
p_dot <- ggplot(dot_df, aes(x = Assay, y = gene)) +
  geom_hline(yintercept = seq_along(gene_order), color = "grey92", linewidth = 0.5) +
  geom_point(
    aes(size = neg_log10_pval, fill = log2FC, color = is_sig),
    shape = 21, stroke = 1
  ) +
  scale_color_manual(values = c("p < 0.05" = "black", "NS (Not Significant)" = "grey70")) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
    name = expression(Log[2]~Fold~Change)
  ) +
  scale_size_continuous(
    range = c(2.5, 8.5), breaks = c(1, 2, 4, 8),
    name = expression(-Log[10]~(italic(P)-value))
  ) +
  guides(
    color = guide_legend(title = "Significance", override.aes = list(size = 4, fill = "white")),
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Osteoclast Markers: RNA Repression without H3K9me3 Gain",
    subtitle = "H3K9me3 peaks exhibit no repressive enrichment at promoters of downregulated osteoclast genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold.italic", size = 11, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    legend.position = "right"
  )

ggsave("Fig_Dotplot_Osteoclast_Markers.pdf", plot = p_dot, width = 6.5, height = 6)
ggsave("Fig_Dotplot_Osteoclast_Markers.png", plot = p_dot, width = 6.5, height = 6, dpi = 300)

 
# 2. DiffBind CUT&Tag Volcano Plot 

# Extract complete peak table from DBA analysis
report_peaks <- dba.report(Samples_analyze, contrast = 1, th = 1, bCounts = FALSE) %>%
  as.data.frame()

p_threshold <- 0.05
fold_threshold <- 0

report_peaks <- report_peaks %>%
  mutate(
    neg_log10_p = -log10(pmax(p.value, 1e-15)),
    Status = case_when(
      p.value < p_threshold & Fold <= -fold_threshold ~ "Down",
      p.value < p_threshold & Fold >=  fold_threshold ~ "Up",
      TRUE                                            ~ "NS"
    )
  )

n_down <- sum(report_peaks$Status == "Down", na.rm = TRUE)
n_up   <- sum(report_peaks$Status == "Up", na.rm = TRUE)
n_ns   <- sum(report_peaks$Status == "NS", na.rm = TRUE)

# Export plot data table
write.csv(report_peaks, "plot_data_tables/SourceData_Volcano_H3K9me3_Peaks.csv", row.names = FALSE)

# Render symmetric Volcano Plot
max_abs_fold <- max(abs(report_peaks$Fold), na.rm = TRUE)
x_limit <- ceiling(max_abs_fold * 1.05)
y_limit <- max(report_peaks$neg_log10_p, na.rm = TRUE) * 1.05

p_volcano <- ggplot(report_peaks, aes(x = Fold, y = neg_log10_p)) +
  geom_vline(xintercept = c(-fold_threshold, fold_threshold), linetype = "dashed", color = "grey65", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey65", linewidth = 0.5) +
  geom_point(data = filter(report_peaks, Status == "NS"),
             color = "#D1D5DB", alpha = 0.4, size = 1.2) +
  geom_point(data = filter(report_peaks, Status == "Down"),
             color = "#3B82F6", alpha = 0.8, size = 1.8) +
  geom_point(data = filter(report_peaks, Status == "Up"),
             color = "#EF4444", alpha = 0.8, size = 1.8) +
  annotate(
    "text", x = -x_limit * 0.92, y = y_limit * 0.95,
    label = paste0("Downregulated Peaks\n(P < 0.05)\nn = ", n_down),
    color = "#1D4ED8", fontface = "bold", size = 4.2, hjust = 0, vjust = 1
  ) +
  annotate(
    "text", x = x_limit * 0.92, y = y_limit * 0.95,
    label = paste0("Upregulated Peaks\n(P < 0.05)\nn = ", n_up),
    color = "#B91C1C", fontface = "bold", size = 4.2, hjust = 1, vjust = 1
  ) +
  scale_x_continuous(limits = c(-x_limit, x_limit), breaks = seq(-x_limit, x_limit, by = 1)) +
  scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = expression(Log[2]~Fold~Change~(Lg~vs~Control)),
    y = expression(-Log[10]~(italic(P)-value)),
    title = "Differential H3K9me3 CUT&Tag Peaks",
    subtitle = paste0("Total consensus peaks: ", nrow(report_peaks))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave("Fig_Volcano_H3K9me3_Clean.pdf", plot = p_volcano, width = 6, height = 5)
ggsave("Fig_Volcano_H3K9me3_Clean.png", plot = p_volcano, width = 6, height = 5, dpi = 300)

 
# 3. Two-Module Osteoclast Repression Heatmap  
neg_regulators <- c(
  "Irf8", "Spen", "Cdkn1a", "Btg3", "Egr1", "Mafb", "Bcl6", 
  "Nrip1", "Rgs18", "Foxo1", "Foxo3", "Foxo4", "Id2"
)

oc_markers_subset <- c(
  "Ctsk", "Acp5", "Nfatc1", "Calcr", "Car2", "Mmp9", "Dcstamp", 
  "Ocstamp", "Atp6v0d2"
)

target_genes_df <- data.frame(
  gene = c(neg_regulators, oc_markers_subset),
  Function = c(
    rep("Negative Regulators / Checkpoints", length(neg_regulators)),
    rep("Osteoclast Markers & Resorption", length(oc_markers_subset))
  )
)

valid_genes <- intersect(target_genes_df$gene, rownames(heat_matrix))
target_genes_df <- target_genes_df %>% filter(gene %in% valid_genes)

# Cluster independently within functional tiers
mat_neg <- heat_matrix[intersect(neg_regulators, valid_genes), ]
mat_oc  <- heat_matrix[intersect(oc_markers_subset, valid_genes), ]

hclust_neg <- hclust(dist(mat_neg))
hclust_oc  <- hclust(dist(mat_oc))

ordered_genes <- c(
  rownames(mat_neg)[hclust_neg$order],
  rownames(mat_oc)[hclust_oc$order]
)

mat_final <- heat_matrix[ordered_genes, ]

# Row and column annotations
annotation_row <- data.frame(
  Function = ifelse(rownames(mat_final) %in% neg_regulators, 
                    "Negative Regulators / Checkpoints", 
                    "Osteoclast Markers & Resorption"),
  row.names = rownames(mat_final)
)

annotation_col <- data.frame(
  Group = factor(c(rep("Control", 3), rep("Lg", 3)), levels = c("Control", "Lg")),
  row.names = colnames(mat_final)
)

ann_colors <- list(
  Group = c("Control" = "#00E5EE", "Lg" = "#FF8282"),
  Function = c(
    "Negative Regulators / Checkpoints" = "#D73027",
    "Osteoclast Markers & Resorption"   = "#313695"
  )
)

# Export plot data table
write.csv(mat_final, "plot_data_tables/SourceData_Heatmap_Two_Modules.csv", row.names = TRUE)

# Render and save Heatmap
pheatmap(
  mat               = mat_final,
  scale             = "row",
  color             = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),
  breaks            = seq(-1.5, 1.5, length.out = 101),
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  gaps_row          = sum(annotation_row$Function == "Negative Regulators / Checkpoints"),
  annotation_col    = annotation_col,
  annotation_row    = annotation_row,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 9.5,
  fontsize_col      = 11,
  main              = "Osteoclast Repression Modules (Control vs Lg)",
  filename          = "Fig_Heatmap_Osteoclast_Two_Modules.pdf",
  width             = 6.5,
  height            = 7
 
# 4. Mirrored Parallel Bar Plot: RNA-seq vs Promoter H3K9me3 

# Literature-curated comprehensive osteoclast regulator dataset
osteoclast_db <- c(
  "Nfatc1", "Fos", "Jun", "Junb", "Jund", "Mitf", "Rela", "Relb", "Nfkb1", "Nfkb2", 
  "Cbfb", "Runx1", "Sp1", "Pparg", "Creb1", "Batf", "Batf3",
  "Irf8", "Mafb", "Bcl6", "Spen", "Btg3", "Egr1", "Cdkn1a", "Foxo1", "Foxo3", 
  "Foxo4", "Id2", "Nrip1", "Rgs18", "Rgs10", "Trim24", "Stat1", "Sox9",
  "Tnfrsf11a", "Csf1r", "Csf1", "Oscar", "Trem2", "Tyrobp", "Fcer1g", 
  "Sirpa", "Cd44", "Syk", "Btk", "Tec", "Traf6", "Traf2", "Traf3", "Mapk1", "Mapk3", 
  "Mapk8", "Mapk14", "Akt1", "Pik3ca", "Src", "Cbl", "Vav3", "Rac1", "Cyld",
  "Dcstamp", "Ocstamp", "Itgb3", "Itgav", "Atp6v0d2", "Cxcr4", "Ccr1", "Sync1",
  "Ctsk", "Acp5", "Mmp9", "Car2", "Clcn7", "Ostm1", "Atp6v1h", "Atp6v0a1", "Atp6v0a3", 
  "Calcr", "Plau", "Slc9a1"
) %>% unique()

# Gene-level representative peak selection (prioritizing max absolute fold change)
promoter_clean <- Samples_peakAnno_matrix_data %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  arrange(desc(abs(Fold)), p.value, abs(distanceToTSS)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

oc_plot_data <- tibble(gene = osteoclast_db) %>%
  left_join(
    res_df %>% dplyr::select(gene, RNA_log2FC = log2FoldChange, RNA_pval = pvalue, RNA_padj = padj),
    by = "gene"
  ) %>%
  left_join(
    promoter_clean %>% dplyr::select(gene = SYMBOL, H3K9_log2FC = Fold, H3K9_pval = p.value, H3K9_FDR = FDR),
    by = "gene"
  ) %>%
  filter(!is.na(RNA_log2FC)) %>%
  mutate(
    H3K9_log2FC = ifelse(is.na(H3K9_log2FC), 0, H3K9_log2FC),
    H3K9_pval   = ifelse(is.na(H3K9_pval), 1, H3K9_pval)
  )

inspect_genes <- c(
  # Inhibitory checkpoints
  "Egr1", "Btg3", "Irf8", "Id2", "Spen", "Mafb", "Nrip1", "Foxo1", "Bcl6",
  # Osteoclast differentiation markers and effectors
  "Car2", "Calcr", "Ocstamp", "Acp5", "Ctsk", "Dcstamp", "Mmp9", "Nfatc1", "Fos", "Atp6v0d2"
)

bar_df <- oc_plot_data %>%
  filter(gene %in% inspect_genes) %>%
  mutate(
    Group = factor(
      ifelse(gene %in% c("Egr1", "Btg3", "Irf8", "Id2", "Spen", "Mafb", "Nrip1", "Foxo1", "Bcl6"),
             "Negative Regulators",
             "Positive Regulators"),
      levels = c("Negative Regulators", "Positive Regulators")
    ),
    Is_Irf8 = (gene == "Irf8")
  ) %>%
  arrange(Group, RNA_log2FC) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

# Highlight Irf8 specifically on Y-axis
label_colors <- ifelse(levels(bar_df$gene) == "Irf8", "#B71C1C", "black")
label_faces  <- ifelse(levels(bar_df$gene) == "Irf8", "bold.italic", "italic")

# Export plot data table
write.csv(bar_df, "plot_data_tables/SourceData_Mirrored_Parallel_Bar.csv", row.names = FALSE)

# Left panel: Transcriptome shift
p_left <- ggplot(bar_df, aes(x = RNA_log2FC, y = gene)) +
  geom_col(
    aes(fill = RNA_log2FC > 0, color = Is_Irf8, linewidth = Is_Irf8),
    width = 0.72
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#313695"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "#B71C1C", "FALSE" = "black"), guide = "none") +
  scale_linewidth_manual(values = c("TRUE" = 0.65, "FALSE" = 0.3), guide = "none") +
  labs(x = expression(Log[2]~Fold~Change~(RNA-seq)), y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 9.5, color = "black"),
    strip.background = element_rect(fill = "#F2F4F7", color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.y = element_text(color = label_colors, face = label_faces, size = 11),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
  ) +
  ggtitle("Transcriptome Shift")

# Right panel: Epigenetic shift (synchronized red/blue encoding)
p_right <- ggplot(bar_df, aes(x = H3K9_log2FC, y = gene)) +
  geom_col(
    aes(fill = H3K9_log2FC > 0, color = Is_Irf8, linewidth = Is_Irf8),
    width = 0.72
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#313695"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "#B71C1C", "FALSE" = "black"), guide = "none") +
  scale_linewidth_manual(values = c("TRUE" = 0.65, "FALSE" = 0.3), guide = "none") +
  labs(x = expression(Log[2]~Fold~Change~(H3K9me3)), y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
  ) +
  ggtitle("Epigenetic Shift (Promoter)")

p_parallel <- p_left + p_right + plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    title = "Direct Parallel Assessment: RNA-seq vs. H3K9me3 Peak Dynamics",
    subtitle = "Comparing expression change with promoter heterochromatin status across osteoclast regulators",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

ggsave("Fig_Mirrored_Parallel_Bar.pdf", plot = p_parallel, width = 8, height = 6.8)
ggsave("Fig_Mirrored_Parallel_Bar.png", plot = p_parallel, width = 8, height = 6.8, dpi = 300)