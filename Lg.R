#### indirect ####
library(dplyr)
library(ggplot2)
# Save complete results
write.csv(
  res_df,
  "DESeq2_all_results.csv",
  row.names = FALSE
)

write.csv(Samples_peakAnno_matrix_data_promoter,file = "Expression_matrix_peakAnno_promoter.csv")


 ---------------------------------------------------------------
# 1. 经典破骨细胞标志基因集合 (涵盖分化因子、融合蛋白、骨吸收酶类)
 ---------------------------------------------------------------
oc_markers <- c(
  "Nfatc1", "Fos", "Jun",               # 核心分化转录因子
  "Acp5", "Ctsk", "Mmp9", "Car2",       # 破骨酶类与功能标志物
  "Dcstamp", "Ocstamp", "Itgb3",        # 融合与黏附
  "Tnfrsf11a", "Oscar"                  # 破骨受体
)

 ---------------------------------------------------------------
# 2. 提取并整合 RNA-seq 与 CUT&Tag 启动子信号
 ---------------------------------------------------------------
# 提取 CUT&Tag 启动子信息（优先取离 TSS 最近或 p 最显著的 peak）
peak_oc <- Samples_peakAnno_matrix_data_promoter %>%
  filter(SYMBOL %in% oc_markers) %>%
  arrange(p.value, abs(distanceToTSS)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(
    gene = SYMBOL,
    chip_log2FC = Fold,
    chip_pval = p.value,
    chip_p = p
  )

# 提取 RNA-seq 信息
rna_oc <- res_df %>%
  filter(gene %in% oc_markers) %>%
  dplyr::select(
    gene,
    rna_log2FC = log2FoldChange,
    rna_pval = pvalue,
    rna_padj = padj
  )

# 合并并整理成长表格格式（便于 ggplot2 分组绘制双轨 Dotplot）
dot_df <- full_join(rna_oc, peak_oc, by = "gene") %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = c("Assay", ".value"),
    names_pattern = "(rna|chip)_(.*)"
  ) %>%
  mutate(
    Assay = ifelse(Assay == "rna", "RNA-seq (Expression)", "H3K9me3 (Promoter Peak)"),
    # 计算显著性打分 (-log10 P-value)
    neg_log10_pval = -log10(pmax(pval, 1e-10)),
    # 标注显著状态：用于加边框或形状区分
    is_sig = ifelse(pval < 0.05, "p < 0.05", "NS (Not Significant)")
  )

# 基因按 RNA 下调幅度排序
gene_order <- rna_oc %>% arrange(rna_log2FC) %>% pull(gene)
dot_df$gene <- factor(dot_df$gene, levels = gene_order)

 ---------------------------------------------------------------
# 3. 绘制双轨 Dotplot
 ---------------------------------------------------------------
p_dot <- ggplot(dot_df, aes(x = Assay, y = gene)) +
  # 背景网格虚线
  geom_hline(yintercept = seq_along(gene_order), color = "grey92", linewidth = 0.5) +
  
  # 气泡点
  geom_point(
    aes(size = neg_log10_pval, fill = log2FC, color = is_sig),
    shape = 21, stroke = 1
  ) +
  
  # 显著性使用黑色边框，NS 使用浅灰边框
  scale_color_manual(values = c("p < 0.05" = "black", "NS (Not Significant)" = "grey70")) +
  
  # 颜色渐变：蓝色下调，红色上调，白色无变化
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = expression(Log[2]~Fold~Change)
  ) +
  
  scale_size_continuous(
    range = c(2, 9),
    breaks = c(1, 2, 4, 8),
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
    subtitle = "H3K9me3 peaks show no significant increase at promoters of repressed osteoclast genes"
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

print(p_dot)
##### 
 
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(pheatmap)

 ---------------------------------------------------------------
# 1. 严格过滤启动子 H3K9me3 下降且 RNA 上调的基因
 ---------------------------------------------------------------
# 启动子 Peak 去重：每个基因保留最显著且离 TSS 最近的 peak
promoter_down <- Samples_peakAnno_matrix_data_promoter %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(Fold < 0 & p.value < 0.05) %>%
  arrange(p.value, abs(distanceToTSS)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(
    gene = SYMBOL,
    chip_log2FC = Fold,
    chip_pval = p.value,
    chip_p = p,
    distanceToTSS,
    annotation
  )

# RNA 表达上调
rna_up <- res_df %>%
  filter(!is.na(gene)) %>%
  filter(log2FoldChange > 0.25 & pvalue < 0.05) %>%
  dplyr::select(
    gene,
    rna_log2FC = log2FoldChange,
    rna_pval = pvalue,
    rna_padj = padj,
    baseMean
  )

# 整合表观去抑制基因集
derepressed_all <- inner_join(rna_up, promoter_down, by = "gene") %>%
  # 计算去抑制效应综合得分（表观下调幅度与转录上调幅度的几何乘积）
  mutate(derepression_score = rna_log2FC * abs(chip_log2FC)) %>%
  arrange(desc(derepression_score))

 ---------------------------------------------------------------
# 2. 注释破骨细胞负调控因子与转录抑制相关分子
 ---------------------------------------------------------------
# 文献与功能明确的单核/巨噬/破骨负调控因子、细胞周期阻滞与共抑制子
known_oc_inhibitors <- c(
  "Irf8",     # 破骨核心负调控 TF，抑制 Nfatc1 与 AP-1
  "Spen",     # 经典转录共抑制因子（SHARP/NCoR 复合物组分）
  "Btg3",     # 细胞周期与增殖负调控蛋白
  "Egr1",     # 调控巨噬细胞分化偏向、竞争性抑制破骨程序
  "Cdkn1a",   # p21，介导生长阻滞并抑制多核化融合
  "Tnfrsf1b", # TNFR2，介导反向保护性信号
  "Mafb",     # 破骨分化抑制 TF（维持单核/巨噬状态）
  "Bcl6",     # 负调控破骨分化晚期基因
  "Nrip1"     # 核受体相互作用共抑制因子
)

# 标记是否为重点破骨抑制相关基因
derepressed_all <- derepressed_all %>%
  mutate(
    is_known_inhibitor = gene %in% known_oc_inhibitors
  )

write.csv(derepressed_all, "Epigenetically_Derepressed_Candidate_Regulators.csv", row.names = FALSE)
# 提取重点关注的抑制因子（若在去抑制集中则直接提取，可按 rna_log2FC 降序排列）
plot_targets <- promoter_down %>%
  filter(gene %in% known_oc_inhibitors | is_known_inhibitor) %>%
  arrange(rna_log2FC)

# 转换为长格式用于双向绘图
mirror_df <- plot_targets %>%
  dplyr::select(gene, `RNA Up (log2FC)` = rna_log2FC, `H3K9me3 Down (log2FC)` = chip_log2FC) %>%
  pivot_longer(cols = -gene, names_to = "Modality", values_to = "FoldChange")

mirror_df$gene <- factor(mirror_df$gene, levels = plot_targets$gene)

p_mirror <- ggplot(mirror_df, aes(x = FoldChange, y = gene, fill = Modality)) +
  geom_col(width = 0.65, position = "identity", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "RNA Up (log2FC)" = "#D73027",
    "H3K9me3 Down (log2FC)" = "#4575B4"
  )) +
  scale_x_continuous(
    breaks = seq(-0.1, 1.5, by = 0.2),
    labels = function(x) abs(x)
  ) +
  labs(
    x = expression(Absolute~Log[2]~Fold~Change~(LG~vs~Control)),
    y = NULL,
    title = "Epigenetic Derepression of Osteoclast Inhibitory Regulators",
    subtitle = "Promoter H3K9me3 loss (Blue, Left) coincides with transcriptional activation (Red, Right)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.text.y = element_text(face = "bold.italic", size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linetype = "dashed")
  )

print(p_mirror)
# 查看关注基因在两个组学中的真实统计值
check_genes <- c("Cdkn1a", "Irf8", "Spen", "Btg3", "Egr1", "Tnfrsf1b", "Wwox", "Stab1")

plot_df %>% 
  filter(gene %in% check_genes) %>%
  dplyr::select(gene, RNA_FC = log2FoldChange, RNA_p = pvalue, 
                H3K9_FC = h3k9_Fold, H3K9_p = h3k9_pval) %>%
  arrange(desc(RNA_FC))
library(dplyr)
library(ggplot2)
library(tidyr)

 ---------------------------------------------------------------
# 1. 全基因组范围提取表观去抑制基因集（解除名单限制）
 ---------------------------------------------------------------
# 提取启动子 H3K9me3 显著下调的 peak
promoter_down <- Samples_peakAnno_matrix_data_promoter %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(Fold < 0 & p.value < 0.05) %>%
  arrange(p.value, abs(distanceToTSS)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(
    gene = SYMBOL,
    h3k9_FC = Fold,
    h3k9_pval = p.value,
    distanceToTSS
  )

# 提取 RNA 表达显著上调的基因
rna_up <- res_df %>%
  filter(!is.na(gene)) %>%
  filter(log2FoldChange > 0.2 & pvalue < 0.05) %>%
  dplyr::select(
    gene,
    rna_FC = log2FoldChange,
    rna_pval = pvalue
  )

# 合并交集：全基因组表观去抑制候选基因
derepressed_all <- inner_join(rna_up, promoter_down, by = "gene") %>%
  # 综合效应打分：转录上调幅度与表观下调幅度的乘积
  mutate(score = rna_FC * abs(h3k9_FC)) %>%
  arrange(desc(score))

 ---------------------------------------------------------------
# 2. 选取 Top 15 个去抑制基因（或叠加关键破骨抑制基因）
 ---------------------------------------------------------------
# 提取综合得分最高的 Top 15 基因
top_candidates <- derepressed_all %>%
  head(15)

# 如果关注的破骨基因（如 Irf8/Spen/Btg3）在交集总表但未进前15，可强制合并展示
focus_genes <- c("Irf8", "Spen", "Btg3", "Tnfrsf1b", "Cdkn1a")
extra_genes <- derepressed_all %>% filter(gene %in% focus_genes)

# 合并并去重
plot_genes_df <- bind_rows(top_candidates, extra_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(rna_FC)

# 导出清单以备论文附表使用
write.csv(derepressed_all, "All_Derepressed_Genes_List.csv", row.names = FALSE)


##### DiffBind_H3K9me3_Volcano ####
library(DiffBind)
library(dplyr)
library(ggplot2)

 -----------------------------------------------------------------------
# 1. 提取差异分析完整结果（包含所有 peak，不设过滤阈值）
 -----------------------------------------------------------------------
# th = 1 确保提取所有 peak，bCounts = TRUE 保留计数
report_peaks <- dba.report(Samples_analyze, contrast = 1, th = 1, bCounts = FALSE) %>%
  as.data.frame()

 ---------------------------------------------------------------
# 2. 定义严格分类与阈值
 --------------------------------------------------------------------
# 设定阈值：p.value < 0.05 且 |Fold| > 0 (或具体设定如 |Fold| > 0.5)
p.value_threshold  <- 0.4
fold_threshold <- 0    # 若需要更大差异幅度可设为 0.25 或 0.5

report_peaks <- report_peaks %>%
  mutate(
    neg_log10_p = -log10(p.value),
    # 区分 Up / Down / NS
    Status = case_when(
      p.value < p.value_threshold & Fold <= -fold_threshold ~ "Down",
      p.value < p.value_threshold & Fold >=  fold_threshold ~ "Up",
      TRUE                                          ~ "NS"
    )
  )

# 统计具体数量用于角标展示
n_down <- sum(report_peaks$Status == "Down", na.rm = TRUE)
n_up   <- sum(report_peaks$Status == "Up", na.rm = TRUE)
n_ns   <- sum(report_peaks$Status == "NS", na.rm = TRUE)

cat(sprintf("DiffBind 峰统计结果:\n显著下调 (Down): %d\n显著上调 (Up): %d\n无显著变化 (NS): %d\n", 
            n_down, n_up, n_ns))

 ----------------------------------------------------------------------
# 3. 绘制紧凑美观的火山图（标明上下调数量）
 --------------------------------------------------------------------
# 确定横轴对称范围，使图形居中工整
max_abs_fold <- max(abs(report_peaks$Fold), na.rm = TRUE)
x_limit <- ceiling(max_abs_fold * 1.05)
y_limit <- max(report_peaks$neg_log10_p, na.rm = TRUE) * 1.05

p_volcano <- ggplot(report_peaks, aes(x = Fold, y = neg_log10_p)) +
  # 辅助虚线
  geom_vline(xintercept = c(-fold_threshold, fold_threshold), linetype = "dashed", color = "grey65", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey65", linewidth = 0.5) +
  
  # 背景 NS 散点（灰色、较小）
  geom_point(data = filter(report_peaks, Status == "NS"),
             color = "#D1D5DB", alpha = 0.4, size = 1.2) +
  
  # 下调峰（经典冷蓝/深蓝）
  geom_point(data = filter(report_peaks, Status == "Down"),
             color = "#3B82F6", alpha = 0.8, size = 1.8) +
  
  # 上调峰（若有，经典暖红）
  geom_point(data = filter(report_peaks, Status == "Up"),
             color = "#EF4444", alpha = 0.8, size = 1.8) +
  
  # 左上角标明下调总数
  annotate(
    "text",
    x = -x_limit * 0.92,
    y = y_limit * 0.95,
    label = paste0("Downregulated Peaks\n(p < 0.4)\nn = ", n_down),
    color = "#1D4ED8",
    fontface = "bold",
    size = 4.2,
    hjust = 0,
    vjust = 1
  ) +
  
  # 右上角标明上调总数
  annotate(
    "text",
    x = x_limit * 0.92,
    y = y_limit * 0.95,
    label = paste0("Upregulated Peaks\n(p < 0.4)\nn = ", n_up),
    color = "#B91C1C",
    fontface = "bold",
    size = 4.2,
    hjust = 1,
    vjust = 1
  ) +
  
  # 坐标轴与范围控制
  scale_x_continuous(
    limits = c(-x_limit, x_limit),
    breaks = seq(-x_limit, x_limit, by = 1)
  ) +
  scale_y_continuous(
    limits = c(0, y_limit),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  labs(
    x = expression(Log[2]~Fold~Change~(KO~vs~WT)),
    y = expression(-Log[10]~(p)),
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

# 输出预览
print(p_volcano)

 ----------------------------------------------------------------------
# 4. 导出黄金比例尺寸（解决“太大/太扁”的问题）
 ----------------------------------------------------------------------
# 标准正方形或 4:3 比例，建议宽 6 英寸、高 5 英寸
ggsave("DiffBind_H3K9me3_Volcano_Clean.pdf", plot = p_volcano, width = 6, height = 5)
ggsave("DiffBind_H3K9me3_Volcano_Clean.png", plot = p_volcano, width = 6, height = 5, dpi = 300)



library(dplyr)

 ---------------------------------------------------------------------
# 方案 A：直接从完整注释矩阵中提取所有显著下调的基因（推荐）
 ---------------------------------------------------------------------
# 1. 全基因组所有下调 Peak 对应的基因（含启动子、增强子、内含子区）
down_genes_all <- Samples_peakAnno_matrix_data %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(Fold < 0 & p.value < 0.4) %>%                  # 若 p 过于严苛，可替换为 p.value < 0.05
  arrange(p.value, Fold) %>%                             # 按显著性与下调幅度排序
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(
    SYMBOL, 
    GENENAME, 
    Fold, 
    p.value, 
    annotation, 
    distanceToTSS,
    seqnames, start, end
  )

# 2. 仅限于启动子区（Promoter）显著下调的基因
down_genes_promoter <- Samples_peakAnno_matrix_data_promoter %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(Fold < 0 & p.value < 0.4) %>%
  arrange(p.value, Fold) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(SYMBOL, GENENAME, Fold, p.value, annotation, distanceToTSS)

table(down_genes_promoter %in% osteoclast_key_genes)

 -------------------------------------------------------------------
# 方案 B：若必须直接使用刚生成的 report_peaks 对象提取
 ---------------------------------------------------------------------
# 将 report_peaks 与注释表的基因坐标对齐合并
if (!"SYMBOL" %in% colnames(report_peaks)) {
  report_peaks_annotated <- report_peaks %>%
    inner_join(
      Samples_peakAnno_matrix_data %>% 
        dplyr::select(seqnames, start, end, SYMBOL, GENENAME, annotation, distanceToTSS),
      by = c("seqnames", "start", "end")
    )
} else {
  report_peaks_annotated <- report_peaks
}

down_genes_from_report <- report_peaks_annotated %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(Fold < 0 & p.value < 0.4) %>%
  arrange(p.value, Fold) %>%
  distinct(SYMBOL, .keep_all = TRUE)

table(down_genes_from_report %in% osteoclast_key_genes)

#######
library(dplyr)
library(tidyr)
library(ggplot2)

 
# 1. 差异阈值设定
  
rna_padj_cut   <- 0.01    # RNA-seq FDR 阈值 (若需放宽可改为 pvalue < 0.05)
rna_lfc_cut    <- 0    # RNA-seq 变化倍数阈值 

chip_pval_cut  <- 0.4    # CUT&Tag 显著性阈值 (基于此前分析，建议使用 p.value < 0.05)
chip_lfc_cut   <- 0.05      # CUT&Tag Fold 阈值 (Fold < 0 为 H3K9me3 丢失，Fold > 0 为获得)

  
# 2. 提取 RNA-seq 差异表达基因 (DEGs)
 
rna_valid <- res_df %>% 
  filter(!is.na(gene) & gene != "")

rna_up_genes <- rna_valid %>%
  filter(padj < rna_padj_cut & log2FoldChange >= rna_lfc_cut) %>%
  pull(gene) %>% unique()

rna_down_genes <- rna_valid %>%
  filter(padj < rna_padj_cut & log2FoldChange <= -rna_lfc_cut) %>%
  pull(gene) %>% unique()

rna_deg_all <- union(rna_up_genes, rna_down_genes)

 
# 3. 提取 CUT&Tag 启动子区差异基因 (Differential Promoter Peaks)
 
# 启动子区按基因折叠去重
promoter_clean <- Samples_peakAnno_matrix_data %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  arrange(p.value, abs(distanceToTSS)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

chip_down_genes <- promoter_clean %>%
  filter(p.value < chip_pval_cut & Fold < -chip_lfc_cut) %>%
  pull(SYMBOL) %>% unique()

chip_up_genes <- promoter_clean %>%
  filter(p.value < chip_pval_cut & Fold > chip_lfc_cut) %>%
  pull(SYMBOL) %>% unique()

chip_diff_all <- union(chip_down_genes, chip_up_genes)

 
# 4. 构建已知报道的破骨调控基因库 (Literature & Database Curated)
 
# 整合已报道的正向促分化因子、转录负调控因子、融合受体及骨吸收效应酶
osteoclast_db <- c(
  # 核心分化转录因子 (Positive TFs)
  "Nfatc1", "Fos", "Jun", "Junb", "Jund", "Mitf", "Rela", "Relb", "Nfkb1", "Nfkb2", 
  "Cbfb", "Runx1", "Sp1", "Pparg", "Creb1", "Batf", "Batf3",
  # 负调控转录因子与共抑制子 (Negative TFs / Repressors / Checkpoints)
  "Irf8", "Mafb", "Bcl6", "Spen", "Btg3", "Egr1", "Cdkn1a", "Foxo1", "Foxo3", 
  "Foxo4", "Id2", "Nrip1", "Rgs18", "Rgs10", "Trim24", "Stat1", "Sox9",
  # RANKL/M-CSF 信号传导与共刺激通路
  "Tnfrsf11a", "Csf1r", "Csf1", "Oscar", "Trem2", "Tyrobp", "Fcer1g", 
  "Sirpa", "Cd44", "Syk", "Btk", "Tec", "Traf6", "Traf2", "Traf3", "Mapk1", "Mapk3", 
  "Mapk8", "Mapk14", "Akt1", "Pik3ca", "Src", "Cbl", "Vav3", "Rac1", "Cyld",
  # 细胞融合与前体迁移 (Fusion, Migration & Adhesion)
  "Dcstamp", "Ocstamp", "Itgb3", "Itgav", "Atp6v0d2", "Cxcr4", "Ccr1", "Sync1",
  # 骨吸收效应分子与质子泵 (Resorption Machinery)
  "Ctsk", "Acp5", "Mmp9", "Car2", "Clcn7", "Ostm1", "Atp6v1h", "Atp6v0a1", "Atp6v0a3", 
  "Calcr", "Plau", "Slc9a1"
) %>% unique()

 
# 5. 计算各维度交集
  
# A. RNA-seq 差异与破骨基因库交集
overlap_rna_all   <- intersect(rna_deg_all, osteoclast_db)
overlap_rna_up    <- intersect(rna_up_genes, osteoclast_db)
overlap_rna_down  <- intersect(rna_down_genes, osteoclast_db)

# B. CUT&Tag 启动子差异与破骨基因库交集
overlap_chip_all  <- intersect(chip_diff_all, osteoclast_db)
overlap_chip_down <- intersect(chip_down_genes, osteoclast_db)  # 表观去抑制位点
overlap_chip_up   <- intersect(chip_up_genes, osteoclast_db)

# C. 三者交集：表观去抑制驱动的破骨候选靶点 (H3K9me3 Down & RNA Up & OC DB)
overlap_derepressed_oc <- intersect(intersect(chip_down_genes, rna_up_genes), osteoclast_db)

# D. 表观顺式直接关闭的破骨候选靶点 (H3K9me3 Up & RNA Down & OC DB)
overlap_silenced_oc    <- intersect(intersect(chip_up_genes, rna_down_genes), osteoclast_db)



##### RNA-seq Heatmap ####
library(pheatmap)
library(RColorBrewer)
library(dplyr)
 
# 1. 严格限定两组核心靶基因列表与功能注释 
neg_regulators <- c(
  "Irf8", "Spen", "Cdkn1a", "Btg3", "Egr1", "Mafb", "Bcl6", 
  "Nrip1", "Rgs18", "Foxo1", "Foxo3", "Foxo4", "Id2"
)

oc_markers <- c(
  "Ctsk", "Acp5", "Nfatc1", "Calcr", "Car2", "Mmp9", "Dcstamp", 
  "Ocstamp", "Atp6v0d2"
)

# 构建精简注释数据框
target_genes_df <- data.frame(
  gene = c(neg_regulators, oc_markers),
  Function = c(
    rep("Negative Regulators / Checkpoints", length(neg_regulators)),
    rep("Osteoclast Markers & Resorption", length(oc_markers))
  )
)

# 仅保留表达矩阵中存在的基因并保持指定顺序
valid_genes <- intersect(target_genes_df$gene, rownames(heat_matrix))
target_genes_df <- target_genes_df %>% filter(gene %in% valid_genes)

# 分别在各组内部进行聚类，使组内表达模式聚集整齐
mat_neg <- heat_matrix[intersect(neg_regulators, valid_genes), ]
mat_oc  <- heat_matrix[intersect(oc_markers, valid_genes), ]

hclust_neg <- hclust(dist(mat_neg))
hclust_oc  <- hclust(dist(mat_oc))

ordered_genes <- c(
  rownames(mat_neg)[hclust_neg$order],
  rownames(mat_oc)[hclust_oc$order]
)

mat_final <- heat_matrix[ordered_genes, ]
 
# 2. 构建行注释、列注释及配色 
# 行注释（仅保留两类）
annotation_row <- data.frame(
  Function = ifelse(rownames(mat_final) %in% neg_regulators, 
                    "Negative Regulators / Checkpoints", 
                    "Osteoclast Markers & Resorption"),
  row.names = rownames(mat_final)
)

# 列注释（Control 与 Lg 分组）
annotation_col <- data.frame(
  Group = factor(c(rep("Control", 3), rep("Lg", 3)), levels = c("Control", "Lg")),
  row.names = colnames(mat_final)
)

# 配色方案：复刻标准发文样式
ann_colors <- list(
  Group = c("Control" = "#00E5EE", "Lg" = "#FF8282"),
  Function = c(
    "Negative Regulators / Checkpoints" = "#D73027",  # 红色系（去抑制上调）
    "Osteoclast Markers & Resorption"   = "#313695"   # 蓝色系（反式下调）
  )
)
 
# 3. 绘制清晰紧凑的两组别热图 
pheatmap(
  mat               = mat_final,
  scale             = "row",
  color             = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),
  breaks            = seq(-1.5, 1.5, length.out = 101),
  cluster_cols      = F,
  cluster_rows      = FALSE,                                # 保持模块间严格分割
  gaps_row          = sum(annotation_row$Function == "Negative Regulators / Checkpoints"), # 仅在两组交界处切开
  annotation_col    = annotation_col,
  annotation_row    = annotation_row,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 9.5,
  fontsize_col      = 11,
  treeheight_col    = 30,
  filename          = "Osteoclast_Filtered_Two_Modules.pdf",
  width             = 6,
  height            = 5
) 
