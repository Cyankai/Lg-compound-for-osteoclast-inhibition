# Integrated Transcriptomic and Epigenetic Analysis of Lg-Mediated Osteoclast Repression
 
This repository contains the complete analytical pipelines, visualization scripts, and source data tables for the multi-omics investigation of the small-molecule compound **Lg** in suppressing osteoclastogenesis. The study integrates bulk RNA-sequencing (RNA-seq) and genome-wide CUT&Tag profiling for repressive histone trimethylation (H3K9me3) in bone marrow-derived macrophages (BMDMs).


## 🚀 Getting Started

### Prerequisites

The analysis was performed in R (version $\ge 4.3.3$). Ensure the following Bioconductor and CRAN packages are installed:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DiffBind", "pheatmap", "RColorBrewer"))
install.packages(c("dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "scales"))
```

---

## 📊 Figure Descriptions

* **Figure 1 (Dotplot)**: Dual-track dotplot illustrating the uncoupling between significant RNA repression and invariant promoter H3K9me3 across canonical osteoclast markers.
* **Figure 2 (Volcano Plot)**: DiffBind-derived volcano plot displaying genome-wide differential H3K9me3 CUT&Tag peaks with exact counts for upregulated, downregulated, and non-significant peaks ($P < 0.05$).
* **Figure 3 (Heatmap)**: Row-standardized ($Z$-score) heatmap partitioned into two functional modules: upstream transcriptional negative checkpoints vs. downstream osteoclast effectors.
* **Figure 4 (Mirrored Bar Plot)**: Bidirectional mirror chart showing parallel transcriptomic ($\text{Log}_2\text{FC}$) and promoter epigenetic ($\text{Log}_2\text{FC}$) changes across osteoclast regulators, with targeted visual emphasis on **Irf8**.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
