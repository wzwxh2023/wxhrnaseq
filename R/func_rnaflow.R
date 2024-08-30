#' Perform RNA-Seq Analysis
#' @name rna_seq_analysis
#' @description
#' This function performs a comprehensive RNA-Seq analysis including data normalization,
#' PCA analysis, differential expression analysis, and visualization.
#'
#' @param count_matrix A matrix of raw count data, with genes in rows and samples in columns.
#' @param sample_groups A vector specifying the group for each sample.
#' @param control_group The name of the control group.
#' @param treatment_group The name of the treatment group.
#' @param output_dir The directory to save output files. Default is "output".
#'
#' @return A list containing:
#'   \item{diff_results}{A data frame of differential expression results}
#'   \item{pca_results}{A list containing PCA results}
#'
#' @export
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst plotPCA DESeq results lfcShrink
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr %>%
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual geom_hline xlim ylim labs theme_bw theme ggsave
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' # Assuming you have your count matrix and sample information
#' result <- rna_seq_analysis(count_matrix, sample_groups, "control", "treatment")
#' }


library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggrepel)

rna_seq_analysis <- function(count_matrix, sample_groups, control_group, treatment_group, output_dir = "output") {

  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE)

  # 准备数据，明确指定因子水平
  meta <- data.frame(
    sample = colnames(count_matrix)[-1],
    group = factor(sample_groups, levels = c(control_group, treatment_group))
  )

  # 创建DESeqDataSet对象
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = meta,
                                design = ~group,
                                tidy = TRUE)

  # 过滤低表达基因
  dds <- dds[rowSums(counts(dds)) > 1,]

  # VST标准化
  vsd <- vst(dds, blind = FALSE)

  # PCA分析
  pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  # 保存PCA结果
  pca_results <- list(
    pca_data = pca_data,
    percentVar = percentVar
  )
  save(pca_results, file = file.path(output_dir, "pca_results.RData"))

  # 绘制PCA图
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group, shape = group)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw()

  ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 8, height = 6)

  # 保存标准化后的表达矩阵
  exprSet_vst <- as.data.frame(assay(vsd))
  save(exprSet_vst, file = file.path(output_dir, "exprSet_vst.RData"))

  # 差异分析
  dds <- DESeq(dds)

  # 设置对比组
  contrast <- c("group", treatment_group, control_group)

  # 获取差异分析结果
  res <- results(dds, contrast = contrast, alpha = 0.05)

  # 创建MA图函数
  create_ma_plot <- function(res_data, title) {
    df <- as.data.frame(res_data) %>%
      mutate(significant = ifelse(padj < 0.05, "Significant", "Not Significant"))

    ggplot(df, aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlim(0, NA) +
      ylim(-5, 5) +
      labs(x = "log10(Mean of normalized counts)",
           y = "log2FoldChange",
           title = title) +
      theme_bw() +
      theme(legend.position = "bottom")
  }

  # MA图
  ma_plot <- create_ma_plot(res, "MA Plot")
  ggsave(file.path(output_dir, "ma_plot.png"), ma_plot, width = 8, height = 6)

  # logFC矫正
  res_shrink <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

  # 矫正后的MA图
  ma_plot_shrink <- create_ma_plot(res_shrink, "MA Plot (Shrunken LFC)")
  ggsave(file.path(output_dir, "ma_plot_shrink.png"), ma_plot_shrink, width = 8, height = 6)

  # 导出差异分析结果
  diff_results <- res_shrink %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")

  write.csv(diff_results, file = file.path(output_dir, "differential_expression_results.csv"), row.names = FALSE)

  return(list(diff_results = diff_results, pca_results = pca_results))
}
