#' Perform Differential Expression Analysis using limma
#'
#' This function conducts differential expression analysis on gene expression data
#' using the limma package. It fits a linear model and applies empirical Bayes statistics.
#'
#' @param exprSet A numeric matrix of expression values, with genes in rows and samples in columns.
#' @param group A vector specifying the group for each sample. Must be the same length as the number of columns in exprSet.
#' @param con_group A character string specifying the name of the control group.
#' @param treat_group A character string specifying the name of the treatment group.
#'
#' @return A data frame containing the differential expression analysis results.
#'   The data frame includes log-fold changes, p-values, and adjusted p-values for each gene.
#'
#' @export
#'
#' @import limma
#'
#' @examples
#' \dontrun{
#' # Assuming you have your expression matrix (exprSet) and group information
#' results <- perform_limma_analysis(exprSet, group, "control", "treatment")
#' }
#'

perform_limma_analysis <- function(exprSet, group, con_group, treat_group) {
  # 确保group是一个factor，并且按照指定的顺序设置levels
  group <- factor(group, levels = c(con_group, treat_group))

  # 构建比较矩阵
  design <- model.matrix(~group)
  colnames(design) <- levels(group)

  # 线性模型拟合
  fit <- lmFit(exprSet, design)

  # 贝叶斯检验
  fit2 <- eBayes(fit)

  # 输出差异分析结果
  allDiff <- topTable(fit2, adjust='fdr', coef=2, number=Inf)

  return(allDiff)
}
