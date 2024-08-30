#' Count Differentially Expressed Genes (DEGs)
#'
#' This function counts the number of up-regulated and down-regulated genes
#' based on given fold change thresholds and p-value threshold.
#'
#' @param diff_results A data frame containing differential expression results.
#'   Must include columns for log fold change and adjusted p-value.
#' @param fc_thresholds A numeric vector of fold change thresholds.
#'   Default is c(1, 1.5, 2).
#' @param p_threshold Numeric. The p-value threshold for significance.
#'   Default is 0.05.
#'
#' @return This function doesn't return a value, but prints the counts of
#'   up-regulated and down-regulated genes for each threshold.
#'
#' @export
#'
#' @importFrom dplyr filter mutate case_when sym
#' @importFrom rlang !!
#'
#' @examples
#' # Assuming diff_results is your differential expression results data frame
#' # count_degs(diff_results)
#' # count_degs(diff_results, fc_thresholds = c(1, 2, 3), p_threshold = 0.01)
#'
count_degs <- function(diff_results, fc_thresholds = c(1, 1.5, 2), p_threshold = 0.05) {
  # Function body remains the same
  # 确定列名
  fc_col <- intersect(c("log2FoldChange", "logFC"), colnames(diff_results))[1]
  p_col <- intersect(c("padj", "adj.P.Val"), colnames(diff_results))[1]

  if(is.na(fc_col) || is.na(p_col)) {
    stop("Cannot find appropriate column names for fold change or adjusted p-value")
  }

  for (fc in fc_thresholds) {
    degs <- diff_results %>%
      filter(!!sym(p_col) < p_threshold) %>%
      mutate(regulation = case_when(
        !!sym(fc_col) >= fc ~ "Up",
        !!sym(fc_col) <= -fc ~ "Down",
        TRUE ~ "Not Significant"
      ))

    up_count <- sum(degs$regulation == "Up")
    down_count <- sum(degs$regulation == "Down")

    cat(sprintf("Threshold: |%s| >= %g, %s < %g\n", fc_col, fc, p_col, p_threshold))
    cat(sprintf("Up-regulated: %d\n", up_count))
    cat(sprintf("Down-regulated: %d\n", down_count))
    cat(sprintf("Total DEGs: %d\n\n", up_count + down_count))
  }
}
