#' Process Expression Data
#'
#' This function processes expression data by mapping probe IDs to gene symbols
#' and selecting the probe with the highest average expression for each gene.
#'
#' @param data A data frame or matrix of expression data, with probe IDs as row names
#'   and samples as columns.
#' @param platform A character string specifying the microarray platform (e.g., "gpl570").
#'
#' @return A data frame of processed expression data, with gene symbols as row names
#'   and samples as columns. For each gene, only the probe with the highest average
#'   expression is retained.
#'
#' @export
#'
#' @import dplyr
#' @import idmap1
#' @import idmap2
#'
#' @examples
#' \dontrun{
#' # Assuming you have your expression data in a data frame called 'expr_data'
#' processed_data <- process_expression_data(expr_data, platform = "gpl570")
#' }
#'
process_expression_data <- function(data, platform) {
  # 获取指定平台的 ID 映射
  ids <- getIDs(platform)

  # 处理输入的表达数据
  expr <- as.data.frame(data)
  expr$probe_id <- rownames(expr)

  # 合并表达数据和 ID 映射
  exprSet <- merge(expr, ids)

  # 移除不需要的列
  exprSet[,c('probe_id','gpl')] <- NULL

  # 计算平均值并选择每个基因表达最高的探针
  exprSet <- as.data.frame(exprSet) %>%
    mutate(avg_value = rowMeans(dplyr::select(., -symbol), na.rm = TRUE)) %>%
    group_by(symbol) %>%
    slice_max(order_by = avg_value) %>%
    ungroup() %>%
    dplyr::select(-avg_value)

  return(exprSet)
}
