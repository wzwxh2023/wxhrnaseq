#' Process Differential Expression Results
#'
#' This function processes differential expression results by converting
#' ENSEMBL IDs to ENTREZ IDs and gene symbols. It also performs various
#' checks on the data and provides informative messages.
#' @name ensemble2symbol
#' @title Convert Ensembl IDs to Gene Symbols
#' @description This function converts Ensembl IDs to gene symbols.
#' @param diff_results A data frame containing differential expression results.
#'   Must include a column named 'gene_id' with ENSEMBL IDs.
#' @param species A character string specifying the species.
#'   Must be either "human" or "mouse". Default is "human".
#'
#' @return A data frame with the original differential expression results,
#'   plus additional columns for ENTREZ IDs and gene symbols.
#'
#' @export
#'
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#'
#' @examples
#' \dontrun{
#' # Assuming diff_results is your differential expression results data frame
#' processed_results <- process_diff_results(diff_results, species = "human")
#' }
#'

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

ensemble2symbol <- function(diff_results, species = "human") {
  # 选择正确的物种数据库
  OrgDb <- switch(tolower(species),
                  "human" = org.Hs.eg.db,
                  "mouse" = org.Mm.eg.db,
                  stop("Unsupported species. Choose 'human' or 'mouse'.")
  )

  # 步骤1：移除ENSEMBL ID中的版本号（"."及其后面的数字）
  diff_results$gene_id <- sub("\\..*", "", diff_results$gene_id)

  # 检查重复的ENSEMBL ID
  duplicate_ids <- diff_results$gene_id[duplicated(diff_results$gene_id)]
  if (length(duplicate_ids) > 0) {
    cat("Warning: Found", length(duplicate_ids), "duplicate ENSEMBL IDs.\n")
    cat("First few duplicates:", head(duplicate_ids), "\n")
  }

  # 检查ENSEMBL ID是否有效
  valid_ensembl <- keys(OrgDb, keytype = "ENSEMBL")
  valid_ids <- diff_results$gene_id %in% valid_ensembl

  cat("Number of valid ENSEMBL IDs:", sum(valid_ids), "\n")
  cat("Number of invalid ENSEMBL IDs:", sum(!valid_ids), "\n")

  if (sum(!valid_ids) > 0) {
    cat("First few invalid IDs:", head(diff_results$gene_id[!valid_ids]), "\n")
  }

  # 步骤2：使用clusterProfiler的bitr函数转换基因ID
  gene_ids <- bitr(unique(diff_results$gene_id),
                   fromType = "ENSEMBL",
                   toType = c("ENTREZID", "SYMBOL"),
                   OrgDb = OrgDb)

  # 检查一对多映射
  multi_map <- gene_ids[duplicated(gene_ids$ENSEMBL) | duplicated(gene_ids$ENSEMBL, fromLast = TRUE), ]
  if (nrow(multi_map) > 0) {
    cat("Warning: Found", nrow(multi_map), "ENSEMBL IDs with multiple mappings.\n")
    print(head(multi_map))
  }

  # 步骤3：使用merge合并转换后的基因ID与原始差异表达结果
  merged_results <- merge(diff_results, gene_ids,
                          by.x = "gene_id", by.y = "ENSEMBL",
                          all.x = TRUE)

  # 步骤4：重新排列列的顺序
  final_results <- merged_results[, c("gene_id", "ENTREZID", "SYMBOL",
                                      setdiff(names(merged_results),
                                              c("gene_id", "ENTREZID", "SYMBOL")))]

  cat("Original number of rows:", nrow(diff_results), "\n")
  cat("Final number of rows:", nrow(final_results), "\n")

  return(final_results)
}
