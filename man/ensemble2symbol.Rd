% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_ensemble2symbol.R
\name{ensemble2symbol}
\alias{ensemble2symbol}
\title{Convert Ensembl IDs to Gene Symbols}
\arguments{
\item{diff_results}{A data frame containing differential expression results.
Must include a column named 'gene_id' with ENSEMBL IDs.}

\item{species}{A character string specifying the species.
Must be either "human" or "mouse". Default is "human".}
}
\value{
A data frame with the original differential expression results,
plus additional columns for ENTREZ IDs and gene symbols.
}
\description{
This function converts Ensembl IDs to gene symbols.
}
\details{
Process Differential Expression Results

This function processes differential expression results by converting
ENSEMBL IDs to ENTREZ IDs and gene symbols. It also performs various
checks on the data and provides informative messages.
}
\examples{
\dontrun{
# Assuming diff_results is your differential expression results data frame
processed_results <- process_diff_results(diff_results, species = "human")
}

}
