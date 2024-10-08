% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_limma_degs.R
\name{perform_limma_analysis}
\alias{perform_limma_analysis}
\title{Perform Differential Expression Analysis using limma}
\usage{
perform_limma_analysis(exprSet, group, con_group, treat_group)
}
\arguments{
\item{exprSet}{A numeric matrix of expression values, with genes in rows and samples in columns.}

\item{group}{A vector specifying the group for each sample. Must be the same length as the number of columns in exprSet.}

\item{con_group}{A character string specifying the name of the control group.}

\item{treat_group}{A character string specifying the name of the treatment group.}
}
\value{
A data frame containing the differential expression analysis results.
The data frame includes log-fold changes, p-values, and adjusted p-values for each gene.
}
\description{
This function conducts differential expression analysis on gene expression data
using the limma package. It fits a linear model and applies empirical Bayes statistics.
}
\examples{
\dontrun{
# Assuming you have your expression matrix (exprSet) and group information
results <- perform_limma_analysis(exprSet, group, "control", "treatment")
}

}
