% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopGenesSVMRMR.R
\name{TopGenesSVMRMR}
\alias{TopGenesSVMRMR}
\title{Selects the top ranked (differentially expressed) genes through SVM-MRMR method.}
\usage{
TopGenesSVMRMR(x, y, method, beta, n)
}
\arguments{
\item{x}{Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.}

\item{y}{px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.}

\item{method}{Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.}

\item{beta}{scalar representing trade-off between SVM and MRMR weights.}

\item{n}{Numeric constant (< N) representing the number of top ranked genes to be selected from the gene expression data.}
}
\value{
A list of differentially expressed specified number of genes through SVM-MRMR method.
}
\description{
This function returns the list of top ranked genes selected through Support Vector Machine-Maximum Relevance and Minimum Redundancy (SVM-MRMR) method.
}
\details{
Selects the top ranked (differentially expressed) genes through SVM-MRMR method.
Take the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs
}
\examples{
x=as.data.frame(matrix(runif(1000), 50))
row.names(x) = paste("Gene", 1:50)
colnames(x) = paste("Samp", 1:20)
y=as.numeric(c(rep(1, 10), rep(-1, 10)))
TopGenesSVMRMR(x, y, method="Linear", beta=0.6, n=5)

}
\author{
Samarendra Das <samarendra4849 at gamil.com>
}
