% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopGenesPvalSVMRMR.R
\name{TopGenesPvalSVMRMR}
\alias{TopGenesPvalSVMRMR}
\title{Selection of genes based on statistical significance values computed through Bootstrap-Support Vector Machine (SVM)-Maximum Relevance and Minimum Redundancy (MRMR) methods.}
\usage{
TopGenesPvalSVMRMR(x, y, method, beta, nboot, p.adjust.method, n)
}
\arguments{
\item{x}{Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.}

\item{y}{px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.}

\item{method}{Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.}

\item{beta}{Scalar representing trade-off between SVM and MRMR weights.}

\item{nboot}{Scalar representing the number of bootsrap samples to be drawn from the data using simple random sampling with replacement (Bootstrap) procedure.}

\item{p.adjust.method}{Character representing the method used for multiple hypothesis correction and computation of adjusted p-values. It can be any method out of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr".}

\item{n}{Numeric constant (< N) representing the number of top ranked genes to be selected from the high dimensional gene expression data.}
}
\value{
A list of differentially expressed specified number of genes through BSM method.
}
\description{
The function selects the top ranked genes from the high dimensional gene expression data using the statistical significance values computed through Bootstrap-Support Vector Machine-Maximum Relevance and Minimum Redundancy (BSM) approach.
}
\details{
Selection of genes based on statistical significance values computed through BSM approach.
Takes the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs.
}
\examples{
x=as.data.frame(matrix(runif(1000), 50))
row.names(x) = paste("Gene", 1:50)
colnames(x) = paste("Samp", 1:20)
y=as.numeric(c(rep(1, 10), rep(-1, 10)))
TopGenesPvalSVMRMR(x, y, method="Linear", beta=0.6, nboot=20, p.adjust.method = "BH", n=5)

}
\author{
Samarendra Das <samarendra4849 at gamil.com>
}
