#' Selection of genes based on statistical significance values computed through Bootstrap-Support Vector Machine (SVM)-Maximum Relevance and Minimum Redundancy (MRMR) methods.
#'
#' @name TopGenesPvalSVMRMR
#' @aliases TopGenesPvalSVMRMR
#' @usage TopGenesPvalSVMRMR(x, y, method, beta, nboot, p.adjust.method, n)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.
#' @param method Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.
#' @param beta Scalar representing trade-off between SVM and MRMR weights.
#' @param nboot Scalar representing the number of bootsrap samples to be drawn from the data using simple random sampling with replacement (Bootstrap) procedure.
#' @param p.adjust.method Character representing the method used for multiple hypothesis correction and computation of adjusted p-values. It can be any method out of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr".
#' @param n Numeric constant (< N) representing the number of top ranked genes to be selected from the high dimensional gene expression data.
#'
#' @return A list of differentially expressed specified number of genes through BSM method.
#'
#' @description The function selects the top ranked genes from the high dimensional gene expression data using the statistical significance values computed through Bootstrap-Support Vector Machine-Maximum Relevance and Minimum Redundancy (BSM) approach.
#'
#' @details Selection of genes based on statistical significance values computed through BSM approach.
#' Takes the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs.
#'
#' @author Samarendra Das <samarendra4849 at gamil.com>
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' TopGenesPvalSVMRMR(x, y, method="Linear", beta=0.6, nboot=20, p.adjust.method = "BH", n=5)
#'
#' @export


TopGenesPvalSVMRMR <- function (x, y, method, beta, nboot, p.adjust.method, n)
{
  this.call = match.call()
  if ((!class(n) == "numeric" & n > nrow(x))) {
    warning("n must be numeric and it should be less than number of rows of x")
  }
  genes.weight <- pvalsvmmrmr (x, y, method, beta, nboot, p.adjust.method, plot=FALSE)[,3]
  id <- sort(genes.weight, decreasing=FALSE, index.return = TRUE)$ix
  TopGenes <- names(genes.weight) [id] [1:n]
  return(TopGenes)
}
############################ TopGenesBootSVMRMR Ends here ##########################################
