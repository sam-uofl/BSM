#' Selects the top ranked (differentially expressed) genes through Bootstrap-SVM-MRMR method.
#'
#' @name TopGenesBootSVMRMR
#' @aliases TopGenesBootSVMRMR
#' @usage TopGenesBootSVMRMR(x, y, method, beta, nboot, n)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.
#' @param method Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.
#' @param beta scalar representing trade-off between SVM and MRMR weights.
#' @param nboot scalar representing the number of bootsrap samples to be drawn from the data using simple random sampling with replacement (Bootstrap) procedure.
#' @param n Numeric constant (< N) representing the number of top ranked genes to be selected from the gene expression data.
#'
#' @return A list of differentially expressed specified number of genes through Boot-SVM-MRMR method.
#'
#' @description This function returns the list of top ranked genes selected through Boot-Support Vector Machine-Maximum Relevance and Minimum Redundancy (SVM-MRMR) method.
#'
#' @details Selects the top ranked (differentially expressed) genes through Boot-SVM-MRMR method.
#' Take the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs
#'
#'
#' @author Samarendra Das <samarendra4849 at gamil.com>
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' TopGenesBootSVMRMR(x, y, method="Linear", beta=0.6, nboot=20, n=5)
#'
#' @export
#'

TopGenesBootSVMRMR <- function (x, y, method, beta, nboot, n)
{
  this.call = match.call()
  if ((!class(n) == "numeric" & n > nrow(x))) {
    warning("n must be numeric and it should be less than number of rows of x")
  }
  genes.weight <-  bootsvmmrmrwt (x, y, method, beta, nboot)$score
  id <- sort(as.vector(genes.weight), decreasing=TRUE, index.return = TRUE)$ix
  TopGenes <- names(genes.weight) [id] [1:n]
  return(TopGenes)
}

###########ends here
