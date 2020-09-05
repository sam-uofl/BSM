#' Selection of the top ranked (differentially expressed) genes through SVM method.
#'
#' @name TopGenesSVM
#' @aliases TopGenesSVM
#' @title Selection of top ranked (differentially expressed) genes through SVM method
#' @usage TopGenesSVM(x, y, n)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control condition respectively.
#' @param n Numeric constant (< N) representing the number of top ranked genes to be selected from the gene expression data.
#'
#' @return A list of differentially expressed specified number of genes through SVM method.
#'
#' @description The function returns the list of top ranked genes of specified size selected through Support Vector Machine (SVM) method.
#'
#' @details Selects the top ranked (differentially expressed) genes of specified size through SVM method.
#' Take the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects
#' (1: case and -1: control) as inputs.
#'
#' @author Samarendra Das <samarendra4849 at gmail.com>
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' TopGenesSVM(x, y, n=5)
#'
#' @export
#'
TopGenesSVM <- function(x, y, n){
  this.call = match.call()
  if ((!class(n) == "numeric" & n > nrow(x))) {
    warning("n must be numeric and it should be less than number of rows of x")
  }
  genes.weight <- weightsvm(x, y)$SV.weight
  id <- sort(genes.weight, decreasing=TRUE, index.return = TRUE)$ix
  TopGenes <- genes.weight [id] [1:n]
  return(TopGenes)
}
###############################TopGenesSVM Ends here###############################
