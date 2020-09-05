#' Computation weights for gene selection through Support Vector Machine  and Maximum Relevance and Minimum Redundancy (MRMR) methods.
#'
#' @name weightsvmmrmr
#' @aliases weightsvmmrmr
#' @usage weightsvmmrmr(x, y, method, beta)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.
#' @param method Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.
#' @param beta scalar representing trade-off between SVM and MRMR weights.
#'
#' @return This function produces a vector weights represents the strength of gene informativeness from gene expression data by integrating scores from SVM and MRMR methods using linear and quadartic techniques of score integration.
#'
#' @description The function computes the gene selection weights for genes through a hybrid approach of Support Vector Machine (SVM) and Maximum Relevance and Minimum Redundancy (MRMR) algorithms using the high dimensional gene expression data.
#'
#' @details Computation weights for gene selection through SVM-MMRMR method.
#' Takes the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs.
#'
#' @references Mundra,P.A. &  Rajapakse, J.C. (2010).SVM-RFE with MRMR filter for gene selection. IEEE Trans. Nanobiosci., 9 (1) 31-37, 10.1109/TNB.2009.2035284.
#'
#' @author Samarendra Das <samarendra4849 at gamil.com>
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' weightsvmmrmr(x, y, method="Linear", beta=0.6)
#'
#' @export
#'

weightsvmmrmr <- function (x, y, method, beta)
{

  this.call = match.call()
  if ((!class(x) == "data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }
  if ((!class(y) == "numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }
  if (!length(y) == ncol(x)) {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if (!class(method)=="character" &
      is.na(match(method, c("Linear", "Quadratic")))==TRUE){
    stop("method of score integration must be either Linear or Quadratic")
  }
  if (!class(beta)=="numeric" & (beta < 0 & beta > 1)) {
    warning("Beta must be numeric and its value must be between 0 and 1")
  }

  if (is.null(beta))
    beta = 0.5
  if (!is.null(beta))
    beta = as.numeric(beta)
  if (is.null(method))
    method = "Linear"

  y <- as.numeric(y)
  genes <- rownames(x)

  if (method=="Linear"){
    MR.weight <- weightmr(x, y)$weight.mr
    SV.weight <- weightsvm(x, y)$SV.weight
    genes <- names(MR.weight) <- names(SV.weight)
    score <- beta*MR.weight + (1-beta)*SV.weight
    names(score) <- genes

  }

  if (method=="Quadratic"){
    MR.weight <- weightmr(x, y)$weight.mr
    SV.weight <- weightsvm(x, y)$SV.weight
    genes <- names(MR.weight) <- names(SV.weight)
    idx1 <- sort(as.numeric(MR.weight), decreasing=FALSE, index.return=TRUE)$ix
    idx2 <- sort(as.numeric(SV.weight), decreasing=FALSE, index.return=TRUE)$ix
    score <- (beta * idx1 * MR.weight + (1-beta) * idx2 * SV.weight)/(beta * idx1 + (1-beta) * idx2)
    names(score) <- genes

  }
  class(score) <- "Combined score"
  return(list(score=score, method=method))
}
