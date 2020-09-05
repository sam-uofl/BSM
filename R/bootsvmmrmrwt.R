#' Computation weights for gene selection through Bootstrap (Boot), Support Vector Machine (SVM)  and Maximum Relevance and Minimum Redundancy (MRMR) methods.
#'
#' @name bootsvmmrmrwt
#' @aliases bootsvmmrmrwt
#' @usage bootsvmmrmrwt(x, y, method, beta, nboot)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.
#' @param method Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.
#' @param beta scalar representing trade-off between SVM and MRMR weights.
#' @param nboot scalar representing the number of bootsrap samples to be drawn from the data using simple random sampling with replacement (Bootstrap) procedure.
#'
#' @return This returns a vector of bootstrap weights for gene selection computed through BSM approach.
#'
#' @description The function computes the integrated bootstrap weights for gene selection using the Bootstrap-Support Vector Machine (SVM)-Maximum Relevance and Minimum Redundancy (MRMR) (BSM) approach from high dimensional gene expression data .
#'
#' @details Computation weights for gene selection through BSM approach.
#' Takes the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs.
#'
#' @author Samarendra Das <samarendra4849 at gamil.com>
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' bootsvmmrmrwt(x, y, method="Linear", beta=0.6, nboot=20)
#'
#' @export
#'
#'
bootsvmmrmrwt <- function (x, y, method, beta, nboot)
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
  if ((!class(nboot) == "numeric")) {
    warning("Number of Bootstrap samples must be numeric")
  }
  if (nboot < 0 & nboot <= 50) {
    warning("Number of Bootstrap samples must be numeric and sufficiently large")
  }
  if (is.null(beta))
    beta = 0.5
  if (!is.null(beta))
    beta = as.numeric(beta)
  if (is.null(method))
    method = "Linear"
  y <- as.numeric(y)
  genes <- rownames(x)
  nn <- nrow(x)
  M <- ncol(x)
  #GeneRankedList <- vector(length = nn)
  boot.mr <- matrix(0, nn, nboot)
  boot.sv <- matrix(0, nn, nboot)
  for (j in 1:nboot) {
    samp <- sample(M, M, replace = TRUE)
    x1 <- x[, samp]
    y1 <- y[samp]
    boot.mr[,j] <- as.numeric(weightmr(x1, y1)$weight.mr)
    boot.sv[,j] <- as.numeric(weightsvm(x1, y1)$SV.weight)
  }
  MR.weight <- as.vector(weightmr(x, y)$weight.mr)
  SV.weight <- as.vector(weightsvm(x, y)$SV.weight)
  MR.weight1 <- MR.weight/(apply(boot.mr, 1, mean))
  SV.weight1 <- SV.weight/(apply(boot.sv, 1, mean))
  names(MR.weight1) <- names(SV.weight1) <- genes
  if (method=="Linear"){

    score <- beta*MR.weight1 + (1-beta)*SV.weight1
    names(score) <- genes
  }

  if (method=="Quadratic"){

    idx1 <- rank(as.numeric(MR.weight), na.last = TRUE, ties.method = "random")
    idx2 <- rank(as.numeric(SV.weight), na.last = TRUE, ties.method = "random")
    score <- (beta * idx1 * MR.weight1 + (1-beta) * idx2 * SV.weight1)/(beta * idx1+(1-beta) * idx2)
    names(score) <- genes
    #return(list(score=score, method="Quadratic"))
  }
  class(score) <- "Combined score through Bootstrap"
  return(list(score=score, method=method))
}

###ends here
