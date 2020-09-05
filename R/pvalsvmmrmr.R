#' Computation statistical significance values for gene selection through Bootstrap-Support Vector Machine (SVM)-Maximum Relevance and Minimum Redundancy (MRMR) methods.
#'
#' @name pvalsvmmrmr
#' @aliases pvalsvmmrmr
#' @usage pvalsvmmrmr(x, y, method, beta, nboot, p.adjust.method, plot)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.
#' @param method Character variable representing either 'Linear' or 'Quadratic' method for integrating the weights/scores computed through SVM and MRMR methods.
#' @param beta Scalar representing trade-off between SVM and MRMR weights.
#' @param nboot Scalar representing the number of bootsrap samples to be drawn from the data using simple random sampling with replacement (Bootstrap) procedure.
#' @param p.adjust.method Character representing the method used for multiple hypothesis correction and computation of adjusted p-values. It can be any method out of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr".
#' @param plot Character representing whether the plots need to be drawn or not. It can take value either TRUE or FALSE.
#'
#' @return This function returns returns a vector of rank scores, p-values and adjusted p-values for all genes in the gene expression data using BSM approach.
#'
#' @description The function computes ths rank scores, statisical significance values, adjusted p-values of gene informativeness to a particular condition using the Bootstrap-Support Vector Machine (SVM)-Maximum Relevance and Minimum Redundancy (MRMR) (BSM) approach from high dimensional gene expression data .
#'
#' @details Computation statistical significance values for gene selection through BSM approach.
#' Takes the gene expression data matrix (rows as genes and coloumns as samples) and vector of class labels of subjects (1: case and -1: control) as inputs.
#'
#' @author Samarendra Das <samarendra4849 at gamil.com>
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' pvalsvmmrmr(x, y, method="Linear", beta=0.6, nboot=20, p.adjust.method = "BH", plot=FALSE)
#'
#' @import graphics
#' @import stats
#'
#' @export

pvalsvmmrmr <- function (x, y, method, beta, nboot, p.adjust.method, plot)
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
  if (!class(p.adjust.method)=="character" & is.na(match(p.adjust.method,
                                                         c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")))==TRUE){
    stop("Choose proper multiple hypotheis correction method")
  }
  if (is.null(beta))
    beta = 0.5
  if (!is.null(beta))
    beta = as.numeric(beta)
  if (is.null(method))
    method = "Linear"
  if(is.null(p.adjust.method))
    p.adjust.method="BH"
  if(is.null(plot))
    plot=FALSE
  y <- as.numeric(y)
  genes <- rownames(x)
  #g <- as.matrix(x)
  nn <- nrow(x)
  M <- ncol(x)
  #GeneRankedList <- vector(length = nn)
  RankSCor <- matrix(0, nn, nboot)
  BootScor <- matrix(0, nn, nboot)
  if(missing(beta)){
    beta=0.5
  }

  for (i in 1:nboot) {
    samp <- sample(M, M, replace = TRUE)
    x1 <- x[, samp]
    y1 <- y[samp]
    aa <- as.numeric(weightmr(x1, y1)$weight.mr)
    bb <- as.numeric(weightsvm(x1, y1)$SV.weight)
    wt.mr <- (aa-min(aa))/(max(aa)-min(aa))
    wt.sv <- (bb-min(bb))/(max(bb)-min(bb))

    if(method=="Linear")
    {
      scor <- beta*as.numeric(wt.mr) + (1-beta)*as.numeric(wt.sv)
      #names(score) <- genes
    }

    if (method=="Quadratic"){
      idx1 <- rank(as.numeric(wt.mr), na.last = TRUE, ties.method = "random")
      idx2 <- rank(as.numeric(wt.sv), na.last = TRUE, ties.method = "random")
      scor <- (beta * idx1 * wt.mr + (1-beta) * idx2 * wt.sv)/(beta * idx1+(1-beta) * idx2)
    }
    #names(scor) <- genes
    BootScor[,i] <- scor
    ranks <- rank(-scor, na.last = TRUE, ties.method = "random")
    rankscor <- (nn+1-ranks)/nn
    RankSCor[,i] <- rankscor
  }

  Q3 <- 0.75
  D <- RankSCor - Q3
  pval <- vector(mode = "numeric", length = nn)
  for (i in 1:nn) {
    Z <- D[i, ]
    Z <- Z[Z != 0]
    n1 <- length(Z)
    r <- rank(abs(Z), na.last = TRUE, ties.method = "random")
    tplus <- sum(r[Z > 0])
    etplus <- n1 * (n1 + 1)/8
    vtplus <- n1 * (n1 + 1) * (2 * n1 + 1)/32
    p.value = pnorm(tplus, etplus, sqrt(vtplus), lower.tail = FALSE)
    pval[i] = p.value
  }

  if (is.na(p.adjust.method)==FALSE){
    methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
    idm <- match( p.adjust.method, methods)
    pval.adj <- p.adjust(pval, methods[idm], length(pval))
  }
  res <- cbind(apply(RankSCor, 1, sum), pval, pval.adj)
  rownames(res) <- genes
  colnames(res) <- c("Rank score", "P-value", "Adjusted P-value")
  if (plot == TRUE) {
    par(mfrow=c(2,2))
    #as.vector(res[, 1])
    ranks <- sort(as.vector(res[, 1]), decreasing = TRUE, index.return = TRUE)$ix
    plot(1:length(ranks), as.vector(res[, 1])[ranks], main = "Bootstrap Score vs. gene ranks",
         type = "p", xlab = "Genes", ylab = "BootScores of genes", col="red", pch=9)
    lines(lowess(as.vector(res[, 1])[ranks]), col="blue")

    #as.vector(res[, 2])
    ranks <- sort(as.vector(res[, 2]), decreasing = TRUE, index.return = TRUE)$ix
    plot(1:length(ranks), as.vector(res[, 2])[ranks], main = "Rank Score vs. gene ranks",
         type = "p", xlab = "Genes", ylab = "RankScores of genes", col="red", pch=8)
    lines(lowess(as.vector(res[, 2])[ranks]), col="blue")

    #as.vector(res[, 3])
    ranks <- sort(as.vector(res[, 3]), decreasing = FALSE, index.return = TRUE)$ix
    plot(1:length(ranks), as.vector(res[, 3])[ranks], main = "Pvalues vs. gene ranks",
         type = "p", xlab = "Genes", ylab = "P-values of genes", col="red", pch=7)
    lines(lowess(as.vector(res[, 3])), col="blue")

    #as.vector(res[, 3])
    ranks <- sort(as.vector(res[, 4]), decreasing = FALSE, index.return = TRUE)$ix
    plot(1:length(ranks), round(as.vector(res[, 4])[ranks], 3), main = "Adjusted values vs. gene ranks",
         type = "p", xlab = "Genes", ylab = "Adj P-values of genes", col="red", pch=3)
    lines(lowess(as.vector(res[, 4])), col="blue")
  }
  return(res)
}
