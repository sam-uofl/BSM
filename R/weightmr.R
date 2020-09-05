##################### weight.mr starts from here ######################
#' Computation of gene ranking weights through Maximum Relevance and Minimum Redundancy (MRMR) method.
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control conditions respectively.
#'
#' @return A vector of weights computed through MRMR method and a complete ranked gene list.
#'
#' @description The function computes weights for features (e.g. genes) using Maximum Relevance and Minimum Redundancy (MRMR) technique and also produce the ranked genelist based on the computed weights.
#'
#' @details Computation of gene ranking weights through Maximum Relevance and Minimum Redundancy (MRMR) method.
#'
#' @references Feature selection based on mutual information: criteria of max-dependency, max-relevance and min-redundancy. IEEE Trans. Pattern Anal. Mach. Intell., 27 (8) 1226-1238. DOI: 10.1109/TPAMI.2005.159.
#'
#' @author Samarendra Das <samarendra4849 at gmail.com>
#'
#' @importFrom stats cor var
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' weightmr(x, y)
#'
#' @export
#'
weightmr <- function (x, y)
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


  #cls <- as.numeric(y)
  genes <- rownames(x)
  x1 <- as.matrix(x)
  nn <- nrow(x1)
  M <- ncol(x1)
  GeneRankedList <- vector(length = nn)

  qsi <- as.vector((apply(abs(cor(t(x1), method = "pearson",
                                  use = "p") - diag(nrow(x1))), 1, sum))/(nrow(x1) -
                                                                            1))
  idx <- which(y == 1)
  idy <- which(y == -1)
  B = vector(mode = "numeric", nn)
  for (i in 1:nrow(x1))
  {
    f.mes <- (((mean(x1[i, idx]) - mean(x1[i, ]))^2) +
                ((mean(x1[i, idy]) - mean(x1[i, ]))^2))/(var(x1[i,
                                                                idx]) + var(x1[i, idy]))
    B[i] <- f.mes
  }
  rsi <- abs(B)
  weight.mr <- rsi/qsi
  names(weight.mr) <- genes
  GeneRankedList <- sort(-weight.mr, index.return = TRUE)$ix
  #rankvalue <- sort(GeneRankedList, index.return = TRUE)$ix
  #rankscore <- (nn + 1 - rankvalue)/(nn)

  #rankingcriteria <- as.vector(rowSums((M1), na.rm = FALSE, dims = 1))

  names(GeneRankedList) <- genes[GeneRankedList]

  #class(weight.mr ) <- "MRMR weights"
  return(list(weight.mr=weight.mr ,  GeneRankedList= GeneRankedList))
}
