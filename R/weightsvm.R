#Selection of top ranked (differentially expressed) genes through Support vector machine method.
#'
#' @name weightsvm
#' @aliases weightsvm
#' @title Computation of weights for gene selection through SVM algorithm
#' @usage weightsvm(x, y)
#'
#' @param x Nxp data frame of gene expression values, where, N represents number of genes and p represents samples/time points generated in a case vs. control gene expression study.
#' @param y px1 numeric vector with entries 1 and -1 representing sample/subject labels, where 1 and -1 represents the labels of subjects/ samples for case and control condition respectively.
#'
#' @return A vector of weights computed through SVM method and a complete ranked gene list.
#'
#' @description computes weights for features (e.g. genes) using Support Vector Machine (SVM) method and also produce the ranked genelist based on the computed weights.
#'
#' @references Guyon et al. (2002). Gene selection for cancer classification using support vector machines. Mach. Learn., 46, 389-422

#' @author Samarendra Das <samarendra4849 at gmail.com>
#'
#' @importFrom e1071 svm
#'
#' @examples
#' x=as.data.frame(matrix(runif(1000), 50))
#' row.names(x) = paste("Gene", 1:50)
#' colnames(x) = paste("Samp", 1:20)
#' y=as.numeric(c(rep(1, 10), rep(-1, 10)))
#' weightsvm(x, y)
#'
#' @export
#'
weightsvm <- function (x, y)
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

  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package \"e1071\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  y <- as.factor(y)
  genes <- rownames(x)
  x1 <- as.matrix(x)
  nn <- nrow(x1)
  M <- ncol(x1)
  GeneRankedList <- vector(length = nn)
  Indexes=seq(1:nn)
  #usethis::use_package("e1071")
  SvmModel <- svm(t(x1[Indexes, ]), y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
  a <- SvmModel$coefs
  b <- SvmModel$SV
  SV.weight <- abs((t(a)%*%b)) # absolute value of SVM weights
  GeneRankedList <- sort(-SV.weight, index.return = TRUE)$ix
  names(SV.weight) <- genes
  #rankvalue <- sort(GeneRankedList, index.return = TRUE)$ix
  #rankscore <- (nn + 1 - rankvalue)/(nn)

  #rankingcriteria <- as.vector(rowSums((M1), na.rm = FALSE, dims = 1))

  names(GeneRankedList) <- genes[GeneRankedList]

  #class(SV.weight) <- "SVM weights"
  return(list(SV.weight=SV.weight ,  GeneRankedList= GeneRankedList))
}
###############################weight.svm Ends here################################
