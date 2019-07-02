##' Equivalence Paired t-Test Summary Function
##'
##' \code{summary.equivpairedt} is a summary function for objects of the class \code{equivpairedt}. This summary function should be paired with the \code{equiv.paired.t.test} function.
##'
##' @seealso \code{\link{equiv.paired.t.test}}
##'
##' @examples
##' # Wellek p 96
##' res = equiv.paired.t.test(diff.mu=.16, diff.sd= 3.99, n=23, alpha = .05, eps_std=.5)
##' summary(res)
##' @export
summary.equivpairedt <- function(res){
  cat("Equivalence paired t-test \n")
  cat(paste("Input type: ", res[[1]][1]),"\n")
  cat(paste("t-test statistic: ", round(res[[4]][1],3)),"\n")
  cat(paste("Equivalence CI: ", round(res[[7]][1],3)),"to",round(res[[7]][2],3),"\n")
  cat(paste("Rejected the null?", ifelse(res[[5]][1]==1,"TRUE","FALSE")),"\n")
  cat(paste("Power of the test: ", round(res[[9]][1],3)),"\n")
}
