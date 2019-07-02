##' Equivalence Fisher-Binomial Test Summary Function
##'
##' \code{summary.equivfisherbinom} is a summary function for objects of the class \code{equivfisherbinom}. This summary function should be paired with the \code{fisher.binom.one.sided} or \code{fisher.binom.two.sided} functions.
##'
##' @seealso \code{\link{fisher.binom.two.sided}}, \code{\link{fisher.binom.one.sided}}
##'
##' @examples
##' # One-sided example from Wellek 2010, pg 176
##' A = c(rep(1, 98), rep(0, 8))
##' B = c(rep(1, 97), rep(0, 10))
##' res = fisher.binom.one.sided(A, B, eps_std=.5, alpha=0.05)
##' summary(res)
##' @export
summary.equivfisherbinom <- function(res){

  # one-sided
  if(length(res)==5){
    cat("Equivalence One-Sided Fisher-Binomial Test \n")
    cat(paste("Input: ", res$intype),"\n")
    cat(paste("Odds Ratio of Treatment Groups A and B: ", round(res$odds.ratio,3)),"\n")
    cat(paste0("Rejected the null? ", paste0(ifelse(res$rej==1,"TRUE","FALSE")), ", p-value of ", round(res$p,3), "\n"))
    cat(paste("Power of the test: ", round(res$power,3)),"\n")

    # two-sided
  } else{
    cat("Equivalence Two-Sided Fisher-Binomial Test \n")
    cat(paste("Input: ", res$intype),"\n")
    cat(paste("Odds Ratio of Treatment Groups A and B: ", round(res$odds.ratio,3)),"\n")
    cat(paste("Equivalence Odds Ratio CI: ", round(res$equiv.CI[1],3)),"to",round(res$equiv.CI[2],3),"\n")
    cat(paste0("Rejected the null? ", paste0(ifelse(res$rej==1,"TRUE","FALSE")), ", p-value of ", round(res$p,3), "\n"))
    cat(paste("Power of the test: ", round(res$power,3)),"\n")

  }
}





