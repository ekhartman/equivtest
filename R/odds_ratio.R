##' Odds Ratio
##'
##' The odds ratio is an appropriate dissimilarity measure for two binomial samples (Wellek 2010). Therefore, this result helps determines our equivalence hypotheses for the \code{\link{fisher.binom.one.sided}} and \code{\link{fisher.binom.two.sided}} tests.
##'
##' @param p1 a (nonempty) numeric value of the average of the binary responses from treatment group A.
##' @param p2 a (nonempty) numeric value of the average of the binary responses from treatment group B.
##'
##' @return Calculates the odds ratio using the outcomes from treatment groups A and B.
##'
##' @seealso \code{\link{fisher.binom.one.sided}}, \code{\link{fisher.binom.two.sided}}
##' @examples
##' # One-sided example (from Wellek 2010 pg 176)
##' A = c(rep(1, 98), rep(0, 8))
##' B = c(rep(1, 97), rep(0, 10))
##' odds.ratio(mean(A),mean(B))
odds.ratio = function(p1, p2) {
  return((p1 * ( 1 - p2)) /( (1 - p1) * p2 ))
}
