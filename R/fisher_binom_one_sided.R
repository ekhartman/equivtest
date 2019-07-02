##' Exact Fisher Binomial One-Sided Test for Two Samples
##'
##' \code{fisher.binom.one.sided} is a one-sided equivalence Fisher type exact test for two binomial distributions with respect to the odds ratio. The alternative hypothesis is: the odds ratio is less than or equal to 1- \code{eps_std}.
##'
##' @param x a (non-empty) numeric vector of binary responses from treatment group A.
##' @param y a (non-empty) numeric vector of binary responses from treatment group B.
##' @param eps_std a numeric value for epsilon to determine the 1-epsilon hypothesis.
##' @param eps_sub a numeric value representing the probability that is substantively considered to demonstrate nonequivalence between the two groups.
##' @param eps_tol a character string specifying the tolerance levels in Wellek (2010), must be either "strict" (0.41) or "liberal" (0.85).
##' @param alpha a numeric value for the significance level, between 0 and 1 (default = 0.05).
##' @param power.out an optional logical vector to output the power calculation.
##'
##' @return Returns a list with the following components:
##' \item{intype}{ the type of epsilon input from user.}
##' \item{odds.ratio}{ odds ratio of treatment groups A and B.}
##' \item{p}{ p-value of the test.}
##' \item{rej}{ 1 if null hypothesis is rejected.}
##'
##' @details To set the equivalence range, enter one of the following variables for epsilon: \code{eps_sub} (preferred), \code{eps_std}, or \code{eps_tol}. Using knowledge of the data, \code{eps_sub} should be the probability that is substantively considered to demonstrate nonequivalence between the two groups (Wellek 2010, p. 16). \code{eps_std} creates an upper bound of 1-epsilon for the odds ratio. Otherwise, use the strict (default) and liberal epsilon tolerances, \code{eps_tol}, for the one-sided Fisher-Binomial test as defined in Wellek (2010).
##'
##' @seealso \code{\link{fisher.binom.two.sided}}
##'
##' @references
##' Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013.
##'
##' Wellek, S. (2010). Testing statistical hypotheses of equivalence and noninferiority. Chapman and Hall/CRC.
##'
##' @examples
##' # One-sided example from Wellek 2010, pg 176
##' A = c(rep(1, 98), rep(0, 8))
##' B = c(rep(1, 97), rep(0, 10))
##' fisher.binom.one.sided(A, B, eps_std=.5, alpha=0.05)
##' @import EQUIVNONINF
##' @export
fisher.binom.one.sided = function(x, y, eps_std, eps_sub, eps_tol, alpha, power.out=FALSE) {

  ################ Error handling

  if(!is.numeric(x) || !is.numeric(y))
    stop("x and y must be numeric vectors.")
  if(anyNA(x) || anyNA(y))
    stop("x and y must not contain NA value(s).")

  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    warning("Incorrect alpha specified, using default alpha of 0.05.")
    alpha=0.05
  }

  if(sum(!missing(eps_sub), !missing(eps_std), !missing(eps_tol)) > 1) {
    stop("Enter at most one variable for epsilon: eps_sub, eps_std, or eps_tol.")
  } else if (missing(eps_sub) && missing(eps_std) && missing(eps_tol)) {
    warning("No eps_sub, eps_std, or eps_tol specified, using a default value of eps_tol = strict (.41).")
    eps_tol <- "strict"
  }

  x.sum <- sum(x == 1)
  m <- length(x)
  p1 <- x.sum / m
  y.sum <- sum(y == 1)
  n <- length(y)
  p2 <- y.sum / n
  odds <- odds.ratio(p1, p2)
  s <- x.sum + y.sum


  # Setting the epsilon value
  if(!missing(eps_sub)){
    epsilon <- log((1+2*eps_sub)/(1-2*eps_sub))
    intype <- "eps_sub"
  } else if(!missing(eps_std)) {
    epsilon <- eps_std
    intype <- "eps_std"
  } else if(!missing(eps_tol)){
    if(eps_tol == "liberal") {
      epsilon <- .85*pooled.sd
      intype <- "eps_tol = liberal"
    } else if(eps_tol=="strict") {
      epsilon <- .41*pooled.sd
      intype <- "eps_tol = strict"
      warning("This strict (default) method is very conservative, use of eps_sub or eps_std is recommended.")
    } else{
      stop("eps_tol accepts inputs 'liberal' or 'strict'.")
    }
  }

  if(length(epsilon) != 1)
    stop("Must enter one epsilon value between 0 and 1.")

  # test results
  p <- fisher.test(matrix(c(x.sum, y.sum, m - x.sum, n - y.sum),nrow = 2,
                         dimnames =list(c("A", "B"),c("Positive", "Negative"))), or = 1 - epsilon, alternative = "greater")$p
  rej <- ifelse(p < alpha,TRUE,FALSE)

  if(power.out==TRUE){
    power <- invisible(capture.output(bi2ste1(m, n, eps=epsilon, alpha, p1, p2)))
    power <- power[4]
    power <- as.numeric(sub(".* POW = *", "\\1", power))
  } else {
    power <- NA
  }

  result = list(intype = intype,
                odds.ratio = odds,
                p = p,
                rej = rej,
                power= power)
  class(result) = "equivfisherbinom"
  return(result)
}

# One-sided example from Wellek 2010, pg 176
# A = c(rep(1, 98), rep(0, 8))
# B = c(rep(1, 97), rep(0, 10))
# fisher.binom.one.sided(A, B, eps_std=.5, alpha=0.05)
# bi2ste1(length(A), length(B), eps=.5, alpha=0.05, sum(A)/length(A), sum(B)/length(B))
