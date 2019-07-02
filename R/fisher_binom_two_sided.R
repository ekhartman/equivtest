##' Exact Fisher Binomial Two-Sided Test for Two Samples
##'
##' \code{fisher.binom.two.sided} is a two-sided equivalence Fisher type exact test for two binomial distributions with respect to the odds ratio. The alternative hypothesis is: the odds ratio is contained within a given equivalence range.
##'
##' @param x a (non-empty) numeric vector of binary responses from treatment group A.
##' @param y a (non-empty) numeric vector of binary responses from treatment group B.
##' @param eps_std a numeric value or vector to determine the bounds of the equivalence odds ratio.
##' @param eps_sub a numeric value or vector representing the probabilities that are substantively considered to demonstrate nonequivalence between the two groups.
##' @param eps_tol a character string specifying the tolerance levels in Wellek (2010), must be either "strict" (0.41) or "liberal" (0.85).
##' @param alpha a numeric value for the significance level, between 0 and 1 (default = 0.05).
##' @param power.out an optional logical vector to output power calculation.
##'
##' @return Returns a list with the following components:
##' \item{intype}{ the type of epsilon input from user.}
##' \item{odds.ratio}{ odds ratio of treatment groups A and B.}
##' \item{p}{ p-value of the test.}
##' \item{rej}{ 1 if null hypothesis is rejected.}
##' \item{equiv.CI}{ equivalence confidence interval.}
##'
##' @details To set the equivalence range, enter one of the following variables for epsilon: \code{eps_sub} (preferred), \code{eps_std}, or \code{eps_tol}. Using knowledge of the data, \code{eps_sub} should be the probability that is substantively considered to demonstrate nonequivalence between the two groups (Wellek 2010, p. 16). \code{eps_std} bounds the odds ratio by epsilon. Otherwise, use the strict (default) and liberal epsilon tolerances, \code{eps_tol}, for the two-sided Fisher-Binomial test as defined in Wellek (2010).
##'
##' @seealso \code{\link{fisher.binom.one.sided}}
##'
##' @references
##' Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013.
##'
##' Wellek, S. (2010). Testing statistical hypotheses of equivalence and noninferiority. Chapman and Hall/CRC.
##'
##' @examples
##' # Two-sided example from Wellek 2010, pg 191
##' A = c(rep(1, 108), rep(0, 117))
##' B = c(rep(1, 63), rep(0, 56))
##' fisher.binom.two.sided(A, B, eps_std=.6667, alpha=0.05)
##' @import BiasedUrn EQUIVNONINF
##' @export
fisher.binom.two.sided = function(x, y, eps_std, eps_tol, eps_sub, alpha, power.out=FALSE) {

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

  if(length(epsilon) != 1 && length(epsilon) != 2)
    stop("Epsilon must be an equivalence range or the bounds of an equivalence range about the odds ratio.")

  if(length(epsilon)==1) {
    epsilon <- ifelse(epsilon < 1, epsilon, 1/epsilon)
    epsilon <- c(epsilon, 1/epsilon)
  }


  p <- max(
    # prob greater than observed with null lower ratio
          sum(dFNCHypergeo(x.sum:m, m, n, s, epsilon[1]))
    # prob less than observed with upper ratio
          , sum(dFNCHypergeo(0:(x.sum), m, n, s, epsilon[2]))
      )

  rej <- ifelse(p < alpha,TRUE,FALSE)

  critical.vals.x <- c(qFNCHypergeo(alpha, m, n, s, epsilon[1], lower.tail = FALSE)
                    , qFNCHypergeo(alpha, m, n, s, epsilon[2]))
  critical.vals.y <- c(qFNCHypergeo(alpha, n, m, s, epsilon[1], lower.tail = FALSE)
                    , qFNCHypergeo(alpha, n, m, s, epsilon[2]))
  critical.vals <- matrix(c(critical.vals.x, critical.vals.y), ncol = 2
                         , dimnames = list(c("lower", "upper"),
                                           c("X", "Y")))

  inverted <- try(uniroot(function(x) max(sum(dFNCHypergeo(x.sum:m, m, n, s, x))
                                         , sum(dFNCHypergeo(0:(x.sum), m, n, s, 1/x))
                                      ) - alpha, c(0.00001, 10), tol = 0.0001)$root)
  if(class(inverted) == "try-error") {
    inverted <- NA
  }

  # calculate odds ratio range
  crit.p1.l <- critical.vals.x[1] / m
  crit.p1.u <- critical.vals.x[2] / m
  crit.p2.l <- critical.vals.y[1] / n
  crit.p2.u <- critical.vals.y[2] / n

  if(power.out==TRUE){
    power <- invisible(capture.output(bi2aeq1(m,n,rho1=epsilon[1],rho2=epsilon[2],alpha,p1,p2)))
    power <- power[3]
    power <- as.numeric(sub(".*POW= *", "\\1", power))
  } else{
    power <- NA
  }

  critical.vals <- c(max(odds.ratio(crit.p1.l, crit.p2.u), odds.ratio(crit.p2.l, crit.p1.u))
                      , min(odds.ratio(crit.p1.u, crit.p2.l), odds.ratio(crit.p2.u, crit.p1.l)))

  # returns an odds ratio, which isn't really the critical value...

  result = list(intype = intype,
                odds.ratio = odds,
                p = p,
                rej = rej,
                equiv.CI = c(inverted, 1 / inverted),
                power = power)
  class(result) = "equivfisherbinom"

  return(result)
}


