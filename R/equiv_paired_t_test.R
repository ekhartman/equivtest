##' Equivalence Paired t-Test
##'
##' \code{equiv.paired.t.test} is a paired equivalence t-test for an asymptotically normal vector of individual differences taken from two treatment groups for symmetric equivalence ranges. Nonsymmetric equivalence ranges are also accepted, however p-values and confidence intervals will not be calculated at this time.
##'
##' @param diff a numeric vector of the differences between paired observations from two treatment groups.
##' @param alpha a numeric value for the significance level, between 0 and 1 (default = 0.05).
##' @param eps_sub a numeric vector or numeric value (for a symmetric equivalence range about 0) on the scale of the variable.
##' @param eps_std a numeric vector or numeric value (for a symmetric equivalence range about 0) of the desired standard deviations to bound the equivalence region.
##' @param eps_tol a character string specifying the tolerance levels in Wellek (2010), must be either "strict" (0.25) or "liberal" (0.50).
##' @param power.out an optional logical vector to output power calculation.
##'
##' @details To set the equivalence range, enter one of the following variables for epsilon: \code{eps_sub} (preferred), \code{eps_std}, or \code{eps_tol}. Using knowledge of the data, \code{eps_sub} should be a substantive equivalence range on the scale of the variable. \code{eps_std} should specify an equivalence range of standard deviations. Otherwise, use the strict (default) and liberal epsilon tolerances, \code{eps_tol}, for the paired equivalence t-test as defined in Wellek (2010).
##'
##' @return Returns a list with the following components:
##' \item{intype}{ the type of epsilon input from user.}
##' \item{t.stat}{ value of the t-statistic.}
##' \item{rej}{ TRUE value if null hypothesis rejected.}
##' \item{CI_std}{ equivalence confidence interval in terms of standard deviations.}
##' \item{CI_sub}{ equivalence confidence interval on the scale of the variable.}
##' \item{power}{ power of the test}
##'
##' @seealso \code{\link{equiv.t.test}}
##'
##' @references
##' Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013.
##'
##' Wellek, S. (2010). Testing statistical hypotheses of equivalence and noninferiority. Chapman and Hall/CRC.
##'
##' @examples
##' diff = c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
##' equiv.paired.t.test(diff, alpha = .05, eps_std=.5)
##' @export
equiv.paired.t.test <- function(diff, alpha = .05, eps_sub, eps_std, eps_tol, power.out=FALSE) {

  ################ Error handling

  if(!is.numeric(diff))
    stop("diff must be a numeric vector")
  if(anyNA(diff))
    stop("diff must not contain NA value(s).")

  if(sum(!missing(eps_sub), !missing(eps_std), !missing(eps_tol)) > 1) {
    stop("Enter at most one variable for epsilon: eps_sub, eps_std, or eps_tol.")
  } else if (missing(eps_sub) && missing(eps_std) && missing(eps_tol)) {
    warning("No eps_sub, eps_std, or eps_tol specified, using a default value of eps_tol = strict (.25).")
    eps_tol <- "strict"
  }

  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    warning("Incorrect alpha specified, using default alpha of 0.05.")
    alpha=0.05
  }



  ################ Calculations

  # Setting the epsilon value
  if(!missing(eps_sub)){
    epsilon <- eps_sub/diff.sd
    intype <- "eps_sub"
  } else if(!missing(eps_std)) {
    epsilon <- eps_std
    intype <- "eps_std"
  } else if(!missing(eps_tol)){
    if(eps_tol == "liberal") {
      epsilon <- .5*diff.sd
      intype <- "eps_tol = liberal"
    } else if(eps_tol=="strict") {
      epsilon <- .25*diff.sd
      intype <- "eps_tol = strict"
      warning("This strict (default) method is very conservative, use of eps_sub or eps_std is recommended.")
    } else{
      stop("eps_tol accepts inputs 'liberal' or 'strict'.")
    }
  }

  if(length(epsilon)==1)
    epsilon=c(-abs(epsilon),abs(epsilon))
  else if(length(epsilon) != 2)
    stop("Epsilon must be an equivalence range or the bound of a symmetric equivalence range about zero.")

  diff.mu <- mean(diff)
  diff.sd <- sd(diff)
  n <- length(diff)

  critconst.pow <- invisible(capture.output(tt1st(n,alpha,theta1=epsilon[1],theta2=epsilon[2], tol=1e-6, itmax=100)))
  critconst2 <- c(str_extract(critconst.pow,pattern="c1 = -?\\d+\\.*\\d*"),str_extract(critconst.pow,pattern="c2 = -?\\d+\\.*\\d*"))
  critconst3 <- as.numeric(c(str_extract(critconst2[1]," -?\\d+\\.*\\d*"),str_extract(critconst2[2]," -?\\d+\\.*\\d*")))

  if(abs(epsilon[1])==abs(epsilon[2])) {
    t.stat <- diff.mu/(diff.sd/sqrt(n))
    rej <- (abs(t.stat) >= critconst3[1] & abs(t.stat) <= critconst3[2])
  } else if (abs(epsilon[1]) != abs(epsilon[2])) {
    t.stat <- diff.mu/(diff.sd/sqrt(n))
    rej <- (t.stat >= critconst3[1] & t.stat <= critconst3[2])
  }

  # only calculate inverted, p for nonsymmetric equivalence ranges
  if(abs(epsilon[1])==abs(epsilon[2])) {
    inverted <- uniroot(function(x) pf(abs(t.stat)^2, 1, n-1, ncp = (n*x^2)) - alpha, c(0,2*abs(t.stat)), tol = 0.0001)$root  #invert test-statistic
    if(class(inverted) == "try-error") inverted = NA
    CI_std <- c(-inverted, inverted)
    CI_sub <- c(-inverted*diff.sd,inverted*diff.sd)
  } else {
    CI_std <- NA
    CI_sub <- NA
  }
  eps_std <- epsilon
  eps_sub <- epsilon*diff.sd

  if(power.out==TRUE){
    power <- c(str_extract(critconst.pow,pattern="POW0 = \\d+\\.*\\d*"))
    power2 <- as.numeric(str_extract(power," \\d+\\.*\\d*"))
  } else {
    power2 <- NA
  }

  result <- list(intype=intype,
                 eps_std=eps_std,
                 eps_sub=eps_sub,
                 test.stat = t.stat,
                 rej = rej,
                 critical.const=critconst3,
                 CI_std=CI_std,
                 CI_sub=CI_sub,
                 power=power2)

  class(result) = "equivpairedt"
  return(result)
}
