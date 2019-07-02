##' Equivalence Two-Sample t-Test
##'
##' \code{equiv.t.test} is a two-sample equivalence t-test for asymptotically normal sample means and symmetric equivalence ranges. Nonsymmetric equivalence ranges are also accepted, however p-values and confidence intervals will not be calculated.
##'
##' @param x a (non-empty) numeric vector of treatment group A.
##' @param y a (non-empty) numeric vector of treatment group B.
##' @param alpha a numeric value for the significance level, between 0 and 1 (default = 0.05).
##' @param eps_sub a numeric vector or numeric value (for a symmetric equivalence range about 0) on the scale of the variable.
##' @param eps_std a numeric vector or numeric value (for a symmetric equivalence range about 0) of the desired standard deviations to bound the equivalence region.
##' @param eps_tol a character string specifying the tolerance levels in Wellek (2010), must be either "strict" (0.36) or "liberal" (0.74).
##' @param power.out an optional logical vector to output power calculation.
##'
##' @details To set the equivalence range, enter one of the following variables for epsilon: \code{eps_sub} (preferred), \code{eps_std}, or \code{eps_tol}. Using knowledge of the data, \code{eps_sub} should be a substantive equivalence range on the scale of the variable. \code{eps_std} should specify an equivalence range of standard deviations. Otherwise, use the strict (default) and liberal epsilon tolerances, \code{eps_tol}, for the two-sample equivalence t-test as defined in Wellek (2010).
##'
##' If a nonsymmetric equivalence range is specified, p-values and equivalence confidence intervals will not be calculated.
##'
##' @return By default, \code{equivtest} returns a list of results with respect to the standard error input.
##'
##' Returns a list of class \code{equiv.t.test} with the following components:
##' \item{intype}{ epsilon input type.}
##' \item{eps_std}{ equivalence range in terms of standard deviations.}
##' \item{eps_sub}{ substantive equivalence range on the scale of the variable.}
##' \item{test.stat}{ value of the test statistic.}
##' \item{rej}{ TRUE value if null hypothesis rejected.}
##' \item{p}{ p-value from hypothesis test.}
##' \item{critical.const}{ critical interval for t-statistic.}
##' \item{SE}{ standard error for given input type.}
##' \item{CI_std}{ equivalence confidence interval in terms of standard deviations.}
##' \item{CI_sub}{ equivalence confidence interval on the scale of the variable.}
##' \item{power}{ power of the test.}
##'
##' @seealso \code{\link{equiv.paired.t.test}}
##'
##' @references
##' Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013.
##'
##' Wellek, S. (2010). Testing statistical hypotheses of equivalence and noninferiority. Chapman and Hall/CRC.
##'
##' @import plm Hmisc stringr sandwich weights
##' @examples
##' # Wellek p 124
##' x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
##' y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
##' eps=c(.5,1)
##' res=equiv.t.test(x,y,eps_sub=eps)
##' summary(res)
##' @export
equiv.t.test <- function(x, y, alpha = .05, eps_sub, eps_std, eps_tol,
                         power.out = FALSE) {

  ################ Error handling for alpha inputs
  if(alpha < 0 | alpha > 1){
    warning("Incorrect alpha specified, using default alpha = 0.05.")
    alpha=0.05
  }

  ################ Error handling for x, y inputs
  # if x, y are not numeric vectors: stop
  if(is.numeric(x)==FALSE | is.numeric(y)==FALSE) stop("x and y must be numeric vectors")
  # drop missing values from x, y
  missing_x = is.na(x)
  missing_y = is.na(y)
  x = x[!missing_x]
  y = y[!missing_y]
  # warning if NAs dropped
  if(sum(missing_x) > 0) warning(paste(sum(missing_x), 'NA value(s) dropped from x'))
  if(sum(missing_y) > 0) warning(paste(sum(missing_y), 'NA value(s) dropped from y'))
  weights.x = rep(1, length(x))
  weights.y = rep(1, length(y))

  ########### Calculations

  dbar <- wtd.mean(x, weights = weights.x) - wtd.mean(y, weights = weights.y) # difference in means
  m <- as.double(length(x)) # no. observations x
  n <- as.double(length(y)) # no. observations y
  N <- m+n                  # total observations
  x.var <-  wtd.var(x, weights = weights.x) # variance x
  y.var <- wtd.var(y, weights = weights.y)  # variance y
  obs_smd = (wtd.mean(x, weights = weights.x) - wtd.mean(y, weights = weights.y)) / sqrt(y.var)
  pooled.sd = sqrt(((m-1)*x.var+(n-1)*y.var)/(m+n-2))
  pooled.se = sqrt((m-1)*x.var + (n-1)*y.var) / sqrt(m*n*(N-2)/N)


  # calculate SE, df
  df = N - 2
  SE = pooled.se
  pooled.sd = SE * sqrt((m*n*(N-2)/N)/(m+n-2))

  ################ Error handling for epsilon inputs
  # if more than one eps is specified: stop
  if(sum(!missing(eps_std), !missing(eps_sub), !missing(eps_tol)) > 1) {
    stop("Enter at most one variable for epsilon: eps_sub, eps_std, or eps_tol.")
  }
  # if no eps is specified: set to strict tolerance
  if(missing(eps_sub) == TRUE & missing(eps_std)==TRUE & missing(eps_tol)==TRUE) {
    warning("No eps_sub, eps_std, or eps_tol specified, using a default value of eps_tol = strict (.36).")
    this="eps_tol = strict"
    eps_tol = "strict"
  }

  ################ Setting the epsilon value
  if(!missing(eps_sub)){
    epsilon=eps_sub/pooled.sd
    this="eps_sub"
  }
  if(!missing(eps_std)) {
    epsilon=eps_std
    this="eps_std"
  }
  if(!missing(eps_tol)){
    if(eps_tol == "liberal") {
      epsilon=.74*pooled.sd
      this="eps_tol = liberal"
    } else if(eps_tol=="strict") {
      epsilon=.36*pooled.sd
      this="eps_tol = strict"
      warning("This strict (default) method is very conservative, use of eps_sub or eps_std is recommended.")
    } else{
      stop("eps_tol accepts inputs 'liberal' or 'strict'.")
    }
  }

  if(length(epsilon)>2 | length(epsilon)<1) stop("Epsilon must be an equivalence range or the bound of a symmetric equivalence range about zero.")
  if(length(epsilon)==1) epsilon=c(-abs(epsilon),abs(epsilon))

  # calculate critical constants
  critconst.pow = invisible(capture.output(tt2st(m,n,alpha,abs(epsilon[1]),epsilon[2],tol=1e-6,itmax = 30)))
  critconst2=c(str_extract(critconst.pow,pattern="c1 = -?\\d+\\.*\\d*"),str_extract(critconst.pow,pattern="c2 = -?\\d+\\.*\\d*"))
  critconst3=as.numeric(c(str_extract(critconst2[1]," -?\\d+\\.*\\d*"),str_extract(critconst2[2]," -?\\d+\\.*\\d*")))

  # calculate t statistic, rejection region
  t.stat <-  dbar / SE
  if(abs(epsilon[1])==abs(epsilon[2])) rej = (abs(t.stat) >= critconst3[1] & abs(t.stat) <= critconst3[2])
  if(abs(epsilon[1]) != abs(epsilon[2])) rej = (t.stat >= critconst3[1] & t.stat <= critconst3[2])

  # calculate inverted, p
  if(abs(epsilon[1])==abs(epsilon[2])) {
     non.cent <- (m*n*epsilon[1]^2)/N
     p = pf(abs(t.stat)^2, 1, df, non.cent)
     inverted <- try(uniroot(function(x) pf(abs(t.stat)^2, 1, N-2, ncp = (m*n*x^2)/N) -
                               ifelse(pf(abs(t.stat)^2, 1, N-2, ncp = (m*n*0^2)/N) < alpha,
                                      pf(abs(t.stat)^2, 1, N-2, ncp = (m*n*obs_smd^2)/N), alpha),
                             c(0,10*abs(t.stat)), tol = 0.0001)$root, silent = TRUE)
     if(class(inverted) == "try-error") inverted = NA
     CI_std = c(-inverted, inverted)
     CI_sub = c(-inverted*pooled.sd,inverted*pooled.sd)
   } else {
     # non.cent1 <- (m*n*epsilon[1]^2)/N
     # non.cent2 <- (m*n*epsilon[2]^2)/N
     p = NA
     CI_std = NA
     CI_sub = NA
  }

  # convert results to appropriate epsilon format and alternate format
  eps_std = epsilon
  eps_sub = epsilon*pooled.sd

  # calculate the power if necessary
  if(power.out==TRUE){
    power=c(str_extract(critconst.pow,pattern="POW0 = \\d+\\.*\\d*"))
    power2=as.numeric(str_extract(power," \\d+\\.*\\d*"))
  } else {
    power2=NA
  }

  # return 'equivtest' class result for summary function
  result <- list(intype=this,
                 eps_std=eps_std,
                 eps_sub=eps_sub,
                 test.stat = t.stat,
                 rej = rej,
                 p = p,
                 critical.const=critconst3,
                 SE=SE,
                 CI_std=CI_std,
                 CI_sub=CI_sub,
                 power=power2)

  class(result) <- "equivttest"

  return(result)
}


# Example 6.1, Wellek 2010 p.124, given epsilon of standardized difference in means
# x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
# y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
# eps=c(-.5,1)
# equiv.t.test(x,y,alpha=.05, eps_std = eps)
# m <- as.double(length(x)) # no. observations x
# n <- as.double(length(y)) # no. observations y
# N <- m+n                  # total observations
# x.var <-  var(x) # variance x
# y.var <- var(y)  # variance y
# pooled.sd = sqrt(((m-1)*x.var+(n-1)*y.var)/(m+n-2))
# res1 = equiv.t.test(x,y,eps_std=eps)
# res2 = equiv.t.test(x,y,eps_sub=eps*pooled.sd)
#
# library(EQUIVNONINF)
# tt2st(m=length(x),n=length(y), alpha=0.05,eps1=.5,eps2=1,tol=1e-6,itmax=100)
