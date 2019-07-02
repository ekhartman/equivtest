##' Two-One-Sided Ratio Test (TOST Ratio)
##'
##' \code{tost.ratio} is a two-one-sided equivalence ratio test for the ratio of the population means of the two treatment groups. The alternative hypothesis tests whether the ratio of the treatment to control population means are within a given equivalence range.
##'
##' @param x a (non-empty) numeric vector of treatment group A.
##' @param y a (non-empty) numeric vector of treatment group B.
##' @param frac a numeric value for the ratio of the population mean of treatment group A to the population mean of treatment group B.
##' @param eps_sub a numeric vector or numeric value (for a symmetric equivalence range about 0) on the scale of the variable.
##' @param eps_std a numeric vector or numeric value (for a symmetric equivalence range about 0) of the desired standard mean differences on the outcome using the pooled standard deviation.
##' @param alpha a numeric value for the confidence level, between 0 and 1 (default = 0.05).
##'
##' @return
##' \item{rej}{ TRUE value if null hypothesis rejected.}
##' \item{tstat}{ the lower and upper test statistics. }
##'
##' @references
##' Berger, R. L., & Hsu, J. C. (1996). Bioequivalence trials, intersection-union tests and equivalence confidence sets. Statistical Science, 11(4), 283-319.
##'
##' Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013.
##'
##' @examples tost.ratio()
##' @export
tost.ratio = function(x, y, frac, alpha = 0.05,  eps_sub, eps_std) {

  ################ Error handling for alpha inputs
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    warning("Incorrect alpha specified, using default alpha of 0.05.")
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

  dbar <- mean(x) - mean(y) # difference in means
  m <- as.double(length(x)) # no. observations x
  n <- as.double(length(y)) # no. observations y
  N <- m+n                  # total observations
  x.var <-  var(x) # variance x
  y.var <- var(y)  # variance y
  pooled.sd = sqrt(((m-1)*x.var+(n-1)*y.var)/(m+n-2))
  pooled.se = sqrt((m-1)*x.var + (n-1)*y.var) / sqrt(m*n * (N-2)/N)
  obs_smd = dbar / pooled.sd

  ################ Error handling for epsilon inputs
  # if more than one eps is specified: stop
  if(sum(!missing(eps_std), !missing(eps_sub)) > 1) {
    stop("Enter at most one variable for epsilon: eps_sub or eps_std.")
  }
  # if no eps is specified: set to strict tolerance
  if(missing(eps_sub) == TRUE & missing(eps_std)==TRUE) {
    warning("No eps_sub or eps_std specified, using a default value of eps_std = 1.")
    this="eps_std"
    eps_std = obs_smd
  }

  ################ Setting the epsilon value
  if(!missing(eps_sub)){
    epsilon=eps_sub
    this="eps_sub"
  }
  if(!missing(eps_std)) {
    epsilon=eps_std*pooled.sd
    this="eps_std"
  }

  if(length(epsilon)>2 | length(epsilon)<1) stop("Epsilon must be an equivalence range or the bound of a symmetric equivalence range about zero.")
  if(length(epsilon)==1) epsilon=c(-abs(epsilon), abs(epsilon))

	if(length(unique(x)) == 2) {
		delta = frac
	} else {
		delta = frac * min(sqrt(var(x)), sqrt(var(y)))
	}

	# t.u = t.test(x, y, mu = delta, alternative = "less")$statistic
  t.u = (mean(y) - mean(x)*epsilon[2])/(pooled.se * sqrt(1/m + epsilon[2]^2 /n))
	t.u.df = t.test(x, y, mu = delta, alternative = "less")$parameter
	# t.l = t.test(x, y, mu = -delta, alternative = "greater")$statistic
	t.l = (mean(y) - mean(x)*epsilon[1])/(pooled.se * sqrt(1/m + epsilon[1]^2 /n))

	t.l.df = t.test(x, y, mu = -delta, alternative = "greater")$parameter
	rej = as.logical(ifelse((t.u < -1 * pt(alpha, t.u.df)) & (t.l > pt(alpha, t.l.df)), TRUE, FALSE))

	result = list(rej = rej, tstat = c(t.u, t.l))
	class(result) <- "equivtostratio"

	return(result)
}


