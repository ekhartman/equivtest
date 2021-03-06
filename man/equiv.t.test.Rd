% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv_t_test.R
\name{equiv.t.test}
\alias{equiv.t.test}
\title{Equivalence Two-Sample t-Test}
\usage{
equiv.t.test(x, y, alpha = 0.05, eps_sub, eps_std, eps_tol,
  power.out = FALSE)
}
\arguments{
\item{x}{a (non-empty) numeric vector of treatment group A.}

\item{y}{a (non-empty) numeric vector of treatment group B.}

\item{alpha}{a numeric value for the significance level, between 0 and 1 (default = 0.05).}

\item{eps_sub}{a numeric vector or numeric value (for a symmetric equivalence range about 0) on the scale of the variable.}

\item{eps_std}{a numeric vector or numeric value (for a symmetric equivalence range about 0) of the desired standard deviations to bound the equivalence region.}

\item{eps_tol}{a character string specifying the tolerance levels in Wellek (2010), must be either "strict" (0.36) or "liberal" (0.74).}

\item{power.out}{an optional logical vector to output power calculation.}
}
\value{
By default, \code{equivtest} returns a list of results with respect to the standard error input.

Returns a list of class \code{equiv.t.test} with the following components:
\item{intype}{ epsilon input type.}
\item{eps_std}{ equivalence range in terms of standard deviations.}
\item{eps_sub}{ substantive equivalence range on the scale of the variable.}
\item{test.stat}{ value of the test statistic.}
\item{rej}{ TRUE value if null hypothesis rejected.}
\item{p}{ p-value from hypothesis test.}
\item{critical.const}{ critical interval for t-statistic.}
\item{SE}{ standard error for given input type.}
\item{CI_std}{ equivalence confidence interval in terms of standard deviations.}
\item{CI_sub}{ equivalence confidence interval on the scale of the variable.}
\item{power}{ power of the test.}
}
\description{
\code{equiv.t.test} is a two-sample equivalence t-test for asymptotically normal sample means and symmetric equivalence ranges. Nonsymmetric equivalence ranges are also accepted, however p-values and confidence intervals will not be calculated.
}
\details{
To set the equivalence range, enter one of the following variables for epsilon: \code{eps_sub} (preferred), \code{eps_std}, or \code{eps_tol}. Using knowledge of the data, \code{eps_sub} should be a substantive equivalence range on the scale of the variable. \code{eps_std} should specify an equivalence range of standard deviations. Otherwise, use the strict (default) and liberal epsilon tolerances, \code{eps_tol}, for the two-sample equivalence t-test as defined in Wellek (2010).

If a nonsymmetric equivalence range is specified, p-values and equivalence confidence intervals will not be calculated.
}
\examples{
# Wellek p 124
x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
eps=c(.5,1)
res=equiv.t.test(x,y,eps_sub=eps)
summary(res)
}
\references{
Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013.

Wellek, S. (2010). Testing statistical hypotheses of equivalence and noninferiority. Chapman and Hall/CRC.
}
\seealso{
\code{\link{equiv.paired.t.test}}
}
