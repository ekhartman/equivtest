context("Wellek book examples")

library(equivtest)
library(EQUIVNONINF)
library(stringr)

test_that("t test Wellek example works", {

  # Example 6.1, Wellek 2010 p.124, given epsilon of standardized difference in means
  x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
  y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
  eps=c(.5,1)
  m <- as.double(length(x)) # no. observations x
  n <- as.double(length(y)) # no. observations y
  N <- m+n                  # total observations
  x.var <-  var(x) # variance x
  y.var <- var(y)  # variance y
  pooled.sd = sqrt(((m-1)*x.var+(n-1)*y.var)/(m+n-2))
  res1 = equiv.t.test(x,y,eps_std=eps, power.out = TRUE)
  res2 = equiv.t.test(x,y,eps_sub=eps*pooled.sd, power.out = TRUE)

  # check works for both eps_std and eps_sub inputs
  expect_equal(res1$eps_std, res2$eps_std)
  expect_equal(res1$eps_sub, res2$eps_sub)
  expect_equal(res1$test.stat, res2$test.stat)
  expect_equal(res1$rej, res2$rej)
  expect_equal(res1$critical.const, res2$critical.const)
  expect_equal(res1$SE, res2$SE)
  expect_equal(res1$power, res2$power)

  # check against Wellek function
  res3 = invisible(capture.output(tt2st(m=length(x),n=length(y), alpha=0.05,eps1=.5,eps2=1,tol=1e-6,itmax=100)))
  eps = c(str_extract(res3,pattern= "eps1 = -?\\d+\\.*\\d*"),str_extract(res3,pattern= "eps2 = -?\\d+\\.*\\d*"))
  eps=c(as.numeric(str_extract(eps[1]," -?\\d+\\.*\\d*")),as.numeric(str_extract(eps[2]," -?\\d+\\.*\\d*")))
  cc = c(str_extract(res3,pattern= "c1 = -?\\d+\\.*\\d*"),str_extract(res3,pattern= "c2 = -?\\d+\\.*\\d*"))
  cc=c(as.numeric(str_extract(cc[1]," -?\\d+\\.*\\d*")),as.numeric(str_extract(cc[2]," -?\\d+\\.*\\d*")))
  pow = c(str_extract(res3,pattern= " POW0 = -?\\d+\\.*\\d*"))
  pow=c(as.numeric(str_extract(pow[1]," -?\\d+\\.*\\d*")))

  # check works for Wellek book function
  expect_equal(res1$eps_std, eps, tolerance=1e-8)
  expect_equal(res1$critical.const, cc, tolerance=1e-8)
  expect_equal(res1$power, pow, tolerance=1e-8)

})


test_that("paired t test Wellek example works", {

  # Wellek p 96, given epsilon of standardized expected change of variable
  res3= invisible(capture.output(tt1st(n=23, alpha=0.05, theta1= -.5, theta2=.5, tol=1e-6, itmax=100)))

  diff = rnorm(23, mean = .16, sd = 3.99)
  res1 = equiv.paired.t.test(diff, alpha = .05, eps_std=.5, power.out = TRUE)

  # check against Wellek function
  eps = c(str_extract(res3,pattern= "theta1 = -?\\d+\\.*\\d*"),str_extract(res3,pattern= "theta2 = -?\\d+\\.*\\d*"))
  eps=c(as.numeric(str_extract(eps[1]," -?\\d+\\.*\\d*")),as.numeric(str_extract(eps[2]," -?\\d+\\.*\\d*")))
  cc = c(str_extract(res3,pattern= "c1 = -?\\d+\\.*\\d*"),str_extract(res3,pattern= "c2 = -?\\d+\\.*\\d*"))
  cc=c(as.numeric(str_extract(cc[1]," -?\\d+\\.*\\d*")),as.numeric(str_extract(cc[2]," -?\\d+\\.*\\d*")))
  pow = c(str_extract(res3,pattern= " POW0 = -?\\d+\\.*\\d*"))
  pow=c(as.numeric(str_extract(pow[1]," -?\\d+\\.*\\d*")))

  # check works for Wellek book function
  expect_equal(res1$eps_std, eps)
  expect_equal(res1$critical.const, cc)
  expect_equal(res1$power, pow)

})


test_that("fisher binomial 1-sided Wellek example works", {

  # One-sided example from Wellek 2010, pg 176
  A = c(rep(1, 98), rep(0, 8))
  B = c(rep(1, 97), rep(0, 10))
  res1=fisher.binom.one.sided(A, B, eps_std=.5, alpha=0.05, power.out = TRUE)
  res3=invisible(capture.output(bi2ste1(length(A), length(B), eps=.5, alpha=0.05, sum(A)/length(A), sum(B)/length(B))))

  # check against Wellek function
  pow = c(str_extract(res3[4],pattern= "POW =  -?\\d+\\.*\\d*"))
  pow=c(as.numeric(str_extract(pow[1]," -?\\d+\\.*\\d*")))

  # check works for Wellek book function
  expect_equal(res1$p, .0499328)
  expect_equal(res1$power, pow)

})

test_that("fisher binomial 2-sided Wellek example works", {

  # Two-sided example from Wellek 2010, pg 191
  A = c(rep(1, 108), rep(0, 117))
  B = c(rep(1, 63), rep(0, 56))
  res1=fisher.binom.two.sided(A, B, eps_std=.6667, alpha=0.05,power.out = TRUE)
  res3a= invisible(capture.output(bi2st(length(A),length(B),rho1=.6667, rho2=1.5000,alpha=.05,s=sum(A,B))))
  res3b= invisible(capture.output(bi2aeq1(length(A),length(B),rho1=.6667, rho2=1.5000,alpha=.05,sum(A)/length(A), sum(B)/length(B))))

  # check against Wellek function
  pow = c(str_extract(res3b[3],pattern= "POW= -?\\d+\\.*\\d*"))
  pow = c(as.numeric(str_extract(pow[1]," -?\\d+\\.*\\d*")))
  cc = c(str_extract(res3a[2],pattern= "C1 = -?\\d+\\.*\\d*"),str_extract(res3a[2],pattern= "C2 = -?\\d+\\.*\\d*"))
  cc=c(as.numeric(str_extract(cc[1]," -?\\d+\\.*\\d*")),as.numeric(str_extract(cc[2]," -?\\d+\\.*\\d*")))

  # check works for Wellek book function
  # expect_equal(res1$critical.const, cc)
  expect_equal(res1$power, pow, tolerance=1e-4)

})


