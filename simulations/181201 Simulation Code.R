####### INSTALL PACKAGE FOLLOWING INSTRUCTIONS ON GITHUB ##############

# library(devtools)
# library(roxygen2)
#
# library(BiasedUrn)
# library(testthat)
# library(Exact)
# library(Hmisc)
# library(plm)
# library(sandwich)
# library(EQUIVNONINF)
# library(stringr)
#
# setwd("C:/Users/Nga Nguy/Dropbox (MIT)/Documents/Works/G4 PML RA/Equivtest/equivtest")
# setwd("D:/Dropbox (MIT)/Documents/Works/G4 PML RA/Equivtest/equivtest")
#
# document()
# setwd("..")
#
# install("equivtest")
# library(equivtest)

#########################################################################


sim_equiv.t.test <- function(n, nsim, mu_x, mu_y, sigma_x, sigma_y) {

  results <- rep(NA, nsim)
  p <- rep(NA, nsim)
  tost_rej <- rep(NA, nsim)

  for(i in 1:nsim) {
    x <- rnorm(n, mu_x, sigma_x)
    y <- rnorm(n, mu_y, sigma_y)

    equiv <- equiv.t.test(x, y, power.out = T, eps_sub = c(-1, 1))
    results[i] <- equiv$rej
    p[i] <- equiv$p

    tost_rej[i] <- tost(x,y,1)$rej
  }

  return(list(results = results, p = p, tost_rej = tost_rej))

}

range <- c(-70.75, seq(-10, 10), 70.75)
test_range <- sapply(range, function(mean_y) {
  test <- sim_equiv.t.test(10000, 500, 0, mean_y, 1, 1)
  power_results <- mean(test$results)
  power_p_left <- mean(test$p < 0.05) # should we output upper-tailed test instead?
  power_p_right <- mean(test$p > 0.95) # should we output upper-tailed test instead?
  power_tost <- mean(test$tost_rej)

  return(c(power_results, power_p_left, power_p_right, power_tost))
})
colnames(test_range) <- range
rownames(test_range) <- c("results", "p_left", "p_right", "tost")
test_range

## THIS WEIRD COMBINATION IS REJECTED HALF THE TIME -- SHOULD HAVE BEEN VERY DIFFICULT TO REJECT!!
set.seed(02142)
x <- rnorm(10000, 0, 1)
y <- rnorm(10000, -70.75, 1)
n <- 10000

equiv.t.test(x,y, eps_sub = c(-1, 1), power.out = T)
tost(x,y,1)

##  THIS IS NEVER REJECTED -- THE REJECTION CALCULATION CURRENTLY DOES NOT TAKE ABS() of t.stat
set.seed(02142)
x <- rnorm(10000, 0, 1)
y <- rnorm(10000, 70.75, 1)
n <- 10000

equiv.t.test(x,y, eps_sub = c(-1, 1), power.out = T)
tost(x,y,1)

## THOUGH THIS IS NEVER REJECTED DESPITE BEING HARDER TO REJECT
set.seed(02142)
x <- rnorm(10000, 0, 1)
y <- rnorm(10000, -20, 1)
n <- 10000

equiv.t.test(x,y, eps_sub = c(-1, 1), power.out = T)
tost(x,y,1)

# The T statistics is exactly what we would get under normal t-testing so it's probably right
T_stat <- ((mean(x) - mean(y)) * sqrt(n*n*(n+n-2)/(n+n))) / sqrt((n-1)*var(x) + (n-1)*var(y))
t.test(x, y)

# I suspect there may be some issues with calculating the pooled standard errors and the epsilon
pooled = sqrt((n - 1) * var(x) + (n - 1) * var(y))/sqrt(n * n * (n + n - 2)/(n + n))
pooled
epsilon <- 1/pooled
epsilon

sqrt(qf(.05, 1, n+n-2, ncp = n*n*epsilon^2/(n+n)))
1 - pf(abs(T_stat)^2, 1, n+n-2, ncp = n*n*epsilon^2/(n+n))

# probably something wrong with using the EQUIVNONINF?

tt2st(10000, 10000, .05, -70.38616, 70.38616, tol=1e-6,itmax = 30)

plot(density(x), xlim=c(-10, 80))
lines(density(y))

plot(density(x/pooled), xlim = c(-400,6000))
lines(density(y/pooled))
