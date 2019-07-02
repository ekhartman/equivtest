equiv.paired.t.test2 <- function(diff, alpha = .05, epsilon = .2) {
  dbar <- mean(diff) #estimate
  d.sd <- sd(diff)
  n <- length(diff)
  non.cent <- n*epsilon^2 #non-centrality parameter for the F distribution
  critical.const <- sqrt(qf(alpha, 1, n-1, non.cent)) #critical constant
  t.stat <-  abs(dbar/(d.sd/sqrt(n)))
  # inverted <- uniroot(function(x) pf(abs(t.stat)^2, 1, n-1, ncp = (n*x^2)) - alpha, c(0,2*abs(t.stat)), tol = 0.0001)$root  #invert test-statistic
  rej = abs(t.stat) <= critical.const
  return(list(t.stat = t.stat, critical.const = critical.const, power = 2*pt(critical.const, n-1)-1, rej = rej, inverted = NA))
}

equiv.t.test2 <- function(x, y, alpha = .05, epsilon = .2, std.err = "nominal", cluster.x = NULL, cluster.y = NULL) {
  # WHAT TO DO WITH NAs?
  missing_x = is.na(x)
  missing_y = is.na(y)

  x = x[!missing_x]
  y = y[!missing_y]

  if(!is.null(cluster.x)) {
    cluster.x = cluster.x[!missing_x]
    cluster.y = cluster.y[!missing_y]
  }

  dbar <- mean(x) - mean(y)
  m <- as.double(length(x))
  n <- as.double(length(y))
  N <- m+n
  x.var <-  var(x)
  y.var <- var(y)

  if(!is.null(cluster.x)) {
    deff.x = 1 + (sum(table(cluster.x))^2/sum(table(cluster.x)) - 1) * max(deff(x, cluster.x)["rho"], 0)
    deff.y = 1 + (sum(table(cluster.y))^2/sum(table(cluster.y)) - 1) * max(deff(y, cluster.y)["rho"], 0)
    m.eff = length(x)/deff.x
    n.eff = length(y)/deff.y
    N.eff = m.eff + n.eff
    non.cent = (m.eff*n.eff*epsilon^2)/N.eff
    critical.const <- sqrt(qf(alpha,1,N.eff,non.cent))
  } else {
    non.cent <- (m*n*epsilon^2)/N
    critical.const <- sqrt(qf(alpha,1,N-2,non.cent))
  }

  # pooled.sd = sqrt(((m-1)*x.var+(n-1)*y.var)/(m+n-2))
  pooled.se = sqrt((m-1)*x.var + (n-1)*y.var) / sqrt(m*n * (N-2)/N)
  # pooled.se = pooled.sd * (1/m + 1/n)

  switch(std.err
         , nominal = {
           # se = sqrt((m-1)*x.var + (n-1)*y.var) / sqrt(m*n * (N-2)/N)
           se = pooled.se
           df = N - 2
         }
         , robust = {
           se = sqrt(vcovHC(lm(out ~ treat, data = data.frame(out = c(x, y), treat = c(rep(1, length(x)), rep(0, length(y))))))[2,2])
           df = N - 2
         }
         , cluster = {
           if(is.null(cluster.x) || is.null(cluster.y)) stop("ERROR: Must enter clustering variable for both x (cluster.x) and y (cluster.y) to do clustered standard errors.")

           t.test.clust = t.test.cluster(y = c(x, y), cluster = c(cluster.x, cluster.y), group = c(rep(1, length(x)), rep(0, length(y))))
           se = t.test.clust["S.E. of Effect",1]
           df = t.test.clust[9,1]
         })
  t.stat <-  dbar / se


  p = pf(abs(t.stat)^2, 1, df , non.cent)
  obs_smd = (mean(x) - mean(y)) / sd(y)




  if(!is.null(cluster.x)) {
    inverted <- try(uniroot(function(x) pf(abs(t.stat)^2, 1, N.eff-2, ncp = (m.eff*n.eff*x^2)/N) - ifelse(pf(abs(t.stat)^2, 1, N.eff-2, ncp = (m.eff*n.eff*0^2)/N.eff) < alpha, pf(abs(t.stat)^2, 1, N.eff-2, ncp = (m.eff*n.eff*obs_smd^2)/N.eff), alpha), c(0,10*abs(t.stat)), tol = 0.0001)$root, silent = TRUE)
  } else {
    inverted <- try(uniroot(function(x) pf(abs(t.stat)^2, 1, N-2, ncp = (m*n*x^2)/N) - ifelse(pf(abs(t.stat)^2, 1, N-2, ncp = (m*n*0^2)/N) < alpha, pf(abs(t.stat)^2, 1, N-2, ncp = (m*n*obs_smd^2)/N), alpha), c(0,10*abs(t.stat)), tol = 0.0001)$root, silent = TRUE)
  }

  if(class(inverted) == "try-error") {
    inverted = NA
  }
  rej = abs(t.stat) <= critical.const
  return(list(t.stat = t.stat, critical.const = critical.const, power = 2*pt(critical.const, ifelse(!is.null(cluster.x),N.eff, N-2))-1, rej = rej, p = p, inverted = inverted))
}

tost2 = function(x, y, frac, alpha = 0.05) {

  if(length(unique(x)) == 2) {
    delta = frac
  } else {
    delta = frac * min(sqrt(var(x)), sqrt(var(y)))
  }

  t.u = t.test(x, y, mu = delta, alternative = "less")$statistic
  t.u.df = t.test(x, y, mu = delta, alternative = "less")$parameter
  t.l = t.test(x, y, mu = -delta, alternative = "greater")$statistic
  t.l.df = t.test(x, y, mu = -delta, alternative = "greater")$parameter
  rej = as.logical(ifelse((t.u < -1 * pt(alpha, t.u.df)) & (t.l > pt(alpha, t.l.df)), TRUE, FALSE))

  return(list(rej = rej))
}

##' odds ratio
##'
##' @export
odds.ratio = function(p1, p2) {
  return((p1 * ( 1 - p2)) /( (1 - p1) * p2 ))
}

##' fisher binom one sided
##'
##' @import BiasedUrn Exact
##' @export
fisher.binom.one.sided2 = function(x, y, eps, alpha) {
  x.sum = sum(x == 1)
  m = length(x)
  p1 = x.sum / m
  y.sum = sum(y == 1)
  n = length(y)
  p2 = y.sum / n
  odds = odds.ratio(p1, p2)
  s = x.sum + y.sum
  #	cat("n:", n, "m:", m, "s:", s, "\n")

  # 	## Exact computation--slow
  # 	calc.p = function(x, m, n, s, eps) {
  # 		num = NULL
  # 		for(i in x:m) {
  # 			num = c(num, choose(m, i) * choose(n, s - i) * (1 - eps)^i)
  # 		}
  #
  # 		denom = NULL
  # 		for(i in max(0, s - n):min(s, m)) {
  # 			denom = c(denom, choose(m, i) * choose(n, s - i) * (1 - eps)^i)
  # 		}
  #
  # 		return(sum(num) / sum(denom))
  # 	}
  #
  # 	if(p1 >= p2) {
  # 		p = calc.p(x.sum, m, n, s, eps)
  # 	} else {
  # 		p = calc.p(y.sum, n, m, s, eps)
  # 	}

  p = fisher.test(matrix(c(x.sum, y.sum, m - x.sum, n - y.sum),
                         nrow = 2,
                         dimnames =
                           list(c("A", "B"),
                                c("Positive", "Negative"))), or = 1 - eps, alternative = "greater")$p

  if(p < alpha) {
    rej = 1
  } else {
    rej = 0
  }

  # x = min(s, m)
  # repeat{
  # 	k.alpha = x >= max(0, s - n) & x <= min(s, m)
  # 	k.alpha.p = calc.p(x + 1, m, n, s, eps) <= alpha
  # 	if(k.alpha + k.alpha.p == 2) {
  # 		x = x - 1
  # 	} else {
  # 		break
  # 	}
  # }
  # k.alpha = x + 1
  #
  # gamma = (alpha - calc.p(k.alpha + 1, m, n, s, eps)) / (calc.p(k.alpha, m, n, s, eps) - calc.p(k.alpha + 1, m, n, s, eps))
  #
  # num = NULL
  # for(j in (k.alpha + 1): min(s, m)) {
  # 	num = c(num, log(choose(m, j) * choose(n, s - j) * odds^j + gamma * choose(m, k.alpha) * choose(n, s - k.alpha) * odds^k.alpha))
  # }
  #
  # denom = NULL
  # for(j in max(0, s - n):min(s, m)) {
  # 	denom = c(denom, log(choose(m, j) * choose(n, s - j) * odds^j))
  # }
  #
  # log.power = logSum(num) - logSum(denom)

  # NOTE: power is off by a little... just not sure why

  return(list(odds.ratio = odds, p = p, rej = rej))
}


# # One-sided example (from Wellek 2010 pg 176)
# A = c(rep(1, 98), rep(0, 8))
# B = c(rep(1, 97), rep(0, 10))
#
# print(mean(A))
# print(mean(B))
#
# fisher.binom.one.sided(A, B, .5, 0.05)
# fisher.binom.one.sided(B, A, .5, 0.05)


##' fisher binom two sample
##'
##' @export
fisher.binom.two.samples2 = function(x, y, equiv.or, alpha) {
  # WHAT TO DO ABOUT NAs?
  x = x[!is.na(x)]
  y = y[!is.na(y)]

  equiv.or = ifelse(equiv.or < 1, equiv.or, 1/equiv.or)

  x.sum = sum(x == 1)
  m = length(x)
  p1 = x.sum / m
  y.sum = sum(y == 1)
  n = length(y)
  p2 = y.sum / n
  odds = odds.ratio(p1, p2)
  s = x.sum + y.sum

  #   p = max(fisher.test(matrix(c(x.sum, y.sum, m - x.sum, n - y.sum),
  #                          nrow = 2,
  #                          dimnames =
  #                            list(c("A", "B"),
  #                                 c("Positive", "Negative"))), or = equiv.or, alternative = "greater")$p
  #           , fisher.test(matrix(c(x.sum, y.sum, m - x.sum, n - y.sum),
  #                                nrow = 2,
  #                                dimnames =
  #                                  list(c("A", "B"),
  #                                       c("Positive", "Negative"))), or = 1 / equiv.or, alternative = "less")$p)

  p = max(
    # prob greater than observed with null lower ratio
    sum(dFNCHypergeo(x.sum:m, m, n, s, equiv.or))
    # prob less than observed with upper ratio
    , sum(dFNCHypergeo(0:(x.sum), m, n, s, 1/equiv.or))
  )

  if(p < alpha) {
    rej = 1
  } else {
    rej = 0
  }

  critical.vals.x = c(qFNCHypergeo(alpha, m, n, s, equiv.or, lower.tail = FALSE)
                      , qFNCHypergeo(alpha, m, n, s, 1 / equiv.or))
  critical.vals.y = c(qFNCHypergeo(alpha, n, m, s, equiv.or, lower.tail = FALSE)
                      , qFNCHypergeo(alpha, n, m, s, 1 / equiv.or))
  critical.vals = matrix(c(critical.vals.x, critical.vals.y), ncol = 2
                         , dimnames = list(c("lower", "upper"),
                                           c("X", "Y")))

  inverted = try(uniroot(function(x) max(sum(dFNCHypergeo(x.sum:m, m, n, s, x))
                                         , sum(dFNCHypergeo(0:(x.sum), m, n, s, 1/x))
  ) - alpha
  , c(0.00001, 10), tol = 0.0001)$root)
  if(class(inverted) == "try-error") {
    inverted = NA
  }


  # calculate odds ratio range
  crit.p1.l = critical.vals.x[1] / m
  crit.p1.u = critical.vals.x[2] / m
  crit.p2.l = critical.vals.y[1] / n
  crit.p2.u = critical.vals.y[2] / n

  # power

  ########################################################################

  #   browser()
  #   tmpdir = getwd()
  #   setwd("~/Dropbox/projects/equivalence/rcode/")
  #   out <- c(m,n,equiv.or,1/equiv.or,alpha,p1,p2)
  #   write(out,file="./aeq_inp",ncolumns=1)
  #   dyn.load("bi2aeq1.so")
  #   is.loaded("bi2aeq1_")
  #   .Fortran("bi2aeq1")
  #   power_temp = as.numeric(strsplit(readLines("aeq1_out")[5], "[= ][ ]")[[1]][c(2,8)])
  #   power = power_temp[2]
  #   power.nr = power_temp[1]
  #   setwd(tmpdir)

  power = NA
  power.nr = NA

  # returns an odds ratio, which isn't really the critical value...
  return(list(odds.ratio = odds, p = p, rej = rej
              , critical.vals = c(max(odds.ratio(crit.p1.l, crit.p2.u), odds.ratio(crit.p2.l, crit.p1.u))
                                  , min(odds.ratio(crit.p1.u, crit.p2.l), odds.ratio(crit.p2.u, crit.p1.l)))
              , inverted = c(inverted, 1 / inverted), power = power, power.nr = power.nr
  ))
}



mann.whit2 = function(alpha, x, y, eps.1, eps.2) {
  m = length(x)
  n = length(y)

  eq.ctr = 0.5 * ( 1 + eps.2 - eps.1)
  eq.length = eps.1 + eps.2

  Wxy = wilcox.test(x = x, y = y)$statistic / (m * n)

  if(m + n < 300) {
    PIxyy = 0
    for (i in 1:m) {
      for (j.1 in 1:(n - 1)) {
        for (j.2 in (j.1 + 1):n) {
          PIxyy = PIxyy + trunc(0.5*(sign(x[i] - max(y[j.1],y[j.2])) + 1))
        }
      }
    }

    PIxxy = 0
    for (i.1 in 1:(m-1)) {
      for (i.2 in (i.1+1):m) {
        for (j in 1:n) {
          PIxxy = PIxxy + trunc(0.5*(sign(min(x[i.1],x[i.2]) - y[j]) + 1))
        }
      }
    }

    PIxxy = PIxxy * 2 / (m * (m-1) * n)
    PIxyy = PIxyy * 2 / (n * (n-1) * m)


    sigma = as.numeric(sqrt((Wxy - (m + n - 1) * Wxy^2 + (m - 1) * PIxxy + (n - 1) * PIxyy) / (m * n)))
  } else {
    sigma = sqrt((n * m * (n + m + 1)) / 12) / (m * n)
  }

  crit = sqrt(qchisq(alpha, 1, ((eps.1 + eps.2)^2/(4 * sigma^2))))

  p = pchisq(((Wxy - eq.ctr)/sigma)^2, 1, ((eps.1 + eps.2)^2/(4 * sigma^2)))

  if (abs((Wxy - eq.ctr) / sigma) >= crit) {
    rej = 0
  }
  if (abs((Wxy - eq.ctr) / sigma) < crit) {
    rej = 1
  }

  if (is.na(sigma) || is.na(crit)) {
    rej = 0
  }

  return(list(Wxy = Wxy, crit = crit, rej = rej, se = sigma, p = p))
}


mann.whit.ties2 = function(alpha, x, y, eps.1, eps.2) {
  m = length(x)
  n = length(y)

  eq.ctr = 0.5 * ( 1 + eps.2 - eps.1)
  eq.length = eps.1 + eps.2

  Wxy = 0
  pi0 = 0
  for(i in 1:m) {
    for(j in 1:n) {
      pi0 = pi0 + 1 - sign(abs(x[i] - y[j]))
      Wxy = Wxy + trunc(0.5*(sign(x[i] - y[j]) + 1))
    }
  }
  Wxy = Wxy / (m * n)
  pi0 = pi0 / (m * n)
  Q = Wxy / (1 - pi0)

  PIxyy0 = 0
  PIxyy1 = 0
  PIxyy2 = 0
  for (i in 1:m) {
    for (j.1 in 1:(n-1)) {
      for (j.2 in (j.1+1):n) {
        PIxyy0 = PIxyy0 + 1 - sign(max(abs(x[i] - y[j.1]), abs(x[i] - y[j.2])))
        PIxyy1 = PIxyy1 + trunc(0.5 * (sign(x[i] - max(y[j.1], y[j.2])) + 1))
        PIxyy2 = PIxyy2 + trunc(0.5 * (sign(x[i] - y[j.1]) + 1)) * (1 - sign(abs(x[i] - y[j.2])))
        PIxyy2 = PIxyy2 + trunc(0.5 * (sign(x[i] - y[j.2]) + 1)) * (1 - sign(abs(x[i]-y[j.1])))
      }
    }
  }

  PIxxy0 = 0
  PIxxy1 = 0
  PIxxy2 = 0
  for (i.1 in 1:(m-1)) {
    for (i.2 in (i.1+1):m) {
      for (j in 1:n) {
        PIxxy0 = PIxxy0 + 1 - sign(max(abs(x[i.1] - y[j]), abs(x[i.2]-y[j])))
        PIxxy1 = PIxxy1 + trunc(0.5 * (sign(min(x[i.1], x[i.2]) - y[j]) + 1))
        PIxxy2 = PIxxy2 + trunc(0.5 * (sign(x[i.1] - y[j]) + 1)) * (1 - sign(abs(x[i.2] - y[j])))
        PIxxy2 = PIxxy2 + trunc(0.5 * (sign(x[i.2] - y[j]) + 1)) * (1 - sign(abs(x[i.1] - y[j])))
      }
    }
  }

  PIxxy0 = 2 / (m * (m - 1) * n) * PIxxy0
  PIxyy0 = 2 / (n * (n - 1) * m) * PIxyy0
  PIxxy1 = 2 / (m * (m - 1) * n) * PIxxy1
  PIxyy1 = 2 / (n * (n - 1) * m) * PIxyy1
  PIxxy2 = 1 / (m * (m - 1) * n) * PIxxy2
  PIxyy2 = 1 / (n * (n - 1) * m) * PIxyy2


  var1 = (pi0 - (m + n - 1) * pi0^2 + (m - 1) * PIxxy0 + (n - 1) * PIxyy0) / (m * n)
  var2 = (Wxy - (m + n - 1) * Wxy^2 + (m - 1) * PIxxy1 + (n - 1) * PIxyy1) / (m * n)
  cov = ((m - 1) * PIxxy2 + (n - 1) * PIxyy2 - (m + n - 1) * Wxy * pi0) / (m * n)

  sigma = sqrt( var2 / (1 - pi0)^2 + Wxy^2 * var1 / (1 - pi0)^4 + 2 * Wxy * cov / (1 - pi0)^3)

  crit = sqrt(qchisq(alpha, 1, (((eps.1 + eps.2)^2)/(4 * sigma^2))))

  if (abs((Q - eq.ctr) / sigma) >= crit) {
    rej = 0
  }
  if (abs((Q - eq.ctr) / sigma) < crit) {
    rej = 1
  }

  if (is.na(sigma) || is.na(crit)) {
    rej = 0
  }

  return(list(Q = Q, crit = crit, rej = rej, se = sigma))
}
