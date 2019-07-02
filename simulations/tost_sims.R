### Erin's 'Figure SI-3: Power Simulations' Code from 2018 paper

library(equivtest)
library(reshape2)
library(parallel)
rm(list = setdiff(ls(), lsf.str()))

set.seed(555)

n_simul = 1000
equiv.range = 0.50

## you can load the results of the simulation
## which take a while to run, or you can
## run for yourself

########## Code for Windows PCs

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("n_simul","equiv.range","tost"))


# if(file.exists("./simul_equiv_vs_90CI_res.RData")) {
#   load("./simul_equiv_vs_90CI_res.RData")
# } else {
  res = as.data.frame(do.call("rbind", parLapply(cl=cl, c(10, 50, 125, 140, 175, 200, 500, 1000, 10000)
    , function(n) {
      res = as.data.frame(do.call("rbind", lapply(c(seq(-1, -.5, 0.1), -.495, -.475, -.45, seq(-.4, .4, .1), .45, .475, .495, seq(.5, 1, 0.1)), function(eps) {
        res = as.data.frame(do.call("rbind", lapply(1:n_simul, function(i){

          X = rnorm(n, eps, 1)
          Y = rnorm(n, 0, 1)

          equiv = tost(X, Y, alpha=0.1, eps_sub = equiv.range)
          e = equiv$rej
          # e_ci = -1 * equiv$CI_std[1] < eps & eps < equiv$CI_std[2]
          # e_pow = equiv$power

          # t = t.test(X, Y, conf.level = 0.9)$conf
          #  print(t)
          # ci = (t[1] > -1 * equiv.range) & (t[2] < equiv.range)
          # ci_0 = t[1] < eps & t[2] > eps
          return(list(e = e)) # e_pow = e_pow, ci = ci, e_ci = e_ci, ci_0 = ci_0))
        })
        ))

        res = lapply(res, function(x) unlist(x))

        return(list(eps = eps, e_pow = mean(res$e))) #,  ci_pow = mean(res$ci), e_inv_cov = mean(res$e_ci, na.rm = TRUE), ci_cov = mean(res$ci_0)))
      })
      ))

      res = as.data.frame(lapply(res, function(x) unlist(x)))
      res$n = n
      return(res)
    })
  ))

  res = as.data.frame(lapply(res, function(x) unlist(x)))

  res2 <- melt(res, id=c("eps","n"))
  res2$n = paste0("Nt = Nc = ", res2$n)
  res2$n = factor(res2$n, levels = paste0("Nt = Nc = ", unique(res$n)))

  # save(res2, file = "./simul_equiv_vs_90CI_res.RData")
# }

stopCluster(cl)

library(ggplot2)

resplot = res2[which(res2$variable=="e_pow"),]

ggplot(data=resplot, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+geom_vline(xintercept = .5) + geom_vline(xintercept=-.5)+
  # geom_vline(aes(xintercept = .37), size = 1, colour = "red")
  ggtitle("TOST: Estimated power over range of mu_x")

