### Erin's 'Figure SI-3: Power Simulations' Code from 2018 paper

library(equivtest)
library(reshape2)
library(parallel)
rm(list = setdiff(ls(), lsf.str()))

set.seed(555)

n_simul = 100
equiv.range = .2

## you can load the results of the simulation
## which take a while to run, or you can
## run for yourself

########## Code for Windows PCs

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("n_simul","equiv.range","equiv.t.test", "melt", "wtd.mean", "wtd.var", "tt2st", "str_extract"))


# if(file.exists("./simul_equiv_vs_90CI_res.RData")) {
#   load("./simul_equiv_vs_90CI_res.RData")
# } else {
  res = as.data.frame(do.call("rbind", parLapply(cl=cl, c(10, 50, 125, 140, 175, 200, 500, 1000, 10000)
    , function(n) {
      res = as.data.frame(do.call("rbind", lapply(c(seq(-2, -.2, 0.1),
                                                        -.1995,-.195,
                                                        seq(-.15,.15,.05),
                                                        .195,.1995,
                                                        seq(.2, 2, 0.1)), function(eps) {
        res = as.data.frame(do.call("rbind", lapply(1:n_simul, function(i){

          X = rnorm(n, eps, 1)
          Y = rnorm(n, 0, 1)

          equiv = equiv.t.test(X, Y, alpha=0.1, eps_std = equiv.range)
          e_ci = equiv$CI_std[1] < eps & eps < equiv$CI_std[2]
          e = equiv$rej


          # equiv2 = equiv.t.test2(X, Y, alpha=0.1, epsilon = equiv.range)
          # e2 = equiv$rej
          # e2 = equiv2$inverted
          # p2 = equiv2$p

          # e_ci2 = -1 * equiv$inverted < eps & eps < inverted
          # e_pow = equiv$power

          # t = t.test(X, Y, conf.level = 0.9)$conf
          # print(t)
          # ci = (t[1] > -1 * equiv.range) & (t[2] < equiv.range)
          # ci_0 = t[1] < eps & t[2] > eps
          # return(list(ci = ci, e_ci = e_ci, ci_0 = ci_0, e = e, e2=e2)) #, p2=p2)) # ,  # , e_pow = e_pow, e2=e2)
          return(list(e = e, e_ci= e_ci))

        })
        ))

        res = lapply(res, function(x) unlist(x))
        return(list(eps = eps, e_pow = mean(res$e), e_inv_cov = mean(res$e_ci, na.rm = TRUE)))

        # return(list(eps = eps, ci_pow = mean(res$ci), e_inv_cov = mean(res$e_ci, na.rm = TRUE),
        #             ci_cov = mean(res$ci_0), e2=mean(res$e2), e_pow = mean(res$e))) # p=mean(res$p), p2 =mean(res$p2), e = mean(res$e),  #,  # e_pow2 = mean(res$e2)))
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

# eps=c(seq(-2, -.2, 0.1),
#       -.1995,-.195,
#       seq(-.15,.15,.05),
#       .195,.1995,
#       seq(.2, 2, 0.1))

resplot = res2[which(res2$variable=="e_pow"),]

ggplot(data=resplot, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+geom_vline(aes(xintercept=-.2))+geom_vline(aes(xintercept = .2))+
  ggtitle("t-Test: Estimated power over range of mu_x-mu_y (new)")

resplot4 = res2[which(res2$variable=="e_inv_cov"),]

ggplot(data=resplot4, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+
  ggtitle("t-Test: Estimated e_inv_cov over range of mu_x-mu_y (new)")


saveRDS(res2, "C:/Users/skahm/OneDrive/Research/Equivalence package/Simulation Results/symm_t.RDS")



library(equivtest)
library(reshape2)
library(parallel)
rm(list = setdiff(ls(), lsf.str()))

set.seed(555)

n_simul = 1000
equiv.range = c(-.5,1)

## you can load the results of the simulation
## which take a while to run, or you can
## run for yourself

########## Code for Windows PCs

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("n_simul","equiv.range","equiv.t.test", "melt", "wtd.mean", "wtd.var", "tt2st", "str_extract"))


# if(file.exists("./simul_equiv_vs_90CI_res.RData")) {
#   load("./simul_equiv_vs_90CI_res.RData")
# } else {
res = as.data.frame(do.call("rbind", parLapply(cl=cl, c(10, 50, 125, 140, 175, 200, 500, 1000, 10000)
     , function(n) {
        res = as.data.frame(do.call("rbind", lapply(c(seq(-2, -.5, 0.1),
                                                      -.4995, -.495,
                                                      seq(-.4,.9,.1),
                                                      .995, .9995,
                                                      seq(1, 2, 0.1)), function(eps) {
        res = as.data.frame(do.call("rbind", lapply(1:n_simul, function(i){

                            X = rnorm(n, eps, 1)
                            Y = rnorm(n, 0, 1)

                            equiv = equiv.t.test(X, Y, alpha=0.1, eps_std = equiv.range)
                            e_ci = equiv$CI_std[1] < equiv.range[1] & equiv.range[2] < equiv$CI_std[2]
                            e = equiv$rej

                            return(list(e = e, e_ci= e_ci))

              })
              ))
        res = lapply(res, function(x) unlist(x))
        return(list(eps = eps, e_pow = mean(res$e), e_inv_cov = mean(res$e_ci, na.rm = TRUE)))

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

eps=c(seq(-2, -.5, 0.1),
      -.4995, -.495,
      seq(-.4,.9,.1),
      .995, .9995,
      seq(1, 2, 0.1))

resplot = res2[which(res2$variable=="e_pow"),]

ggplot(data=resplot, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+geom_vline(aes(xintercept=-.5))+geom_vline(aes(xintercept = 1))+
  ggtitle("t-Test: Estimated power over range of mu_x-mu_y (new)")

resplot4 = res2[which(res2$variable=="e_inv_cov"),]

ggplot(data=resplot4, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+
  ggtitle("t-Test: Estimated e_inv_cov over range of mu_x-mu_y (new)")

saveRDS(res2, "C:/Users/skahm/OneDrive/Research/Equivalence package/Simulation Results/nonsymm_t.RDS")
