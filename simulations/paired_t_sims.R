### Erin's 'Figure SI-3: Power Simulations' Code from 2018 paper

library(equivtest)
library(reshape2)
library(parallel)
rm(list = setdiff(ls(), lsf.str()))

set.seed(555)

n_simul = 1000
equiv.range = c(-.2, .3)
# equiv.range = c(0, 0.2)

## you can load the results of the simulation
## which take a while to run, or you can
## run for yourself

########## Code for Windows PCs

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("n_simul","equiv.range","equiv.paired.t.test","melt"))


# if(file.exists("./simul_equiv_vs_90CI_res.RData")) {
#   load("./simul_equiv_vs_90CI_res.RData")
# } else {
res = as.data.frame(do.call("rbind", parLapply(cl=cl, c(10, 50, 125, 140, 175, 200, 500, 1000, 10000)
    , function(n) {
        res = as.data.frame(do.call("rbind", lapply(c(seq(-.5, -.3, 0.1),
                                                      -.295, -.275,
                                                      seq(-.2, .2, 0.1),
                                                      .275, .295,
                                                      seq(.3, .5, 0.1)), function(eps) {
           res = as.data.frame(do.call("rbind", lapply(1:n_simul, function(i){

             x = rnorm(n, mean=eps, sd=1)
             y = rnorm(n)

             diff = x-y

             equiv = try(equiv.paired.t.test(diff=diff, alpha=0.1, eps_std = equiv.range, power.out=FALSE),silent=TRUE)

             if(class(equiv) == "try-error") {
               equiv = NA
               e = NA
               e_ci = NA
             } else {
               e = equiv$rej
               e_ci = equiv$CI_std[1] < eps & eps < equiv$CI_std[2]
             }


             return(list(e = e, e_ci = e_ci))
           })
           ))

           res = lapply(res, function(x) unlist(x))

           return(list(eps = eps, e_pow = mean(res$e, na.rm = TRUE), e_inv_cov = mean(res$e_ci, na.rm = TRUE)))
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
beep(5)



library(ggplot2)

# eps=seq(-.5, .5, 0.05)

resplot = res2[which(res2$variable=="e_pow"),]

ggplot(data=resplot, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+geom_vline(aes(xintercept=-.2))+geom_vline(aes(xintercept = .3))+
  # geom_vline(aes(xintercept = .37), size = 1, colour = "red")+
  ggtitle("Paired t-Test: Estimated power over range of diff_mu (new)")


resplot4 = res2[which(res2$variable=="e_inv_cov"),]

ggplot(data=resplot4, aes(x=eps, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("mu_x")+
  geom_line(size=.5)+
  ggtitle("Paired t-Test: Estimated e_inv_cov over range of mu_x-mu_y (new)")


saveRDS(res2, "C:/Users/skahm/OneDrive/Research/Equivalence package/Simulation Results/paired_symm_t.RDS")

