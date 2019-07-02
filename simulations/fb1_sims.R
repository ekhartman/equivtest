###### Simulations for Fisher-Binomial 1-sided

# A = c(rep(1, 98), rep(0, 8))
# B = c(rep(1, 97), rep(0, 10))
# fisher.binom.one.sided(A, B, eps_sub=.5, alpha=0.05)
# bi2ste1(length(A), length(B), eps=.5, alpha=0.05, sum(A)/length(A), sum(B)/length(B))

library(equivtest)
library(reshape2)
library(parallel)
library(ggplot2)

rm(list = setdiff(ls(), lsf.str()))

set.seed(555)

n_simul = 1000
eps = .4

## you can load the results of the simulation
## which take a while to run, or you can
## run for yourself

########## Code for Windows PCs
########## See code in original "Power_Simulations" file to switch back to Mac-friendly version

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("n_simul","eps", "fisher.binom.one.sided","melt", "odds.ratio"))


# if(file.exists("./fb_1_simulations.RData")) {
#  load("./fb_1_simulations.RData")
# } else {
  res = as.data.frame(do.call("rbind", parLapply(cl=cl, c(10, 50, 125, 140, 175, 200, 500, 1000, 10000)
     , function(n) {
         res = as.data.frame(do.call("rbind", lapply(c(seq(.1, .6, by=.05), .605, .615, .625, .635, seq(.65, 1.5, by= .05)), function(ORatio) {
           res = as.data.frame(do.call("rbind", lapply(1:n_simul, function(i){

             A_prop = .5*ORatio/(.5+.5*ORatio)

             A <- sample(c(0,1), size=n, prob=c(1-A_prop, A_prop), replace=T)
             B <- sample(c(0,1), size=n, prob=c(.5, .5), replace=T)

             equiv = equivtest::fisher.binom.one.sided(A, B, alpha=0.1, eps_std = eps)
             e = equiv$rej

             # equiv2 = fisher.binom.one.sided2(A, B, eps=eps, alpha=.1)
             # e2 = equiv2$rej

             return(list(e = e)) #, e2 = e2))

          })
        ))

        res = lapply(res, function(x) unlist(x))

        return(list(ORatio = ORatio, e_pow = mean(res$e))) #, e_pow2 = mean(res$e2)))
      })
      ))

      res = as.data.frame(lapply(res, function(x) unlist(x)))
      res$n = n
      return(res)
     })
  ))

  res = as.data.frame(lapply(res, function(x) unlist(x)))

  res2 <- melt(res, id=c("ORatio","n"))
  res2$n = paste0("Nt = Nc = ", res2$n)
  res2$n = factor(res2$n, levels = paste0("Nt = Nc = ", unique(res$n)))

stopCluster(cl)

beep(5)

resplot = res2[which(res2$variable=="e_pow"),]

ggplot(data=resplot, aes(x=ORatio, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("Odds Ratio")+
  geom_line(size=.5)+geom_vline(aes(xintercept=.6))+
  ggtitle("FB1: Estimated power over range of odds ratio (new)")

