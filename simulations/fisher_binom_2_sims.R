### Erin's 'Figure SI-3: Power Simulations' Code from 2018 paper

library(equivtest)
library(reshape2)
library(parallel)
rm(list = setdiff(ls(), lsf.str()))

set.seed(555)

n_simul = 1000
equiv.range = .4

## you can load the results of the simulation
## which take a while to run, or you can
## run for yourself

########## Code for Windows PCs

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("n_simul","equiv.range","odds.ratio", "fisher.binom.two.sided","melt"))

  res = as.data.frame(do.call("rbind", parLapply(cl=cl, c(10, 50, 125, 140, 175, 200, 500, 1000, 10000)
    , function(n) {
      res = as.data.frame(do.call("rbind", lapply(c(c(seq(.1, .4, by=.1), .415, .435, seq(.5, .9,by=.1)),
                                                    sort(1/c(seq(.3, .4, by=.1), .415, .435, seq(.5, .9,by=.1)))), function(ORatio) {
        res = as.data.frame(do.call("rbind", lapply(1:n_simul, function(i){

          A_prop = .5*ORatio/(.5+.5*ORatio)

          A <- sample(c(0,1), size=n, prob=c(1-A_prop, A_prop), replace=T)
          B <- sample(c(0,1), size=n, prob=c(.5, .5), replace=T)

          equiv = equivtest::fisher.binom.two.sided(x=A, y=B, alpha=0.1, eps_std = .4)
          e = equiv$rej
          e_ci = equiv$equiv.CI[1] < ORatio & ORatio < equiv$equiv.CI[2]

          return(list(e = e, e_ci= e_ci))

        })
        ))

        res = lapply(res, function(x) unlist(x))
        return(list(ORatio = ORatio, e_pow = mean(res$e), e_inv_cov = mean(res$e_ci, na.rm = TRUE)))

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

library(ggplot2)

resplot = res2[which(res2$variable=="e_pow"),]

ggplot(data=resplot, aes(x=ORatio, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("odds ratio")+
  geom_line(size=.5)+geom_vline(aes(xintercept=.4))+geom_vline(aes(xintercept = 2.5))+
  ggtitle("FB2: Estimated power over range of odds ratio (new)")

resplot4 = res2[which(res2$variable=="e_inv_cov"),]

ggplot(data=resplot4, aes(x=ORatio, y=value, group = as.factor(n), colour = as.factor(n)))+xlab("odds ratio")+
  geom_line(size=.5)+geom_vline(aes(xintercept=.4))+geom_vline(aes(xintercept = 2.5))+
  ggtitle("FB2: Estimated e_inv_cov over range of odds ratio (new)")

