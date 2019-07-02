# check t-tests work with critical constants from Table 6.1

library(equivtest)
library(EQUIVNONINF)
library(stringr)

### new code doesn't match

eps = c(.25,.5,.75,1)
n = seq(10,75,by=5)
alpha=.05
res = matrix(NA, ncol=4, nrow=length(n))

for(i in 1:length(eps)){
  for(j in 1:length(n)){

    x = rnorm(n=n[j])
    y = rnorm(n=n[j])

    res[j,i] = (equiv.t.test(x,y,alpha, eps_std=eps[i])$critical.const[2]-
                  equiv.t.test(x,y,alpha, eps_std=eps[i])$critical.const[1])/2

  }
}


### original code matches

eps = c(.25,.5,.75,1)
n = seq(10,75,by=5)
alpha=.05
res2 = matrix(NA, ncol=4, nrow=length(n))

for(i in 1:length(eps)){
  for(j in 1:length(n)){

    x = rnorm(n=n[j])
    y = rnorm(n=n[j])

    res2[j,i] = equiv.t.test2(x,y, alpha, epsilon=eps[i])$critical.const

  }
}


### plain tt2st doesn't match 6.1, matches new equiv.t.test

eps = c(.25,.5,.75,1)
n = seq(10,75,by=5)
alpha=.05
res3 = matrix(NA, ncol=4, nrow=length(n))

for(i in 1:length(eps)){
  for(j in 1:length(n)){


    critconst.pow = invisible(capture.output(tt2st(n[j],n[j],alpha,eps[i],eps[i],tol=1e-6,itmax = 5)))
    critconst2=c(str_extract(critconst.pow,pattern="c1 = -?\\d+\\.*\\d*"),str_extract(critconst.pow,pattern="c2 = -?\\d+\\.*\\d*"))
    critconst3=as.numeric(c(str_extract(critconst2[1]," -?\\d+\\.*\\d*"),str_extract(critconst2[2]," -?\\d+\\.*\\d*")))


    res3[j,i] = (critconst3[1]) #-critconst3[1])/2

  }
}





