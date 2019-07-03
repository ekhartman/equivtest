---
title: "equivtest"
output:
  rmarkdown::html_vignette:
    toc: true
    keep_md: true
vignette: >
  %\VignetteIndexEntry{equivtest}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(equivtest)
```


The **equivtest** package implements the equivalence tests discussed in its companion paper, *Hartman, E., & Hidalgo, F. D. (2018). An Equivalence Approach to Balance and Placebo Tests. American Journal of Political Science, 62(4), 1000-1013*.

For each equivalence test in this package, we have included excerpts of the relevant discussion from the companion paper with a worked through example from *Wellek, S. (2010). Testing statistical hypotheses of equivalence and noninferiority. Chapman and Hall/CRC*. 

In general, the equivalence tests have up to three suitable epsilon inputs. The simplest, **eps_tol**, uses 'strict' or 'liberal' equivalence ranges for each test as defined in Wellek 2010. **eps_sub** allows the user to define the equivalence range based upon substantive knowledge of the data. For example, in the equivalence t-test **eps_sub** allows the user to input acceptable ranges for the raw difference in means. Finally, **eps_std** uses a standardized measure to set the equivalence range. For more information, see the help pages for the relevant equivalence testing functions. 

### Equivalence t-tests

As discussed in Hartman & Hidalgo (2018):

*The equivalence range for the t-test for equivalence is typically defined in standardized differences rather than the raw difference in means between the two groups, but researchers can easily map their substantive ranges to standardized differences by scaling by the standard deviation in the covariate. The standardized difference is a useful metric when testing for equivalence because, given some difference between the means of the two distributions, the two groups are increasingly indistinguishable as the variance of the distributions grows towards infinity, and increasingly disjoint as the variance of the distributions shrinks towards zero (Wellek 2010). We also recommend the t-test for equivalence because it is the uniformly most powerful invariant (UMPI) test for two normally distributed variables (Wellek 2010, pg. 120).*

For the equivalence t-test, we are interested in $H_1: \epsilon_L < \frac{\mu_T-\mu_C}{\sigma} < \epsilon_U$. For more information concerning acceptable $\epsilon$ inputs, refer to the *equiv.t.test* documentation or Hartman & Hidalgo (2018).

We now implement Example 6.1 from Wellek (2010). In summary, we wish to compare two treatments using a nonsymmetric equivalence range. For parameter choices refer to the code below. The results interpretation from Wellek (2010) is as follows:

*The data do not allow to reject the null hypothesis of nonequivalence of (treatment A) to (treatment B).*


```r
# Wellek p 124
x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
res=equiv.t.test(x,y,eps_std=c(.5,1), alpha = .05)
summary(res)
#> Equivalence t-test 
#> Input: eps_std,  SE = 2.793
#> T-statistic critical interval: 0.28 to 0.931 
#> Substantive equivalence CI: NA to NA 
#> Standardized equivalence CI: NA to NA 
#> Reject the null hypothesis? FALSE, p-value of NA
```



### Equivalence paired t-test

For the paired equivalence t-test, we are interested in $H_1: \theta_1 < \delta/\sigma_D < \theta_2$, where our differences, $D_i$, are distributed $D_i \sim N(\delta, \sigma_D^2)$. For more information concerning acceptable $\epsilon$ inputs, refer to the *equiv.paired.t.test* documentation or Wellek (2010).

We now implement Example 5.3 from Wellek (2010). In summary, we wish to evaluate the stability in a vector of observed pre-treatment differences. For parameter choices refer to the code below. The results interpretation from Wellek (2010) is as follows:

*We can conclude that at the 5% level, the experimental data of the present example contain sufficient evidence in support of the hypothesis that (the variable) does not change to a relevant extent over the (pre-treatment time interval).*


```r
# Wellek p 96
x= c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
res=equiv.paired.t.test(diff= x, alpha = .05, eps_std=.5, power.out = TRUE)
summary(res)
#> Equivalence paired t-test 
#> Input type:  eps_std 
#> t-test statistic:  2.051 
#> Equivalence CI:  -1.097 to 1.097 
#> Rejected the null? FALSE 
#> Power of the test:  0.215
```



### Fisher-Binomial Tests of Equivalence

As discussed in Hartman & Hidalgo (2018):

*The Fisher type exact test is well adapted to equivalence between two groups with binary outcomes. This test is based on the odds ratio as opposed to the mean difference between the two groups ... We will call the rate of units with a response value of 1 in the treatment condition $p_T$ and the rate of units with a response value of 1 in the control condition $p_C$. The test statistic is the odds ratio of the two groups, $\rho = p_T(1−p_T)/p_C(1−p_C)$.*

#### Fisher-Binomial One-Sided

In the one-sided case, we are interested in $H_1: \rho \le 1- \epsilon$ with $0 < \epsilon < 1$. For more information concerning acceptable $\epsilon$ inputs, refer to the *fisher.binom.one.sided* documentation.

We now implement Example 6.4 from Wellek (2010). In summary, we wish to compare treatment group A with treatment group B using an equivalence test to evaluate the one-sided equivalency between the two treatments. For treatment vectors and parameter choices refer to the code below. The results interpretation from Wellek (2010) is as follows:

*At the 5% level of significance, the data of the trial allow the conclusion that the odds of a favorable response to (treatment A) do not fall by more than 50% as compared to those of a positive outcome after administering (treatment B).*


```r
# One-sided example from Wellek 2010, Example 6.4 pg 176
A = c(rep(1, 98), rep(0, 8))
B = c(rep(1, 97), rep(0, 10))
res = fisher.binom.one.sided(A, B, eps_std= .5, 
                             alpha= 0.05, power.out = TRUE)
summary(res)
#> Equivalence One-Sided Fisher-Binomial Test 
#> Input:  eps_std 
#> Odds Ratio of Treatment Groups A and B:  1.263 
#> Rejected the null? TRUE, p-value of 0.05
#> Power of the test:  0.58
```



#### Fisher-Binomial Two-Sided

In the two-sided case, we are interested in $H_1: \epsilon_L < \rho < \epsilon_U$ with $\epsilon_L < 1 < \epsilon_U$. For more information concerning acceptable $\epsilon$ inputs, refer to the *fisher.binom.two.sided* documentation.

We now implement Example 6.5 from Wellek (2010). In summary, we wish to compare the previous treatment status of subjects (A) with the response to trial medication (B) using an equivalence test to evaluate two-sided equivalency. For treatment vectors and parameter choices refer to the code below. The results interpretation from Wellek (2010) is as follows:

*The (data) do not allow to reject the null hypothesis that there are relevant differences between patients with and without (previous treatment) with respect to the probability of a positive response to the study medication.*


```r

# Two-sided example from Wellek 2010, pg 191
A = c(rep(1, 108), rep(0, 117))
B = c(rep(1, 63), rep(0, 56))
res = fisher.binom.two.sided(A, B, eps_std=.41, alpha=0.05, power.out = TRUE)
summary(res)
#> Equivalence Two-Sided Fisher-Binomial Test 
#> Input:  eps_std 
#> Odds Ratio of Treatment Groups A and B:  0.821 
#> Equivalence Odds Ratio CI:  0.551 to 1.815 
#> Rejected the null? TRUE, p-value of 0.002
#> Power of the test:  0.914
```



### Two One-Sided Test (TOST) 

As discussed in Hartman & Hidalgo (2018):

*Rather than studying the standardized difference, as is used in the equivalence t-test ... researchers may wish to conduct a test for equivalence of the raw mean difference. This can be accomplished using a Two-One-Sided-Test (TOST) (Berger and Hsu 1996). The TOST test is conducted using two one sided t-tests centered around the bounds of the equivalence range. One advantage of the TOST is that it allows for the researcher to define the equivalence range on the scale of the variable of interest as opposed to standardizing substantive ranges.*

In the traditional case, we are interested in $H_1: \epsilon_L < \mu_T-\mu_C < \epsilon_U$.  For more information, refer to the *tost* documentation or Hartman & Hidalgo (2018).

Here we use the equivalence t-test dataset from before, Wellek (2010) Example 6.1, refer to the code below for parameter choices. The observed mean difference is -3.033. Using a population mean difference of 10, we reject the null hypothesis. 


```r
# Wellek p 124
x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
res1=tost(x,y,eps_sub=10, alpha=.1)
summary(res1)
#> Two-One-Sided Test 
#> Lower and Upper Test Statistics: -4.667, 2.495
#> Rejected the null? TRUE
```

### Two One-Sided Ratio Test (TOST Ratio) 

As discussed in Hartman & Hidalgo (2018):

*The TOST test can also be adapted to test for equivalence of the ratio of the means of the two groups, instead of the raw difference between the means. The TOST ratio test has the advantage of having an absolute scale that is independent of the scale of the variable of interest.*

For the TOST ratio test, we are interested in $H_1: \epsilon_L < \frac{\mu_T}{\mu_C} < \epsilon_U$. For more information, refer to the *tost.ratio* documentation or Hartman & Hidalgo (2018).

Here we use the equivalence t-test dataset from before, Wellek (2010) Example 6.1, refer to the code below for parameter choices. The observed ratio is 0.579. Using a ratio of population mean X to population mean Y of 2, we reject the null hypothesis. 


```r

x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
res1=tost.ratio(x,y,frac=2)
#> Warning in tost.ratio(x, y, frac = 2): No eps_sub or eps_std specified,
#> using a default value of eps_std = 1.
summary(res1)
#> Two-One-Sided Ratio Test 
#> Lower and Upper Test Statistics: -2.119, 7.718
#> Rejected the null? TRUE
```



