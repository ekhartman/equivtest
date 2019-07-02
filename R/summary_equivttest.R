##' Equivalence t-Test Summary Function
##'
##' \code{summary.equivttest} is a summary function for objects of the class \code{equivttest}. This summary function should be paired with the \code{equiv.t.test} function.
##'
##' @seealso \code{\link{equiv.t.test}}
##'
##' @examples
##' # Wellek p 124
##' x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
##' y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
##' eps=c(.5,1)
##' res=equiv.t.test(x,y,eps_sub=eps)
##' summary(res)
##' @export
summary.equivttest <- function(res){

    cat("Equivalence t-test \n")
    cat(paste0("Input: ", res$intype, ", ", res$std.err.type, " SE = ", round(res$SE,3), "\n"))
    cat(paste("T-statistic critical interval:",  round(res$critical.const[1],3),"to",  round(res$critical.const[2],3)), "\n")
    cat(paste("Substantive equivalence CI:",  round(res$CI_sub[1],3),"to",  round(res$CI_sub[2],3)), "\n")
    cat(paste("Standardized equivalence CI:",  round(res$CI_std[1],3),"to",  round(res$CI_std[2],3)), "\n")
    cat(paste0("Reject the null hypothesis? ", paste0(res$rej), ", p-value of ",round(res$p,3), "\n"))
    if(!is.na(res$power)) cat(paste("Power of the test =", round(res$power,3)), "\n")

}
