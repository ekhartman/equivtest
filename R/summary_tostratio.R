##' Equivalence Two-One-Sided Ratio Test Summary Function
##'
##' \code{summary.equivtostratio} is a summary function for objects of the class \code{equivtostratio}. This summary function should be paired with the \code{tost.ratio} function.
##'
##' @seealso \code{\link{equiv.tost}}
##'
##' @examples
##' # Wellek p 124
##' x=c(10.3,11.3,2,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
##' y=c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3,6,3.5,18.7,9.6)
##' res=tost.ratio(x,y,frac=2)
##' summary(res)
##' @export
summary.equivtostratio <- function(res){

    cat("Two-One-Sided Ratio Test \n")
    cat(paste0("Lower and Upper Test Statistics: ", round(res$tstat[1],3),", ", round(res$tstat[2],3),"\n"))
    cat(paste0("Rejected the null? ", paste0(ifelse(res$rej==1,"TRUE","FALSE")), "\n"))

}





