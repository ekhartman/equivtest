# equivtest

This project was generously supported by a Methods Standards seed grant from [EGAP](https://www.egap.org).


Equivalence testing procedures

Once you have the package files in R:

```{r]

### delete old versions of package out of Packages window in RStudio
### refresh the project for the package

library(devtools)
library(roxygen2)

# if you want to run individual functions instead of loading the package, 
# you may need to load these required packages
# library(BiasedUrn)
# library(testthat)
# library(Exact)
# library(Hmisc)
# library(plm)
# library(sandwich)
# library(EQUIVNONINF)
# library(stringr)

# set to package file directory
# setwd("C:/Users/Sydney/OneDrive/Research/Equivalence package/equivtest")

document()
setwd("..")
library(roxygen2)
install("equivtest")
library(equivtest)

### running unit tests

library(testthat)
# setwd("C:/Users/Sydney/OneDrive/Research/Equivalence package/equivtest")
devtools::test()
```

If installing after downloading the zip-file from Github:

- Download zip-file 
- Copy the 'equivtest-master' folder to computer
- Rename as 'equivtest'
- Use the code above, setting working directory appropriately
