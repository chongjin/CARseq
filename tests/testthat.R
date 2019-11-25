library(testthat)
library(CARseq)
library(nloptr)
library(MASS)
library(zetadiv)

# nonnegative functions used a fitter in glm.nb
# do NOT move to unit tests, or the test cannot find the method
glm.fit.cons.nonneg=function(...) {zetadiv::glm.fit.cons(cons=1, ...)}

test_check("CARseq")
