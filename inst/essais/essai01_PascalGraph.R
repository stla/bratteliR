library(bratteli)

setwd("C:/SL/MyPackages/bratteli/inst/essais")

# the Pascal graph ####
Pascal <- function(n) {
  M <- matrix(0, nrow = n+1, ncol = n+2)
  for(i in 1:(n+1)) {
    M[i, ][c(i, i+1L)] <- 1
  }
  M
}

bratteliGraph("Pascal.tex", Pascal, 3)
