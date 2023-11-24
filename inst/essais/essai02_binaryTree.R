library(bratteli)

setwd("C:/SL/MyPackages/bratteli/inst/essais")

# the binary tree ####
tree <- function(n) {
  M <- matrix(0, nrow = 2^n, ncol = 2^(n+1))
  for(i in 1:nrow(M)) {
    M[i, ][c( 2*(i-1)+1, 2*(i-1)+2 )] <- 2
  }
  M
}

bratteliGraph("binaryTree2.tex", tree, 3)
