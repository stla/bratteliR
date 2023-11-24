#' @title Bratteli dimensions
#' @description Dimensions of the vertices of a Bratteli graph.
#'
#' @param Mn a function returning for each integer \code{n} the incidence
#' matrix between levels \code{n} and \code{n+1}; the matrix \code{Mn(0)}
#' must have one and only one row
#' @param N the level up to which the dimensions are wanted
#'
#' @return The dimensions of the vertices in a list.
#' @export
#' @importFrom gmp as.bigz
#'
#' @examples
#' # the Pascal graph ####
#' Pascal <- function(n) {
#'  M <- matrix(0, nrow = n+1, ncol = n+2)
#'  for(i in 1:(n+1)) {
#'    M[i, ][c(i, i+1L)] <- 1
#'  }
#'  M
#' }
#' bratteliDimensions(Pascal, 4)
#'
#' # the Euler graph ####
#' Euler <- function(n) {
#'   M <- matrix(0, nrow = n+1, ncol = n+2)
#'   for(i in 1:(n+1)) {
#'     M[i, ][c(i, i+1L)] <- c(i, n+2-i)
#'   }
#'   M
#' }
#' bratteliDimensions(Euler, 4)
bratteliDimensions <- function(Mn, N) {
  Dims <- vector("list", N)
  M <- Mn(0)
  if(nrow(M) != 1L) stop("Mn(0) must have one and only one row.")
  dims0 <- as.bigz(M)
  Dims[[1L]] <- as.character(M)
  for(k in seq_len(N-1)){
    dims0 <- gmp::`%*%`(dims0, as.bigz(Mn(k)))
    Dims[[k+1L]] <- as.character(c(dims0))
  }
  Dims
}

#' @title Bratteli kernels
#' @description Central kernels of a Bratteli graph.
#'
#' @param Mn a function returning for each integer \code{n} the incidence
#' matrix between levels \code{n} and \code{n+1}; the matrix \code{Mn(0)}
#' must have one and only one row
#' @param N the level up to which the kernels are wanted
#'
#' @return The kernels in a list.
#' @export
#' @importFrom gmp as.bigz
#'
#' @examples
#' # the Pascal graph ####
#' Pascal <- function(n) {
#'  M <- matrix(0, nrow = n+1, ncol = n+2)
#'  for(i in 1:(n+1)) {
#'    M[i, ][c(i, i+1L)] <- 1
#'  }
#'  M
#' }
#' bratteliKernels(Pascal, 4)
#'
#' # the Euler graph ####
#' Euler <- function(n) {
#'   M <- matrix(0, nrow = n+1, ncol = n+2)
#'   for(i in 1:(n+1)) {
#'     M[i, ][c(i, i+1L)] <- c(i, n+2-i)
#'   }
#'   M
#' }
#' bratteliKernels(Euler, 4)
bratteliKernels <- function(Mn, N) {
  Kernels <- vector("list", N)
  M <- Mn(0)
  m <- nrow(M)
  n <- ncol(M)
  if(m != 1L) stop("Mn(0) must have one and only one row.")
  dims0 <-  c(as.bigz(M))
  Kernels[[1L]] <- matrix(as.character(dims0), dimnames = list(1L:n, 1L:m))
  for(k in 1L:(N-1)) {
    M <- as.bigz(Mn(k))
    m <- nrow(M)
    n <- ncol(M)
    S <- lapply(1L:ncol(M), function(i) which(M[, i] != 0L))
    dims <- c(gmp::`%*%`(dims0, M))
    P <- lapply(1L:n, function(i) {
      as.character(dims0[S[[i]]] * M[S[[i]], i] / dims[i])
    })
    Kernels[[k+1L]] <-
      matrix("0", nrow = n, ncol = m, dimnames = dimnames(t(M)))
    for(i in 1L:n){
      Kernels[[k+1L]][i, ][S[[i]]] <- P[[i]]
    }
    dims0 <- dims
  }
  Kernels
}

#' @title Intrinsic distances
#' @description Intrinsic distances on a Bratteli graph
#'
#' @param Mn a function returning for each integer \code{n} the incidence
#' matrix between levels \code{n} and \code{n+1}; the matrix \code{Mn(0)}
#' must have one and only one row
#' @param N the level up to which the distances are wanted
#'
#' @return The distance matrices in a list.
#' @export
#' @importFrom gmp as.bigq
#' @importFrom kantorovich kantorovich
#'
#' @examples
#' # the Pascal graph ####
#' Pascal <- function(n) {
#'  M <- matrix(0, nrow = n+1, ncol = n+2)
#'  for(i in 1:(n+1)) {
#'    M[i, ][c(i, i+1L)] <- 1
#'  }
#'  M
#' }
#' bratteliDistances(Pascal, 4)
#'
#' # the Euler graph ####
#' Euler <- function(n) {
#'   M <- matrix(0, nrow = n+1, ncol = n+2)
#'   for(i in 1:(n+1)) {
#'     M[i, ][c(i, i+1L)] <- c(i, n+2-i)
#'   }
#'   M
#' }
#' bratteliDistances(Euler, 4)
bratteliDistances <- function(Mn, N){
  ckernels <- bratteliKernels(Mn, N)
  RHO <- lapply(ckernels, function(kernel) {
    matrix("", nrow = nrow(kernel), ncol = nrow(kernel))
  })
  RHO[[1L]] <- (diag(nrow(ckernels[[1L]])) + 1) %% 2
  storage.mode(RHO[[1L]]) <- "character"
  n <- length(ckernels) - 1L
  for(k in 1L:n){
    diag(RHO[[k+1L]]) <- "0"
    K <- nrow(RHO[[k+1L]])
    kernel <- ckernels[[k+1L]]
    D <- unname(RHO[[k]])
    for(i in 1L:(K-1L)){
      for(j in (i+1L):K){
        RHO[[k+1L]][i, j] <- RHO[[k+1L]][j, i] <- as.character(
          kantorovich(as.bigq(kernel[i, ]), as.bigq(kernel[j, ]), dist = D)
        )
      }
    }
    M <- Mn(k)
    dimnames(RHO[[k+1L]]) <- list(colnames(M), colnames(M))
  }
  RHO
}
