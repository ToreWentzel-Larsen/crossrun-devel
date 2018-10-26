#' Joint Probabilities for Crossings and Longest Run,
#' based on the empirical median
#' @description Joint probability distribution for the number of crossings
#' C and the longest run L in a sequence of n independent continuous observations
#' classified as above or below the empirical median. To enhance precision,
#' results are stored in mpfr arrays and the probabilities are multiplied by
#' \eqn{choose(n,m)/2} where m=n/2, even n assumed. The probabilities are
#' integers in this representation.
#'
#' @param nmax max sequence length.
#' @param prec mpft precision.
#' @param printn logical for including progress output.
#' @return nfi, nfn, list of number of subsets of size m, 0<=m<=n
#' with first element included (nfi) and not included (nfn).
#' @examples
#' em10 <- crossrunem(nmax=10,printn=TRUE)
#' @export
crossrunem <- function(nmax = 100, prec = 120,
                       printn = FALSE) {
  # conditioning of S= first value included in the subset (S=1)
  # or not (S=0).
  # nfi: number of subsets, first value included,
  # nfn: number of subsets, first value not included.
  # separate code by brute force for low n (n=1,2,3,4):
  nfi <- list(Rmpfr::mpfrArray(0, prec, dim = c(1, 1, 2)))
  nfn <- list(Rmpfr::mpfrArray(0, prec, dim = c(1, 1, 2)))
  dimnames(nfi[[1]]) <- list("c=0","l=1",c("m=0","m=1"))
  dimnames(nfn[[1]]) <- list("c=0","l=1",c("m=0","m=1"))
  nfn[[1]][1,1,1] <- 1 # m=0
  nfi[[1]][1,1,2] <- 1 # m=1
  nfi[[2]] <- Rmpfr::mpfrArray(0, prec, dim = c(2, 2, 3))
  nfn[[2]] <- Rmpfr::mpfrArray(0, prec, dim = c(2, 2, 3))
  dimnames(nfi[[2]]) <- list(paste0("c=",0:1),paste0("l=",1:2),
                             paste0("m=",0:2))
  dimnames(nfn[[2]]) <- list(paste0("c=",0:1),paste0("l=",1:2),
                             paste0("m=",0:2))
  nfn[[2]][1,2,1] <- 1 # m=0
  nfi[[2]][2,1,2] <- 1 # m=1
  nfn[[2]][2,1,2] <- 1
  nfi[[2]][1,2,3] <- 1 # m=2
  nfi[[3]] <- Rmpfr::mpfrArray(0, prec, dim = c(3, 3, 4))
  nfn[[3]] <- Rmpfr::mpfrArray(0, prec, dim = c(3, 3, 4))
  dimnames(nfi[[3]]) <- list(paste0("c=",0:2),paste0("l=",1:3),
                             paste0("m=",0:3))
  dimnames(nfn[[3]]) <- list(paste0("c=",0:2),paste0("l=",1:3),
                             paste0("m=",0:3))
  nfn[[3]][1,3,1] <- 1 # m=0
  nfi[[3]][2,2,2] <- 1 # m=1
  nfn[[3]][2,2,2] <- 1
  nfn[[3]][3,1,2] <- 1
  nfi[[3]][2,2,3] <- 1 # m=2
  nfi[[3]][3,1,3] <- 1
  nfn[[3]][2,2,3] <- 1
  nfi[[3]][1,3,4] <- 1 # m=3
  nfi[[4]] <- Rmpfr::mpfrArray(0, prec, dim = c(4, 4, 5))
  nfn[[4]] <- Rmpfr::mpfrArray(0, prec, dim = c(4, 4, 5))
  dimnames(nfi[[4]]) <- list(paste0("c=",0:3),paste0("l=",1:4),
                             paste0("m=",0:4))
  dimnames(nfn[[4]]) <- list(paste0("c=",0:3),paste0("l=",1:4),
                             paste0("m=",0:4))
  nfn[[4]][1,4,1] <- 1 # m=0
  nfi[[4]][2,3,2] <- 1 # m=1
  nfn[[4]][2,3,2] <- 1
  nfn[[4]][3,2,2] <- 2
  nfi[[4]][2:3,2,3] <- 1 # m=2
  nfi[[4]][4,1,3] <- 1
  nfn[[4]][2:3,2,3] <- 1
  nfn[[4]][4,1,3] <- 1
  nfi[[4]][2,3,4] <- 1 # m=3
  nfi[[4]][3,2,4] <- 2
  nfn[[4]][2,3,4] <- 1
  nfi[[4]][1,4,5] <- 1 # m=4
  # iterative procedure for higher n (>=5):
  if (nmax>4) for (nn in 5:nmax) {
    nfi[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    nfn[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    dimnames(nfi[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    dimnames(nfn[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    # separate computation for m=0,n:
    nfn[[nn]][1,nn,1] <- 1 # m=0
    nfi[[nn]][1,nn,nn+1] <- 1 # m=n
    # iterative procedure nfi, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:mm) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + nfn[[1]][1,1,1]
        else
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),,mm-gg+1])[,nn-gg]
      }
      if (gg<nn-gg) {
        if (gg==1) {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1]
        }
        else {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1])[,gg]
        }
        nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfn[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm-gg+1]
      } # end low g
    } # end iterative procedure nfi
    # iterative procedure nfn, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:(nn-mm)) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + nfi[[1]][1,1,2]
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),,mm+1])[,nn-gg]
        }
      } # end high g
      if (gg<nn-gg) {
        if (gg==1) {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1]
        }
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1])[,gg]
        }
        nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfi[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm+1]
      } # end low g
    } # end iterative procedure nfi
    if (printn==TRUE) print(nn)
    if (printn==TRUE) print(Sys.time())
  } # end for nn
  return(list(nfi = nfi, nfn = nfn))
} # end function crossrunem

#'
#' Continuation of an existing sequence of joint probabilities for
#' crossings and lomgest run, based on the empirical median.
#' @description  Continuation of an existing sequence of the number of
#' crossings C and the longest run L in a sequence of n independent
#' continuous observations classified as above or below the empirical
#' median. To enhance precision, results are stored in mpfr arrays and
#' the probabilities are multiplied by \eqn{choose(n,m)/2} where m=n/2,
#' even n assumed. The probabilities are integers in this representation.
#'
#' @param emstart existing sequence
#' @param n1 sequence length for the first new case addedc
#' @param nmax max sequence length.
#' @param prec mpft precision.
#' @param printn logical for including progress output.
#' @return nfi, nfn, list of number of subsets of size m, 0<=m<=n
#' with first element included (nfi) and not included (nfn).
crossrunemcont <- function(emstart, n1=61, nmax = 100,
                           prec = 120, printn = FALSE) {
  nfi <- emstart$nfi
  nfn <- emstart$nfn
  for (nn in n1:nmax) {
    nfi[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    nfn[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    dimnames(nfi[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    dimnames(nfn[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    # separate computation for m=0,n:
    nfn[[nn]][1,nn,1] <- 1 # m=0
    nfi[[nn]][1,nn,nn+1] <- 1 # m=n
    # iterative procedure nfi, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:mm) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + nfn[[1]][1,1,1]
        else
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),,mm-gg+1])[,nn-gg]
      }
      if (gg<nn-gg) {
        if (gg==1) {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1]
        }
        else {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1])[,gg]
        }
        nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfn[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm-gg+1]
      } # end low g
    } # end iterative procedure nfi
    # iterative procedure nfn, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:(nn-mm)) {
      if (gg>=nn-gg) {
        if (nn-gg==1)
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + nfi[[1]][1,1,2]
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),,mm+1])[,nn-gg]
        }
      } # end high g
      if (gg<nn-gg) {
        if (gg==1) {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1]
        }
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] +
            cumsumm(nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1])[,gg]
        }
        nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <-
          nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] +
          nfi[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm+1]
      } # end low g
    } # end iterative procedure nfi
    if (printn==TRUE) print(nn)
    if (printn==TRUE) print(Sys.time())
  } # end for nn
  return(list(nfi = nfi, nfn = nfn))
} # end function crossrunemcont
