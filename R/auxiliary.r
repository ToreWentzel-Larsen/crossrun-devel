#' Row-wise Cumulative Sums
#' @description Row-wise Cumulative Sums in mpfr Array.
#'
#' @param mtrx mpfr two-dimensional array.
#' @return mpfr array with row-wise cumulative sums,
#' same dimension as the original array.
#' @examples
#' nill <- Rmpfr::mpfr(0, 120)
#' one <- Rmpfr::mpfr(1, 120)
#' two <- Rmpfr::mpfr(2, 120)
#' contents <- c(one,nill,nill, one,one,one, two,two,two)
#' mtrx3 <- Rmpfr::mpfr2array(contents, dim = c(3, 3))
#' print(mtrx3)
#' print(cumsumm(mtrx3))
#' @export
cumsumm <- function(mtrx) {
    mtrs <- mtrx
    nrw <- nrow(mtrx)
    for (rw in 1:nrw) mtrs[rw, ] <- cumsum(mtrx[rw, ])
    mtrs
}

#' Column-Wise Cumulative Sums
#' @description Column-wise cumulative sums in
#' mpfr array.
#'
#' @param mtrx mpfr two-dimensional array.
#' @return mpfr array with column-wise cumulative sums,
#' same dimension as the original array.
#' @examples
#' nill <- Rmpfr::mpfr(0, 120)
#' one <- Rmpfr::mpfr(1, 120)
#' two <- Rmpfr::mpfr(2, 120)
#' contents <- c(one,nill,nill, one,one,one, two,two,two)
#' mtrx3 <- Rmpfr::mpfr2array(contents, dim = c(3, 3))
#' print(mtrx3)
#' print(cumsummcol(mtrx3))
#' @export
cumsummcol <- function(mtrx) {
    mtrs <- mtrx
    ncl <- ncol(mtrx)
    for (cl in 1:ncl) mtrs[, cl] <- cumsum(mtrx[, cl])
    mtrs
}

#' Box Cumulative Sums
#' @description A box cumulative sum is defined as the
#' cumulative sum over a lower left rectangle. This function
#' is primarily for use when the components are point
#' probabilities for the number of crossings C and the longest
#' run L, then component (c,l) in the result is the
#' probability \eqn{P(C \ge c, L \le l)}.
#'
#' @param mtrx mpfr array
#' @return mpfr array
#' @examples
#' nill <- Rmpfr::mpfr(0, 120)
#' one <- Rmpfr::mpfr(1, 120)
#' two <- Rmpfr::mpfr(2, 120)
#' contents <- c(one,nill,nill, one,one,one, two,two,two)
#' mtrx3 <- Rmpfr::mpfr2array(contents, dim = c(3, 3))
#' print(mtrx3)
#' print(boxprobt(mtrx3))
#' @export
boxprobt <- function(mtrx) {
  mtrs <- mtrx
  rw <- nrow(mtrx)
  cl <- ncol(mtrx)
  for (row in 1:rw) mtrs[row, ] <- cumsum(mtrx[row, ])
  for (col in 1:cl) mtrs[, col] <- rev(cumsum(rev(mtrs[, col])))
  return(mtrs)
}

#' Exact Joint Probabilities for Low n
#' @description Exact joint probabilities, for low n,
#' of the number of crossings C and the longest run L
#' in n independent Bernoulli observations with success
#' probability p. Probabilites are multiplied by \eqn{2^{n-1}}.
#' @param  n number, length of seqience, at most 6.
#' @param  p success probability.
#' @param prec precision in mpfr calculations.
#' Default 120.
#' @return mpfr array
#' @examples
#' exactbin(n=6)
#' exactbin(n=5, p=0.6)
#' @export
exactbin <- function(n, p = 0.5, prec = 120) {
    nill <- Rmpfr::mpfr(0, prec)
    one <- Rmpfr::mpfr(1, prec)
    two <- Rmpfr::mpfr(2, prec)
    pm <- Rmpfr::mpfr(p, prec)
    qm <- one - pm
    res <- Rmpfr::mpfr2array(one, dim = c(1, 1))
    if (n == 2)
        res <- Rmpfr::mpfr2array(c(nill, two * pm * qm, pm^2 + qm^2, nill),
            dim = c(2, 2))
    if (n == 3)
        res <- Rmpfr::mpfr2array(c(nill, nill, pm * qm, nill, two * pm * qm,
            nill, pm^3 + qm^3, nill, nill), dim = c(3, 3))
    if (n == 4)
        res <- Rmpfr::mpfr2array(c(rep(nill, 3), 2 * pm^2 * qm^2, nill, 2 *
            pm^2 * qm^2, two * pm * qm * (1 - pm * qm), nill, nill,
            two * pm * qm * (pm^2 + qm^2), nill, nill, pm^4 + qm^4,
            rep(nill, 3)), dim = c(4, 4))
    if (n == 5)
        res <- Rmpfr::mpfr2array(c(rep(nill, 4), pm^2 * qm^2, nill, nill, pm *
            qm * (1 - pm * qm), 4 * pm^2 * qm^2, nill, nill, 2 * pm^2 *
            qm^2, pm * qm * (2 * pm^3 + pm * qm + 2 * qm^3), nill, nill,
            nill, 2 * pm * qm * (pm^3 + qm^3), rep(nill, 3), pm^5 +
                qm^5, rep(nill, 4)), dim = c(5, 5))
    if (n == 6)
        res <- Rmpfr::mpfr2array(c(rep(nill, 5), 2 * pm^3 * qm^3, nill, nill,
            pm^4 * qm^2 + pm^2 * qm^4, 2 * pm^4 * qm^2 + 8 * pm^3 *
                qm^3 + 2 * pm^2 * qm^4, 3 * pm^4 * qm^2 + 4 * pm^3 *
                qm^3 + 3 * pm^2 * qm^4, nill, nill, 2 * pm^3 * qm^3,
            2 * pm^5 * qm + 2 * pm^4 * qm^2 + 4 * pm^3 * qm^3 + 2 *
                pm^2 * qm^4 + 2 * pm * qm^5, 4 * pm^4 * qm^2 + 4 * pm^2 *
                qm^4, nill, nill, nill, 2 * pm^4 * qm^2 + 2 * pm^2 *
                qm^4, 2 * pm^5 * qm + pm^4 * qm^2 + pm^2 * qm^4 + 2 *
                pm * qm^5, rep(nill, 3), nill, 2 * pm^5 * qm + 2 * pm *
                qm^5, rep(nill, 4), pm^6 + qm^6, rep(nill, 5)), dim = c(6,
            6))
    rownames(res) <- c(0:(n - 1))
    colnames(res) <- c(1:n)
    return(res)
}

#' Number of Crossings and Longest Run
#' @description Auxiliary function for simclbin, computing
#' the number of crossings (type=0) or longest run (type=2)
#' in a sequence of independent normal observations. Crossings
#' and runs are related to whether the observations are above
#' a shift.
#' @param seri numeric; seri a sequence of random draws
#' @param shift numeric; shift for the observatoobs
#' @param type numeric; 0 number of crossings, 1 longest run
#' @return number of crossings or longest run, numeric
#' @export
clshift <- function(seri, shift = 0, type = 0) {
    rle.sh <- rle(seri + shift > 0)$lengths
    if (type == 0)
        res <- length(rle.sh) - 1
    if (type == 1)
        res <- max(rle.sh)
    return(res)
}

#' Simulation of Independent Bernoulli Observations
#' @description Simulation of a sequence of independent Bernoulli
#' Observations. To reduce the amount of random draws, each
#' simulation is based on a sequence of standard normal
#' variables, and whether each observation is above a shift
#' defined by the binomial probabilities assumed.
#' @param nser length of sequence simulated
#' @param nsim number of simulations
#' @param  probs binomial probabilites
#' @return a data frame with the number of crossings and
#' longest run for each probability. For instance
#' the variables nc0.5 and lr0.5 are the number of
#' crossings and the longest run for success probability
#' 0.5. One row for each simulation.
#' @examples
#' cl30simbin <- simclbin(nser=30, nsim=100)
#' mean(cl30simbin$nc0.5) # mean number of crossings, p=0.5
#' mean(cl30simbin$lr0.9) # mean longest run, p=0.9
#' @export
simclbin <- function(nser = 100, nsim = 1e+05, probs = c(0.5, 0.6, 0.7,
    0.8, 0.9)) {
    nprob <- length(probs)
    shifts <- stats::qnorm(probs)
    series <- data.frame(matrix(stats::rnorm(nser * nsim), nrow = nsim))
    res <- data.frame(matrix(rep(NA, 2 * nsim * nprob), nrow = nsim))
    names(res) <- paste(rep(c("nc", "lr"), nprob), rep(probs, rep(2,
        nprob)), sep = "")
    for (shnr in 1:nprob) {
        res[, 1 + 2 * (shnr - 1)] <- apply(series, 1, clshift, shift = shifts[shnr],
            type = 0)
        res[, 2 * shnr] <- apply(series, 1, clshift, shift = shifts[shnr],
            type = 1)
    }
    return(res)
}

#' Simulation of observations in n=2m independent and
#' identical continuous variables classified as above or
#' below the median.
#' @description Simulation of observations in n=2m independent and
#' identical continuous variables classified as above or
#' below the median. The observations are drawn from a standard
#' normal distruibution, but the results do not depend on the
#' distribution as long as it is continuous.
#' @param m1, half the number of observations
#' @param nsim number of simulations
#' @return a data frame with the number of crossings cs and
#' longest run ls in each simulation.
#' @examples
#' simcl14 <- simclem(m1=14, nsim=100)
#' mean(simcl14$cs) # mean number of crossings
#' mean(simcl14$ls) # mean longest run
#' @export
simclem <- function(m1=14, nsim = 1e+05) {
  n1 <- 2*m1
  cs <- rep(NA, nsim)
  ls <- rep(NA, nsim)
  for (sim in 1:nsim) {
    series <- stats::rnorm(n1)
    med <- stats::median(series)
    above <- as.numeric(series>med)
    rleabove <- rle(above)$lengths
    cs[sim] <- length(rleabove) - 1
    ls[sim] <- max(rleabove)
  } # end for sim
  return(data.frame(cs=cs,ls=ls))
} #end function simclem

