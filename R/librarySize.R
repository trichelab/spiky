#' adapted from code in Caleb Lareau's mtscATAC paper, which adapts Picard's 
#' 
#' @param nT      total number of genomic fragments 
#' @param nU      unique genomic fragments 
#' 
#' @return        estimated number of fragments in the source library
#' 
#' @import        Rsamtools 
#' 
#' @export 
librarySize <- function(nT, nU) {

  f <- function(x, c, n) return(c / x - 1 + exp( -n / x ))

  m = 1
  M = 100

  nDup <- (nT - nU) + 1 # charity to handle only unique reads observed

  if (nU > nT | (f(m*nU, nU, nT) < 0) | nU < 0 | nT < 0 | nDup < 0) {
    message("Library size returns 0 -- invalid inputs")
    message("Check this cell more closely!")
    return(0)
  }

  while (f(M * nU, nU, nT) > 0) M <- M*10.0

  # shady?
  for(i in seq(0, 40)) {
    r <-  (m + M) / 2.0
    u <- f(r * nU, nU, nT);
    if (u == 0) {
      break
    } else if (u > 0) {
      m = r
    } else if (u < 0) {
      M = r
    }
  }

  return(round(nU * (m + M) / 2.0))

}

# todo: add a wrapper to do this for spikes alongside genomic contigs.
# todo: obtain a correction term that equilibrates complexity ests by spikes.
# todo: use the resulting complexity ratio and the spike frag ratio to equalize
