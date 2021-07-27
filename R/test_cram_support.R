#' a function to determine if CRAM support is enabled
#' 
#' At present, CRAM support is not enabled by default in Rsamtools. To fix this,
#' it may be necessary to use BiocManager::install() on GitHub branches, i.e. 
#'   BiocManager::install(c(
#'     "Bioconductor/Rhtslib@enable-cram", 
#'     "Bioconductor/Rsamtools@cram"
#'   ))
#' 
#' Until such time as CRAM support is included by default, though, we need to
#' know whether access to CRAMs is available.  This function tests that.  If 
#' the function detects that CRAM support is NOT installed, it instructs the 
#' user on the best known way to enable it. (Ultimately, we want it in RELEASE)
#' 
#' @param   cram    a CRAM filename (default is example.spike.cram in extdata)
#' @param   reads   how many reads to fetch, i.e. the yieldSize (1000)
#' 
#' @return  TRUE or FALSE, perhaps with a message suggesting how to make it TRUE
#'
#' @examples
#' cram <- system.file("extdata", "example.spike.cram", package="spiky",
#'                   mustWork=TRUE)
#' \dontrun{
#'   test_cram_support(cram)
#' }
#' 
#' @import GenomicAlignments
#' 
#' @export 
test_cram_support <- function(cram=NULL, reads=1000L) { 

  if (is.null(cram)) {
    cram <- system.file("extdata", "example.spike.cram", package="spiky", 
                        mustWork=TRUE)
  }

  cf <- BamFile(cram, yieldSize=reads, asMates=TRUE)
  gal <- try(readGAlignments(cf)) 
  CRAMsupport <- (length(gal) == reads)
  if (CRAMsupport) message("You seem to have CRAM support.")
  if (!CRAMsupport) message("You do not seem to have CRAM support.")
  if (!CRAMsupport) message("Try the instructions in ?test_cram_support.")
  return(CRAMsupport)

}
