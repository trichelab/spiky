#' A handful of methods that I've always felt were missing
#'
#' Particularly, simple methods to plot coverage results.
#' 
#' @param  x                an Rle or RleList, usually
#' @param  y                not usedan Rle or RleList, usually
#' @param  ...              other params such as `ylim` passed to `barplot`
#' 
#' @return                  invisibly, the plot details
#' 
#' @details
#' selectMethod("plot", "Rle") and also
#' selectMethod("plot", "RleList") too.
#' 
#' @importFrom graphics barplot par text title
#' 
#' @rdname spiky-methods
#' @name   spiky-methods
NULL


#' @rdname spiky-methods
#' @export 
setMethod("plot", "Rle", 
          function(x, y, ...) {
            titlex <- sum(width(x)) * 0.9
            titley <- ifelse(exists("ylim"), ylim[2] * 0.9, max(x) * 0.9)
            barplot(runValue(x), width(x), las=2, inside=FALSE, col="blue", ...)
            if (exists("title")) text(titlex, titley, title)
          })


#' @rdname spiky-methods
#' @export 
setMethod("plot", "SimpleRleList", 
          function(x, y, ...) {
            yl <- c(0, max(max(x)))
            idx <- which(!all(x == 0))
            par(mfrow=c(length(idx), 1))
            par("mai"=c(0.05, 0.5, 0.05, 0.05))
            if (is.null(names(idx))) names(idx) <- paste0("run", idx)
            for (ni in names(idx)) plot(x[[idx[ni]]], ylim=yl, title=ni)
          })
