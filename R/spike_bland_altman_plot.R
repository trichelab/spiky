#' Bland-Altman plot for cfMeDIP spike standards
#'
#' @param   fit             a model fit, from predict_pmol (?)
#'
#' @return                  a ggplot2 object
#'
#' @import  BlandAltmanLeh
#' @import  ggplot2
#' @import  scales
#'
#' @examples
#'
#' data(spike_res)
#' data(spike, package="spiky")
#' fit <- model_glm_pmol(covg_to_df(spike_res, spike=spike),spike=spike)
#' ba_plot <- spike_bland_altman_plot(fit)
#'
#' @export
spike_bland_altman_plot <- function(fit) {

  if (!is.data.frame(fit$x)) fit$x <- attr(fit, "data")
  fit$x$pred_conc <- predict(fit, fit$x)

  fit$x$fraglen <- as.factor(paste0(fit$x$fraglen, "bp"))
  BA <- bland.altman.plot(fit$x$conc,
                          fit$resid,
                          conf.int = .95,
                          col.points = fit$x$fraglen,
                          pch = 19,
                          graph.sys = "ggplot2")

  BA <- BA +
    aes(color = fit$x$fraglen,shape = fit$x$fraglen) +
    theme_bw(base_line_size = 0.25) +
    scale_color_manual(values = c("80bp" = "black",
                                  "160bp" = "darkgray",
                                  "320bp" = "navy")) +
    theme(legend.title = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 14),
          axis.title.y = element_text(size = 20),
          axis.text.x  = element_text(size = 16),
          axis.text.y = element_text(vjust = 0.5, size = 16),
          axis.title.x = element_text(size = 20),
          panel.grid.minor.x = element_blank()) +
    xlab("Mean of measurements (picomoles)") +
    ylab("Difference (picomoles)")
  return(BA)

}
