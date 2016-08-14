#' Plot delta F over F
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' plot_dF_F0()
#'
plot_dF_F0 <- function(intensity, goodfr){
  Fintfr <- fridstim[fintfr]:(fridstim[fintfr]+stimplotlen-1)
  Fintfr <- Fintfr[which(Fintfr < nframesfl)]
  Fintfrg <- intersect(goodfr, Fintfr)
  Fintfr[!Fintfr%in%goodfr] <- NA
  F0int <- mean(intensity[Fintfrg[1:5]])
  deltaFint <- intensity[Fintfr] - F0int
  dFF0int <- deltaFint/F0int * 100
  dat <- data.frame(x=(1:length(dFF0int)), y=dFF0int, d=flydist[Fintfr])
  p <- ggplot(data=dat, aes(x=x, y=y)) +
    geom_smooth(method="loess", span = 0.4, level=0.95) +
    ylim(-5, 10) +
    geom_line(data=dat, aes(x=x, y=d))
  ggsave(filename = paste0(dir, prefix, "_dFF0int_", fintfr, ".pdf"), width = 8, height = 8)

}
