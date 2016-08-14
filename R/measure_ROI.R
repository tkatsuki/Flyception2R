#' Create delta F over F plot for ROI
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' measureROI()
#'

measureROI <- function(img, x, y, w, h){
  ROI <- img[x:(x+w-1),y:(y+h-1),]
  meanint <- colMeans(ROI, dim=2)
  ROIF0 <- mean(ROI[1:5])
  deltaROIF <- ROI - ROIF0
  ROIdFF0 <- deltaROIF/ROIF0 * 100
  dat <- data.frame(x=wFfr, y=ROIdFF0, d=flydist[frida[Ffr]])
  p <- ggplot(data=dat, aes(x=x, y=y)) +
    geom_smooth(method="loess", span = 0.4, level=0.95) +
    ylim(-5, 10) +
    geom_line(data=dat, aes(x=x, y=d))
  ggsave(filename = paste0(dir, prefix, "_dFF0int_ROI_", ffr, ".pdf"), width = 8, height = 8)
  return(meanint)

}
