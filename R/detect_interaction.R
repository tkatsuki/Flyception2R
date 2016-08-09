#' Detect fly-fly interactions
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' detect_interaction()
#'

detect_interaction <- function(){
  # If flies are interacting
  if(interacting==T){
    flydist <- sqrt((trja[,1] - trja[,3])^2 + (trja[,2] - trja[,4])^2)
    png(file=paste0(dir, prefix, "_flydist.png"), width=400, height=400)
    plot(flydist, ylab="Distance between flies (mm)")
    dev.off()
    closefr <- which(flydist < dist_thresh)
    interaction <- c(closefr[1], closefr[which(diff(closefr)>stimlen)+1])
    if(is.na(stimfr)){
      fridstim <- NA
    } else {
      fridstim <- sapply(stimulusa, function(x) which.min(abs(frida-x)))
      fridstim <- fridstim - 100
      fridstim[which(fridstim < 1)] <- 1
    }
  } else {
    flydist <- rep(0, nframesfl)
  }
}
