#' Analyze trajectory of flies
#'
#'
#' @param dir Input directory.
#' @param output Output path for saving results.
#' @param fpsfv fly-view frame fate.
#' @param interaction logical. If True, measures distance between two flies.
#' @export
#' @examples
#' analyze_trajectories()
#'

analyze_trajectories <- function(dir, output, fpsfv, interaction=F){
  fvtrj <- read.table(paste0(dir, list.files(dir, pattern="fv-traj-")))
  trj <- fvtrj[,c(2,3)]
  trja <- read.table(paste0(dir, list.files(dir, pattern="av-traj-")), colClasses = "character")
  trjancol <- ncol(trja)
  trja <- trja[,2:trjancol]
  trja <- as.data.frame(sapply(trja,gsub,pattern="\\[",replacement=""), stringsAsFactors=F)
  trja <- as.data.frame(sapply(trja,gsub,pattern="\\]",replacement=""), stringsAsFactors=F)
  trja <- data.frame(sapply(trja, as.numeric))
  data(map, package = "Flyception2R")
  trjaint <- round(trja)
  matchedrow <- row.match(trjaint[,1:2], map[,3:4])
  matchedrow[which(is.na(matchedrow))] <- sapply(which(is.na(matchedrow)), function(x) which.min(abs((map[,3]-trjaint[x, 1])^2) + (map[,4]-trjaint[x, 2])^2))
  trjav <- map[matchedrow, 1:2]
  
  headpos <- fvtrj[,c(4,5)]
  distance <- dipr::trackDistance(trj)
  distance <- zoo::rollmedian(distance, k=5)
  speed <- zoo::rollsum(distance, k=200)/200*fpsfv
  error <- sqrt(diff(headpos[,1])^2 + diff(headpos[,2])^2)
  png(file=paste0(output, "_error.png"), width=400, height=400)
  plot(error)
  dev.off()

  if(interaction==T){
    flydist <- sqrt((trja[,1] - trja[,3])^2 + (trja[,2] - trja[,4])^2)
    png(file=paste0(output, "_flydist.png"), width=400, height=400)
    plot(flydist, ylab="Distance between flies (mm)")
    dev.off()
  } else {
    flydist <- NA
  }
  return(list("speed"=speed, "error"=error, "trja"=trja, "flydist"=flydist))
  
}
