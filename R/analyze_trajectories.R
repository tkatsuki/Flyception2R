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
  trja <- data.frame(sapply(trja, as.numeric)) # trajectory is in pixel coordinate
  data(map, package = "Flyception2R") # load the pixel to degree map

  if(interaction == F){
    matchedrow <- sapply(1:nrow(trja), function(x) which.min(abs((map[,3]-trja[x, 1])^2) + (map[,4]-trja[x, 2])^2))
    trjadeg <- map[matchedrow, 1:2]
    trjaalpha <- pi*trjadeg[,1]/180
    a <- 15.174 # distance between x and y mirrors
    b <- 68.167 # distance between y center and arena surface at the center
    c <- 65.167 # distance between y center and top flat surface
    trjahr <- (a+b)*cos(2*trjaalpha) - a
    trjaxr <- a*tan(2*trjaalpha) + trjahr*tan(2*trjaalpha) # height should be h not b
    trjabeta <- pi*trjadeg[,2]/180
    trjyfr <- function(alpha, beta) ((a+b)*cos(2*alpha) - a)*sin(2*beta)
    trjayr <- trjyfr(trjaalpha, trjabeta)
    trjamm <- cbind(trjaxr, trjayr) # trajectory in mm
  } else {
    matchedrow <- sapply(1:nrow(trja), function(x) which.min(abs((map[,3]-trja[x, 1])^2) + (map[,4]-trja[x, 2])^2))
    trjadeg <- map[matchedrow, 1:2]
    trjaalpha <- pi*trjadeg[,1]/180
    a <- 15.174 # distance between x and y mirrors
    b <- 68.167 # distance between y center and arena surface at the center
    c <- 65.167 # distance between y center and top flat surface
    trjahr <- (a+b)*cos(2*trjaalpha) - a
    trjaxr <- a*tan(2*trjaalpha) + trjahr*tan(2*trjaalpha) # height should be h not b
    trjabeta <- pi*trjadeg[,2]/180
    trjyfr <- function(alpha, beta) ((a+b)*cos(2*alpha) - a)*sin(2*beta)
    trjayr <- trjyfr(trjaalpha, trjabeta)
    trjamm <- cbind(trjaxr, trjayr) # trajectory in mm
    
    matchedrow <- sapply(1:nrow(trja), function(x) which.min(abs((map[,3]-trja[x, 3])^2) + (map[,4]-trja[x, 4])^2))
    trjadeg <- map[matchedrow, 1:2]
    trjaalpha <- pi*trjadeg[,1]/180
    a <- 15.174 # distance between x and y mirrors
    b <- 68.167 # distance between y center and arena surface at the center
    c <- 65.167 # distance between y center and top flat surface
    trjahr <- (a+b)*cos(2*trjaalpha) - a
    trjaxr <- a*tan(2*trjaalpha) + trjahr*tan(2*trjaalpha) # height should be h not b
    trjabeta <- pi*trjadeg[,2]/180
    trjyfr <- function(alpha, beta) ((a+b)*cos(2*alpha) - a)*sin(2*beta)
    trjayr <- trjyfr(trjaalpha, trjabeta)
    trjamm <- cbind(trjamm, trjaxr, trjayr) # trajectory in mm
  }
  
  headpos <- fvtrj[,c(4,5)]
  distance <- dipr::trackDistance(trj)
  distance <- zoo::rollmedian(distance, k=5)
  speed <- zoo::rollsum(distance, k=200)/200*fpsfv
  error <- sqrt(diff(headpos[,1])^2 + diff(headpos[,2])^2)
  png(file=paste0(output, "_error.png"), width=400, height=400)
  plot(error)
  dev.off()
  
  if(interaction==T){
    flydist <- sqrt((trjamm[,1] - trjamm[,3])^2 + (trjamm[,2] - trjamm[,4])^2)
    png(file=paste0(output, "_flydist.png"), width=400, height=400)
    plot(flydist, ylab="Distance between flies (mm)")
    dev.off()
  } else {
    flydist <- NA
  }
  return(list("speed"=speed, "error"=error, "trja"=trjamm, "flydist"=flydist))
  
}
