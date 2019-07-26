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
  fvtrj <- read.table(paste0(dir, list.files(dir, pattern="^fv-traj-")))
  trj <- fvtrj[,c(2,3)]
  trja <- read.table(paste0(dir, list.files(dir, pattern="^av-traj-")), colClasses = "character")
  trjancol <- ncol(trja)
  
  # Read trajactory coordinates from av-trj file
  if(trjancol==3|trjancol==4){
    trja <- trja[,2:3]
  }else if(trjancol==6|trjancol==7){
    trja <- trja[,c(2,3,5,6)]
  }else {
    stop("Unknown av-trj format")
  }
  
  trja <- as.data.frame(sapply(trja,gsub,pattern="\\[",replacement=""), stringsAsFactors=F)
  trja <- as.data.frame(sapply(trja,gsub,pattern="\\]",replacement=""), stringsAsFactors=F)
  trja <- data.frame(sapply(trja, as.numeric)) # trajectory is in pixel coordinate
  data(map, package = "Flyception2R") # load the pixel to degree map
  
  # Offset needs to be taken into account
  #define XVOLTPERDEGREE 0.55
  #define YVOLTPERDEGREE 0.525
  #define XOFFSET -0.25
  #define YOFFSET -0.315
  #Therefore x offset in degree is -0.592 and y offset is -1.192
  map[,1] <- map[,1] - 0.592
  map[,2] <- map[,2] + 1.192
  
  
  ## Creating a finer map by loess fit
  names(map)[1] <- "ax"
  names(map)[2] <- "ay"
  names(map)[3] <- "x"
  names(map)[4] <- "y"
  mapcl <- map[which(map[,3]!=0),]
  # create x surface
  mapfitx <- loess(ax ~ y*x, mapcl, degree=2, span=0.25, normalize=F)
  mappoints <- list(x=seq(1, 512, 1), y=seq(1, 512, 1))
  mapsurfacex <- predict(mapfitx, expand.grid(mappoints), se=F)
  mapsurfacesx <- list(mappoints$x, mappoints$y,
                       matrix(mapsurfacex, length(mappoints$x), length(mappoints$y)))
  names(mapsurfacesx) <- c("x", "y", "z")
  # create y surface
  mapfity <- loess(ay ~ y*x, mapcl, degree=2, span=0.25, normalize=F)
  mapsurfacey <- predict(mapfity, expand.grid(mappoints), se=F)
  mapsurfacesy <- list(mappoints$y, mappoints$y,
                       matrix(mapsurfacey, length(mappoints$x), length(mappoints$y)))
  names(mapsurfacesy) <- c("x", "y", "z")
  
  get_world_traj_from_angle <- function(xangle, yangle) {
    a     <- 15.174 # distance between x and y mirrors
    b     <- 68.167 # distance between y center and arena surface at the center
    c     <- 65.167 # distance between y center and top flat surface
    alpha <- pi*xangle/180 # x mirror angle in radians
    beta  <- pi*yangle/180 # y mirror angle in radians

    hr <- (a+b)*cos(2*alpha) - a
    xr <- a*tan(2*alpha) + hr*tan(2*alpha)
    yr <- hr*sin(2*beta)
    
    return(cbind(xr,yr))
  }
  
  
  if(interaction == F){
    
    # Use flyview trajectory for the tracked fly's angle
    xanglefv <- fvtrj[,6]
    yanglefv <- fvtrj[,7]
    
    # Lookup fly1 angle from arena view interoplated map using arena trajectory
    xangle1 <- array(NA,length(trja[,1]))
    xangle1[which(trja[,1] != 0)] <- mapsurfacesx$z[as.matrix(round(trja[which(trja[,1] != 0),1:2]))]
    yangle1 <- array(NA,length(trja[,2]))
    yangle1[which(trja[,2] != 0)] <- mapsurfacesy$z[as.matrix(round(trja[which(trja[,2] != 0),1:2]))] 
    
    trjfvmm <- get_world_traj_from_angle(xanglefv,yanglefv)
    trjamm <- get_world_traj_from_angle(xangle1,yangle1)

    
  } else {
    
    # TODO: Save both fv and av traj for nn identity tracking.
    
    
    # Use flyview trajectory for the tracked fly
    xanglefv <- fvtrj[,6]
    yanglefv <- fvtrj[,7]
    
    # Lookup fly1 angle from arena view interoplated map using arena trajectory
    xangle1 <- array(NA,length(trja[,1]))
    xangle1[which(trja[,1] != 0)] <- mapsurfacesx$z[as.matrix(round(trja[which(trja[,1] != 0),1:2]))]
    yangle1 <- array(NA,length(trja[,2]))
    yangle1[which(trja[,2] != 0)] <- mapsurfacesy$z[as.matrix(round(trja[which(trja[,2] != 0),1:2]))] 
    
    #  Lookup fly2 angle from arena view interoplated map using arena trajectory
    xangle2 <- array(NA,length(trja[,3]))
    xangle2[which(trja[,3] != 0)] <- mapsurfacesx$z[as.matrix(round(trja[which(trja[,3] != 0),3:4]))]
    yangle2 <- array(NA,length(trja[,4]))
    yangle2[which(trja[,4] != 0)] <- mapsurfacesy$z[as.matrix(round(trja[which(trja[,4] != 0),3:4]))]  
    
    trjfvmm <- get_world_traj_from_angle(xanglefv,yanglefv)
    trjamm <- get_world_traj_from_angle(xangle1,yangle1)
    trjamm <- cbind(trjamm,get_world_traj_from_angle(xangle2,yangle2))
    
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
    flydist <- rep(0, nrow(trjamm))
  }
  return(list("speed"=speed, "error"=error, "trja"=trjamm, "trjfv"=trjfvmm, "flydist"=flydist))
  
}
