#' Detect flash
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' detect_flash()
#'

detect_flash <- function(input, output, type=c("fluo", "fly", "arena"), flash_thresh, reuse=F){
  ## TODO
  ## Adaptive ROI

  if(type=="fluo"){
    if(file.exists(paste0(output, "_flimgint.RDS"))==T & reuse==T){
      message("Loading from RDS file")
      flimgint <- readRDS(paste0(output, "_flimgint.RDS"))
    } else{
    message(sprintf("Reading %s", input))
    message(sprintf("Flash thresh is %f", flash_thresh))
    flimgint <- readTIFF2(input, intensity=T)
    png(file=paste0(output, "_flflash.png"), width=400, height=400)
    plot(flimgint)
    dev.off()
    nframesfl <- dim(flimg)[3]
    message(paste0("Number of frames in fluo-view: ", nframesfl))
    rm(flimg)
    }

    flflashes <- which(flimgint > flash_thresh)
    flimgflash <- min(flflashes)
    if(flimgflash==Inf) stop("Flash was not detected in fluo-view.")
    message(sprintf("Flash was detected in fluo-view frames: %s", paste(flflashes, collapse=" ")))
    return(list("flflashes"=flflashes, "nframesfl"=nframesfl))
  }

  # Detect flash in fly-view
  if(type=="fly"){
    if(file.exists(paste0(output, "_fvimgsubint.RDS"))==T & reuse==T){
      message("Loading from RDS file")
      fvimgsubint <- readRDS(paste0(output, "_fvimgsubint.RDS"))
    } else{
      message(sprintf("Reading %s", input))
      message(sprintf("Flash thresh is %f", flash_thresh))
      # Load only diagonal ROIs
      nframesfv <- dipr::readFMF(input, getFrames=T)
      fvimgsub1 <- dipr::readFMF(input, crop=c(5,10,5,10))
      fvimgsub2 <- dipr::readFMF(input, crop=c(220,225,220,225))
      fvimgsubint1 <- colMeans(fvimgsub1, dim=2)
      fvimgsubint2 <- colMeans(fvimgsub2, dim=2)
      rm(fvimgsub1)
      rm(fvimgsub2)
      fvimgsubint <- ifelse(fvimgsubint1 > fvimgsubint2, fvimgsubint1, fvimgsubint2)
      png(file=paste0(output, "_fvflash.png"), width=400, height=400)
      plot(fvimgsubint)
      dev.off()
      saveRDS(fvimgsubint, file=paste0(output, "_fvimgsubint.RDS"))
    }

    fvflashes <- which(fvimgsubint > flash_thresh)
    fvimgflash <- min(fvflashes)
    if(fvimgflash==Inf) stop("Flash was not detected in fly-view.")
    message(sprintf("Flash was detected in fly-view frames: %s", paste(fvflashes, collapse=" ")))
    return(list("fvflashes"=fvflashes, "nframesfv"=nframesfv))
  }

  # Detect flash in arena-view
  if(type=="arena"){
    if(file.exists(paste0(output, "_avimgsubint.RDS"))==T & reuse==T){
      message("Loading from RDS file")
      avimgsubint <- readRDS(paste0(output, "_avimgsubint.RDS"))
     }else{
      message(sprintf("Reading %s", input))
      message(sprintf("Flash thresh is %f", flash_thresh))
      nframesav <- dipr::readFMF(input, getFrames=T)
      avimgsub <- dipr::readFMF(input, crop=c(5,10,5,10))
      avimgsubint <- colMeans(avimgsub, dim=2)
      rm(avimgsub)
      png(file=paste0(output, "_avflash.png"), width=400, height=400)
      plot(avimgsubint)
      dev.off()
      saveRDS(avimgsubint, file=paste0(output, "_avimgsubint.RDS"))
    }

    avflashes <- which(avimgsubint > flash_thresh)
    avimgflash <- min(avflashes)
    if(avimgflash==Inf) stop("Flash was not detected in arena-view.")
    message(sprintf("Flash was detected in arena-view frames: %s", paste(avflashes, collapse=" ")))
    return(list("avflashes"=avflashes, "nframesav"=nframesav))
  }
}
