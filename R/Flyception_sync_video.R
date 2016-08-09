sync_frame <- function(rdir, dir, prefix, reuse=T, fpsfl=100){
  require(RImageBook)
  require(Rcpp)
  require(zoo)
  require(RNiftyReg)
  require(plotrix)
  require(ggplot2)
  require(plyr)
  require(reshape2)
  require(rlogging)
  source(paste0(rdir, "sfeatures.R"))
  source(paste0(rdir, "rollmeanimg.R"))
  source(paste0(rdir, "rollmedianimg.R"))
  source(paste0(rdir, "sweepC.R"))
  source(paste0(rdir, "readFMF2.R"))
  source(paste0(rdir, "plotter.R"))
  source(paste0(rdir, "drawtext.R"))
  source(paste0(rdir, "fmf2tif.R"))
  source(paste0(rdir, "pseudoColor3.R"))
  
  # Start logging
  SetLogFile(base.file=paste0(prefix, "_sync_frame_log.txt"), folder=dir)
  message(dir)
  
  # Load fluorescence movie and detect flash
  flfile <- paste0(dir, list.files(dir, pattern="ome\\.tif$"))
  message(sprintf("Reading %s", flfile))
  flimg <- readImage(flfile)
  message("Detecting flash in flview")
  flimgint <- colMeans(flimg, dim=2)   
  png(file=paste0(dir, prefix, "_flflash.png"), width=400, height=400)
  plot(flimgint)
  dev.off()
  flimgintdif <- diff(flimgint)
  flflashes <- which(flimgintdif > 0.005) + 1
  flimgflash <- min(flflashes)
  if(flimgflash==Inf) stop("Flash was not detected in flcamera.")
  nframesfl <- dim(flimg)[3] 
  message(paste0("Number of frames in flview: ", nframesfl))
  flimgrt <- rotate(flip(flimg), -90)
  rm(flimg)
  message(sprintf("Flash was detected in fluo-view frames: %s", paste(flflashes, collapse=" ")))
  
  # Detect flash in flyview
  message("Detecting flash in flyview")
  if(file.exists(paste0(dir, prefix, "_fvimgsubint.RDS"))==T & reuse==T){
    message("Loading from RDS file")
    fvimgsubint <- readRDS(paste0(dir, prefix, "_fvimgsubint.RDS"))
  }else{
    fvimgsub1 <- readFMF(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), crop=c(5,10,5,10))
    fvimgsub2 <- readFMF(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), crop=c(220,225,220,225))
    fvimgsubint1 <- colMeans(fvimgsub1, dim=2)
    fvimgsubint2 <- colMeans(fvimgsub2, dim=2)
    rm(fvimgsub1)
    rm(fvimgsub2)
    fvimgsubint <- ifelse(fvimgsubint1 > fvimgsubint2, fvimgsubint1, fvimgsubint2)
    png(file=paste0(dir, prefix, "_fvflash.png"), width=400, height=400)
    plot(fvimgsubint)
    dev.off()
    saveRDS(fvimgsubint, file=paste0(dir, prefix, "_fvimgsubint.RDS"))
  }
  
  fvimgflash <- min(which(fvimgsubint > 135))
  if(fvimgflash==Inf) stop("Flash was not detected in fvcamera.")
  message(sprintf("Flash was detected in fly-view frames: %s", paste(which(fvimgsubint > 135), collapse=" ")))

  # Detect flash in arenaview
  message("Detecting flash in arenaview")
  if(file.exists(paste0(dir, prefix, "_avimgsubint.RDS"))==T & reuse==T){
    message("Loading from RDS file")
    avimgsubint <- readRDS(paste0(dir, prefix, "_avimgsubint.RDS"))
  }else{
    avimgsub <- readFMF(paste0(dir, list.files(dir, pattern="^av.*fmf$")), crop=c(5,10,5,10))
    avimgsubint <- colMeans(avimgsub, dim=2)
    rm(avimgsub)
    png(file=paste0(dir, prefix, "_avflash.png"), width=400, height=400)
    plot(avimgsubint)
    dev.off()
    saveRDS(avimgsubint, file=paste0(dir, prefix, "_avimgsubint.RDS"))
  }
  avimgflash <- min(which(avimgsubint > 100))
  if(avimgflash==Inf) stop("Flash was not detected in avcamera.")
  message(sprintf("Flash was detected in Arena-view frames: %s", paste(which(avimgsubint > 100), collapse=" ")))

  # Start time
  message("Analyzing metadata")
  metadata <-scan(paste0(dir, list.files(dir, pattern="metadata\\.txt$")), what=character(),sep="")
  log <- scan(paste0(dir, list.files(dir, pattern="fv-log-")), what=character(),sep="")
  avlog <- scan(paste0(dir, list.files(dir, pattern="av-log-")), what=character(),sep="")
  starttimefl <- metadata[which(metadata == "Time")[1]+2]
  
  # Exposure and binning
  exposure <- substr(metadata[which(metadata == "Exposure-ms")[1]+2], 1, nchar(metadata[which(metadata == "Exposure-ms")[1]+2])-1)
  binning <- metadata[which(metadata == "Binning")[1]+2]
  message(sprintf("flview exposure: %s", exposure))
  message(sprintf("flview binning: %s", binning))
  
  # Elapsed time (in ms) of each frame from the fluorescence camera relative to the flash
  elapsedtimefl <- metadata[grep("ElapsedTime-ms", metadata)+2]
  elapsedtimefl <- as.numeric(substr(elapsedtimefl, 1, nchar(elapsedtimefl)-1))
  #fpsfl <- round(length(elapsedtimefl)/tail(elapsedtimefl, n=1)*1000, 2)
  message(paste("fluorescence camera:", fpsfl, "fps"))
  elapsedtimeflflash <- elapsedtimefl - elapsedtimefl[flimgflash]
  elapsedtimeflflashdiff <- diff(elapsedtimeflflash)
  
  # Elapsed time (in ms) of each frame from the fly view camera
  timestampusec <- as.numeric(log[grep("TimeStamp", log)+1])
  elapsedtimefv <- (timestampusec - timestampusec[1])/1000
  elapsedtimefvflash <- elapsedtimefv - elapsedtimefv[fvimgflash]
  fpsfv <- round(length(elapsedtimefv)/((tail(elapsedtimefv, n=1) - elapsedtimefv[1])/1000), 2)
  message(paste("flyview camera:", fpsfv, "fps")) 
  elapsedtimefvflashdiff <- diff(elapsedtimefvflash)
  # Only used for plotting purpose
  elapsedtime <- elapsedtimeflflash
  
  # Elapsed time (in ms) of each frame from the arenaview camera
  avtimestampcyclesec <- avlog[grep("TimeStamp", avlog)+1]
  avtimestampcyclesec <- as.numeric(avtimestampcyclesec)
  cnt <- 0
  avtimestampsec <- rep(0, length(avtimestampcyclesec))
  avtimestampsec[1] <- avtimestampcyclesec[1]
  for(t in 2:length(avtimestampsec)){
    if(avtimestampcyclesec[t-1]==127 & avtimestampcyclesec[t]==0) cnt <- cnt + 1
    avtimestampsec[t] <- avtimestampcyclesec[t] + 127*cnt
  }
  avtimestampsec <- avtimestampsec*1000
  avtimestampcnt <- avlog[grep("TimeStamp", avlog)+2]
  avtimestampmsec <- 1/8000*as.numeric(avtimestampcnt)*1000
  elapsedtimeav <- (avtimestampsec + avtimestampmsec) - (avtimestampsec[1] + avtimestampmsec[1])
  elapsedtimeavflash <- elapsedtimeav - elapsedtimeav[avimgflash]
  fpsav <- round(length(elapsedtimeav)/((tail(elapsedtimeav, n=1) - elapsedtimeav[1])/1000), 2)
  message(paste("arenaview camera:", fpsav, "fps")) 
  
  # Align frames between flyview and flview cameras
  message("Aligning frames between flyview and flview")
  if(file.exists(paste0(dir, prefix, "_frid.RDS"))==T & reuse==T){
    message("Loading RDS file")
    frid <- readRDS(paste0(dir, prefix, "_frid.RDS"))
  }else{
    frameratio <- round(fpsfv/fpsfl)
    message(paste0("fv/fl frame ratio: ", frameratio))
    # Hypothetical trigger
    frid <- seq(fvimgflash-(flimgflash-1)*frameratio, 
                fvimgflash+frameratio*(nframesfl-flimgflash), frameratio)
    # Interval between frames in ms
    framediff <- elapsedtimefvflashdiff
    # Convert ms to trigger count
    framediff[which(framediff<2)] <- 1
    framediff[which(framediff>2 & framediff<3)] <- 2
    framediff[which(framediff>3)] <- 3
    # Generate sequence of triggers for each frame
    fvfrsum <- cumsum(framediff)
    froffset <- fvfrsum[fvimgflash] - fvimgflash
    fvfrsum <- fvfrsum - froffset
    # Find matching and nearest frames to the hypothetical trigger
    frid2 <- match(frid, fvfrsum)
    frid2[which(is.na(frid2))] <- sapply(frid[which(is.na(frid2))], function(x) which.min(abs(fvfrsum-x)))
    frid <- unlist(frid2)
    # Check if two flashes match
    fvflashesfrid <- frid[flflashes]
    fvflashes <- which(fvimgsubint > 135)
    message(sprintf("Hypothetical flash for Fly-view: %s", paste(fvflashesfrid, collapse=" ")))
    message(paste0("Did flash frames match in flyview and flview? ", all.equal(fvflashesfrid, fvflashes)))
    if(length(flflashes)!=length(fvflashes)) {
      warning("Number of flash detected did not match between flyview and flview.")
    }
    # Load flyview camera images
    fvimgmaxfr <- readFMF(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), getFrames=T)
    frid <- frid[which(frid<=fvimgmaxfr)]
    if(length(frid) < nframesfl){
      frid <- c(frid, rep(tail(frid, 1), nframesfl-length(frid)))
    }
    saveRDS(frid, paste0(dir, prefix, "_frid.RDS"))
  }
  fvimgl <- readFMF(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), frames=frid)
  nframesfv <- dim(fvimgl)[3]
  
  # Align frames between flview and arenaview cameras
  message("Aligning frames between flview and arenaview")
  if(file.exists(paste0(dir, prefix, "_frida.RDS")) & reuse==T){
    message("Loading RDS file")
    frida <- readRDS(paste0(dir, prefix, "_frida.RDS"))
  }else{
    frameratio2 <- round(fpsav/fpsfl)
    message(paste0("av/fl frame ratio: ", frameratio2))
    # Hypothetical trigger
    frida <- seq(avimgflash-(flimgflash-1)*frameratio2, 
                 avimgflash+frameratio2*(nframesfl-flimgflash), frameratio2)
    # Check if two flashes match
    avflashesfridav <- frida[flflashes]
    avflashes <- which(avimgsubint > 100)
    message(sprintf("Hypothetical flash for Arena-view: %s", paste(avflashesfridav, collapse=" ")))
    message(paste0("Did flash frames match in arenaview and flview? ", all.equal(avflashesfridav, avflashes)))
    if(length(flflashes)!=length(avflashes)) {
      message("Number of flash detected did not match between areanview and flfiew.")
    }
  }
  
}