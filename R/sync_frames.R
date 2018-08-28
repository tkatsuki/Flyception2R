#' Synchronize camera frames
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' sync_frames()

sync_frames <- function(dir, fluo_flash, fly_flash, arena_flash, output, reuse=F){

  # Start time
  message("Analyzing metadata")
  metadata <-scan(paste0(dir, list.files(dir, pattern="metadata\\.txt$")), what=character(),sep="")
  log <- scan(paste0(dir, list.files(dir, pattern="fv-log-")), what=character(),sep="")
  avlog <- scan(paste0(dir, list.files(dir, pattern="av-log-")), what=character(),sep="")
  starttimefl <- metadata[which(metadata == "Time")[1]+2]

  # Exposure and binning
  exposure <- substr(metadata[which(metadata == "Exposure-ms")[1]+2], 1, nchar(metadata[which(metadata == "Exposure-ms")[1]+2])-1)
  message(sprintf("fluo-view exposure: %s ms", exposure))

  # Elapsed time (in ms) of each frame from the fluorescence camera relative to the flash
  # elapsedtimefl <- metadata[grep("ElapsedTime-ms", metadata)+2] # Switched to a more accurate time stamp info PVCAM-TimeStampBOF
  # elapsedtimefl <- as.numeric(substr(elapsedtimefl, 1, nchar(elapsedtimefl)-1))
  elapsedtimefl <- as.numeric(metadata[grep("PVCAM-TimeStampBOF", metadata)+2])
  elapsedtimefl <- elapsedtimefl/10
  fpsfl <- round(length(elapsedtimefl)/tail(elapsedtimefl, n=1)*1000)
  message(paste("fluo-view fps:", fpsfl))
  elapsedtimefl <- elapsedtimefl - elapsedtimefl[fluo_flash$flflashes[1]]
  elapsedtimediff <- diff(elapsedtimefl)
  png(file=paste0(output, "_elapsedtimefldiff.png"), width=400, height=400)
  plot(1:length(elapsedtimediff), elapsedtimediff)
  dev.off()
  
  # Elapsed time (in ms) of each frame from the fly view camera
  timestampusec <- as.numeric(log[grep("TimeStamp", log)+1])
  elapsedtimefv <- (timestampusec - timestampusec[1])/1000
  elapsedtimefv <- elapsedtimefv - elapsedtimefv[fly_flash$fvflashes[1]]
  fpsfv <- round(length(elapsedtimefv)/((tail(elapsedtimefv, n=1) - elapsedtimefv[1])/1000))
  message(paste("fly-view fps:", fpsfv))
  framediff <- diff(elapsedtimefv)

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
  elapsedtimeavflash <- elapsedtimeav - elapsedtimeav[arena_flash$avflashes[1]]
  fpsav <- round(length(elapsedtimeav)/((tail(elapsedtimeav, n=1) - elapsedtimeav[1])/1000))
  message(paste("arena-view fps:", fpsav))

  # Align frames between fly-view and fluo-view
  message("Aligning frames between fly-view and fluo-view")
  if(file.exists(paste0(output, "_frid.RDS"))==T & reuse==T){
    message("Loading RDS file")
    frid <- readRDS(paste0(output, "_frid.RDS"))
  }else{
    frameratio <- round(fpsfv/fpsfl)
    message(paste0("fv/fl frame ratio: ", frameratio))
    frid <- match(timestampfl, elapsedtimefv)
    frid[which(is.na(frid))] <- sapply(timestampfl[which(is.na(frid))], function(x) which.min(abs(elapsedtimefv-x)))
    
    # Check if two flashes match
    fvflashesfrid <- frid[fluo_flash$flflashes]
    message(sprintf("Hypothetical flash for fly-view: %s", paste(fvflashesfrid, collapse=" ")))
    message(sprintf("Actual flash for fly-view: %s", paste(fly_flash$fvflashes, collapse=" ")))
    if(all.equal(fvflashesfrid, fly_flash$fvflashes)==T){
      message(paste0("Hypothetical flashes match actual flashes in fly-view."))
    }else{
      message(paste0("Hypothetical flashes didn't match actual flashes in fly-view. Syncing failed?"))
    }
    if(length(fluo_flash$flflashes)!=length(fly_flash$fvflashes)) {
      warning("Number of flash in fly-view and fluo-view didn't match.")
    }
    frid <- frid[which(frid<=fly_flash$nframesfv)]
    if(length(frid) < fluo_flash$nframesfl){
      frid <- c(frid, rep(tail(frid, 1), fluo_flash$nframesfl-length(frid)))
    }
    saveRDS(frid, paste0(output, "_frid.RDS"))
  }

  # Align frames between fluo-view and arena-view
  message("Aligning frames between fluo-view and arena-view")
  if(file.exists(paste0(output, "_frida.RDS")) & reuse==T){
    message("Loading RDS file")
    frida <- readRDS(paste0(output, "_frida.RDS"))
  }else{
    frid <- match(timestampfl, elapsedtimefv)
    frid[which(is.na(frid))] <- sapply(timestampfl[which(is.na(frid))], function(x) which.min(abs(elapsedtimefv-x)))
    
    frameratio2 <- round(fpsav/fpsfl)
    message(paste0("av/fl frame ratio: ", frameratio2))
    frida <- match(timestampfl, elapsedtimeav)
    frida[which(is.na(frida))] <- sapply(timestampfl[which(is.na(frida))], function(x) which.min(abs(elapsedtimeav-x)))

        # Check if two flashes match
    avflashesfridav <- frida[fluo_flash$flflashes]
    message(sprintf("Hypothetical flash for arena-view: %s", paste(avflashesfridav, collapse=" ")))
    message(sprintf("Actual flash for arena-view: %s", paste(arena_flash$avflashes, collapse=" ")))
    if(all.equal(avflashesfridav, arena_flash$avflashes)==T){
      message(paste0("Hypothetical flashes match actual flashes in arena-view."))
    }else{
      message(paste0("Hypothetical flashes didn't match actual flashes in arena-view. Syncing failed?"))
    }
    if(length(fluo_flash$flflashes)!=length(arena_flash$avflashes)) {
      message("Number of flash in arena-view and fluo-fiew did not match.")
    }
    frida <- frida[which(frida<=arena_flash$nframesav)]
    if(length(frida) < fluo_flash$nframesfl){
      frida <- c(frida, rep(tail(frida, 1), fluo_flash$nframesfl-length(frida)))
    }
    saveRDS(frida, paste0(output, "_frida.RDS"))

  }
  return(list("fpsfv"=fpsfv, "elapsedtimeavflash"=elapsedtimeavflash, "frid"=frid, "frida"=frida))
}
