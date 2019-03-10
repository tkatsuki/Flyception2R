#' Synchronize camera frames
#'
#'
#' @param dir Input directory.
#' @param fluo_flash Result of detect_flash() function for fluo-view.
#' @param fly_flash Result of detect_flash() function for fly-view.
#' @param arena_flash Result of detect_flash() function for arena-view.
#' @param output Output path for saving results.
#' @param reuse True if one wants to reuse previously saved results.
#' @param hypothetical True if one wants to use hypothetically generated 20.0ms, 19.9ms alternating timestamp for fluo-view.
#' @export
#' @examples
#' sync_frames()

sync_frames <- function(dir, fluo_flash, fly_flash, arena_flash, output, reuse=F, hypothetical=F){

  # Start time
  message("Analyzing metadata")
  metadata <-scan(paste0(dir, list.files(dir, pattern="metadata\\.txt$")), what=character(),sep="")
  log <- scan(paste0(dir, list.files(dir, pattern="fv-log-")), what=character(),sep="")
  avlog <- scan(paste0(dir, list.files(dir, pattern="av-log-")), what=character(),sep="")

  # Exposure and binning
  exposure <- substr(metadata[which(metadata == "Exposure-ms")[1]+2], 1, nchar(metadata[which(metadata == "Exposure-ms")[1]+2])-1)
  message(sprintf("fluo-view exposure: %s ms", exposure))

  # Elapsed time (in ms) of each frame from the fluorescence camera relative to the flash
  elapsedtimefl <- as.numeric(metadata[grep("PVCAM-TimeStampBOF", metadata)+2])
  elapsedtimefl <- elapsedtimefl/10
  fpsfl <- round(length(elapsedtimefl)/tail(elapsedtimefl, n=1)*1000)
  message(paste("fluo-view fps:", fpsfl))
    # Add 5 ms of lead as the flash fires 5 ms into a frame
  elapsedtimefl <- elapsedtimefl - elapsedtimefl[fluo_flash$flflashes[1]] + 5
  elapsedtimefldiff <- diff(elapsedtimefl)
  png(file=paste0(output, "_elapsedtimefldiff.png"), width=400, height=400)
  plot(1:length(elapsedtimefldiff), elapsedtimefldiff)
  dev.off()
  
  if(hypothetical==T){
    # Generate hypothetical timestamp to ignore software jitter
    # Check if the sequence is 
    oddmean <- mean(elapsedtimefldiff[seq(from=1,to=length(elapsedtimefldiff), by=2)])
    if(oddmean > 19.99){
      elapsedtimefl <- cumsum(rep(c(20,19.9),length(elapsedtimefl)/2))
      elapsedtimefl <- elapsedtimefl - elapsedtimefl[fluo_flash$flflashes[1]] + 5
    }else{
      elapsedtimefl <- cumsum(rep(c(19.9,20),length(elapsedtimefl)/2))
      elapsedtimefl <- elapsedtimefl - elapsedtimefl[fluo_flash$flflashes[1]] + 5
    }
  }

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
    avtimestampsec[t] <- avtimestampcyclesec[t] + 128*cnt
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
    frid <- match(elapsedtimefl, elapsedtimefv)
    frid[which(is.na(frid))] <- sapply(elapsedtimefl[which(is.na(frid))], function(x) which.min(abs(elapsedtimefv-x)))
    
    # Check if two flashes match
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
    frameratio2 <- round(fpsav/fpsfl)
    message(paste0("av/fl frame ratio: ", frameratio2))
    frida <- match(elapsedtimefl, elapsedtimeavflash)
    frida[which(is.na(frida))] <- sapply(elapsedtimefl[which(is.na(frida))], function(x) which.min(abs(elapsedtimeavflash-x)))

    # Check if two flashes match
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
