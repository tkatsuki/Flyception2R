#' FlyceptionR main script
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' Flyception()
#'

Flyception <- function(rdir, dir, prefix, autopos=T, video_out=F, interaction=T, stimulus=100, stimlen=200,
                       stimplotlen=800, reuse=T, fmf2tif=T, fpsfl=100, zoom=0.85, FOI=F, binning=1){

  source(paste0(rdir, "plotter.R"))

  # Note that some environment dependent variables need to be correctly set
  # flimgrt: camera orientation
  # flydist: distance threashold for interaction
  # flflash_thresh: flash level might be different from setup to setup

  # TO DO
  # Output: frid, frida, registered images, good frames
  # Memory usage
  # Analyze specified range only
  # Create frame
  # Parallelism
  # FOI for readTiff

  ## Part 0. Initialization
  # Start logging
  rlogging::SetLogFile(base.file=paste0(prefix, "_log.txt"), folder=dir)
  message(dir)

  # Prepare filenames
  output_prefix <- paste0(dir, prefix)
  fluo_view_tif <- paste0(dir, list.files(dir, pattern="ome\\.tif$"))
  fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
  arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))

  # Set thresholds
  fluo_flash_thresh <- 0.01
  fv_flash_thresh <- 135
  av_flash_thresh <- 100
  dist_thresh <- 4
  zoom <- 1
  rotate_camera <- -90

  ## Part 1. Detect flash
  message("Detecting flash in fluo-view")
  fluo_flash <- detect_flash(input=fluo_view_tif, type="fluo", output=output_prefix,
                             flash_thresh=fluo_flash_thresh, reuse=reuse)
  message("Detecting flash in fly-view")
  fly_flash <- detect_flash(input=fly_view_fmf, type="fly", output=output_prefix,
                            flash_thresh=fv_flash_thresh, reuse=reuse)
  message("Detecting flash in arena-view")
  arena_flash <- detect_flash(input=arena_view_fmf, type="arena", output=output_prefix,
                              flash_thresh=av_flash_thresh, reuse=reuse)

  ## Part 2. Syncing frames and generate frame IDs
  syncing <- sync_frames(dir, fluo_flash=fluo_flash, fly_flash=fly_flash, arena_flash=arena_flash,
                         output=output_prefix, reuse=reuse)

  ## Part 3. Analyze trajectories
  trj_res <- analyze_trajectories(dir=dir, output=output_prefix, fpsfv=syncing$fpsfv, interaction=interaction)

  ## Part 4. Detect stimulus
  message("Detecting stimulus")
  fvtrj <- read.table(paste0(dir, list.files(dir, pattern="fv-traj-")))
  stimulus <- which(fvtrj[,10]==1)
  if(length(stimulus)==0){
    fridstim <- NA
    message(paste0("No stimulus was detected."))
  } else {
    stimfr <- sapply(stimulus, function(x) which.min(abs(syncing$frid-x)))
    message(paste0("Stimuli were given at the following frames:"))
    message(stimfr)
    dfstim <- data.frame(flview=stimfr, flyview=syncing$frid[stimfr], arenaview=syncing$frida[stimfr])
    write.table(dfstim, paste0(dir, prefix, "_fridstim.txt"))
  }

  ## Part 5. Detect interaction
  if(interaction==T){
    message("Detecting interaction")
    closefr <- which(trj_res$flydist < dist_thresh)
    closefrid <- sapply(closefr, function(x) which.min(abs(syncing$frida-x)))
    write.table(closefrid, paste0(dir, prefix, "_closefrid.txt"))
  }

  ## Part 6. Load images
  message(sprintf("Reading %s", fluo_view_tif))
  flimg <- EBImage::readImage(fluo_view_tif)
  flref <- EBImage::normalize(EBImage::rotate(EBImage::flip(flimg[,,fluo_flash$flflashes[1]]), rotate_camera))

  # Analyze only part of the movie?
  if(FOI!=F && length(FOI)==2){
    flimg <- flimg[,,FOI[1]:FOI[2]]
    message(sprintf("Fluo-view frames from %d to %d will be analyzed.", FOI[1], FOI[2]))
    frid <- syncing$frid[FOI[1]:FOI[2]]
    frida <- syncing$frida[FOI[1]:FOI[2]]
  }else{
    message("All frames will be analyzed.")
    frid <- syncing$frid
    frida <- syncing$frida
  }
  # Crop fluo-view movie for speed
  if(dim(flimg)[1] > 130){
    flimg <- flimg[round((dim(flimg)[1] - 128)/2):(round((dim(flimg)[1]/2+128/2))-1),
                   round((dim(flimg)[2] - 128)/2):(round((dim(flimg)[2]/2+128/2))-1),]
  }
  flimgrt <- EBImage::rotate(EBImage::flip(flimg), rotate_camera)
  # Load fly-view camera images
  fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)
  # Load arena-view camera images
  avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
  EBImage::writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
  rm(avimgl)

  ## Part 7. Detect window on the head
  fvimgbwbrfh <- detect_window(fvimgl=fvimgl, output=output_prefix, reuse=reuse)

  ## Part 8. Position calibration
  fvref <- dipr::readFMF(fly_view_fmf, frames=fly_flash$fvflashes[1])[,,1]/255
  center <- align_cameras(flref=flref, fvref=fvref, output=output_prefix, center=c(0, 0), zoom=0.95, autopos=F)

  ## Part 9. Image registration
  registered_images <- register_images(fvimgl=fvimgl, flimgrt=flimgrt, fvimgbwbrfh=fvimgbwbrfh,
                                       angles=trj_res$angles, zoom=zoom, center=center,
                                       output=output_prefix, reuse=reuse)

  ## Part 10. Look for good frames based on size, position, motion, and focus
  goodfr <- find_goodframes(window_mask=fvimgbwbrfh, fvimgl=fvimgl, output=output_prefix, reuse=reuse)

  ## Part 11. Calculate fluorescence intensity in the brain window
  message("Measuring fluorescence intensity...")
  if(file.exists(paste0(output_prefix, "_intensity_br.RDS"))==T & reuse==T){
    message("Using RDS file")
    intensity_br <- readRDS(paste0(output_prefix, "_intensity_br.RDS"))
  }else{
    intensity_br <- colSums(registered_images$fvimgbwbrfhregimg*registered_images$flimgreg, dims=2)/as.integer(goodfr$objsize)
    saveRDS(intensity_br, paste0(dir, prefix, "_intensity_br.RDS"))
  }
  intensity <- zoo::na.approx(intensity_br)

  ## Part 12. Plot delta F over F0
  message("Creating dF/F0 plots")
  F0int <- mean(intensity[1:5])
  deltaFint <- intensity - F0int
  dFF0int <- deltaFint/F0int * 100
  dat <- data.frame(x=(1:length(dFF0int)), y=dFF0int, d=trj_res$flydist[frida])
  p <- ggplot2::ggplot(data=dat, ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_smooth(method="loess", span = 0.4, level=0.95) +
    ggplot2::ylim(-5, 10) +
    ggplot2::geom_line(data=dat, ggplot2::aes(x=x, y=d))
  ggplot2::ggsave(filename = paste0(output_prefix, "_dFF0int.pdf"), width = 8, height = 8)

  ## Part 13. Create delta F over F0 pseudocolor representation
  message("Calculating dF/F0 images...")
  dF_F0_image(flimgreg, fvimgbwbrfhregimg, regimgi, colmax, cmin, cmax)

  # Create trajectory of the flies at the time of stimulus delivery
  message("Creating trajectory of the flies...")
  for(tj in 1:length(fridstim)){
    Ffr <- fridstim[tj]:(fridstim[tj]+stimplotlen-1)
    pdf(file= paste0(dir, prefix, "_trackResult_", tj, ".pdf"), width = 4.4, height = 4, bg = "white")
    par(plt = c(0, 1, 0, 1), xaxs = "i", yaxs = "i")
    plot(trja[frida[Ffr[1]:(Ffr[1]+stimlen-1)],1]*10, -trja[frida[Ffr[1]:(Ffr[1]+stimlen-1)],2]*10,
         type = "l", lty = 1, pch=tj, col="red",
         axes = F, xlim = c(-240, 240), ylim = c(-220, 220))
    par(new=T)
    plot(trja[frida[(Ffr[1]+stimlen):tail(Ffr, n=1)],1]*10, -trja[frida[(Ffr[1]+stimlen):tail(Ffr, n=1)],2]*10,
         type = "l", lty = 2, pch=tj, col="red",
         axes = F, xlim = c(-240, 240), ylim = c(-220, 220))

    if(interacting==T){
      par(new=T)
      plot(trja[frida[Ffr[1]:(Ffr[1]+stimlen-1)],3]*10, -trja[frida[Ffr[1]:(Ffr[1]+stimlen-1)],4]*10,
           type = "l", lty = 1, pch=tj, col="blue",
           axes = F, xlim = c(-240, 240), ylim = c(-220, 220))
      par(new=T)
      plot(trja[frida[(Ffr[1]+stimlen):tail(Ffr, n=1)],3]*10, -trja[frida[(Ffr[1]+stimlen):tail(Ffr, n=1)],4]*10,
           type = "l", lty = 2, pch=tj, col="blue",
           axes = F, xlim = c(-240, 240), ylim = c(-220, 220))
    }
    par(new=T)
    draw.ellipse(0,0,11.0795*20,10*20)
    dev.off()

  }

  rm(fvimgl)
  rm(flimgrt)
  rm(regimgi)

  # Plot intensity, speed, window size, etc. for both all frames and good frames
  message("Plotting results...")
  plotter(dir, prefix, intensity*16384, speed, error, objsize, quantcnt,
          elapsedtime, elapsedtimefvflash, elapsedtimeavflash, exposure, binning, fpsfl,
          starttimefl, type="all", stim = fridstim, stimlen = stimlen, flydist=flydist)
  intensity_g <- intensity
  intensity_g[-goodfr] <- NA
  speed_g <- speed
  speed_g[-frid[goodfr]] <- NA
  error_g <- error
  error_g[-frid[goodfr]] <- NA
  objsize_g <- objsize
  objsize_g[-goodfr] <- NA
  quant_g <- quantcnt
  quant_g[-goodfr] <- NA
  flydist_g <- flydist
  flydist_g[-goodfr] <- NA
  elapsedtime_g <- elapsedtime
  elapsedtime_g[-goodfr] <- NA
  elapsedtimefvflash_g <- elapsedtimefvflash
  elapsedtimefvflash_g[frid[-goodfr]] <- NA
  elapsedtimeavflash_g <- elapsedtimeavflash
  elapsedtimeavflash_g[frida[-goodfr]] <- NA
  plotter(dir, prefix, intensity_g, speed_g, error_g, objsize_g, quant_g,
          elapsedtime_g, elapsedtimefvflash_g, elapsedtimeavflash_g, exposure, binning, fpsfl,
          starttimefl, type="goodfr", stim = fridstim, stimlen = stimlen, flydist=flydist_g)

  rm(dFF0finmaskfly)

  # Convert fmf to tif format
  if(fmf2tif==T){
    fmf2tif(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), skip=10)
    fmf2tif(paste0(dir, list.files(dir, pattern="^av.*fmf$")), skip=2)
  }

  message("Finished processing!")
  gc()
}
