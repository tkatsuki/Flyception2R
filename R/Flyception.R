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
  fluo_flash <- detect_flash(input=fluo_view_tif, type="fluo", output=output_prefix, flash_thresh=fluo_flash_thresh, reuse=F)
  message("Detecting flash in fly-view")
  fly_flash <- detect_flash(input=fly_view_fmf, type="fly", output=output_prefix, flash_thresh=fv_flash_thresh, reuse=F)
  message("Detecting flash in arena-view")
  arena_flash <- detect_flash(input=arena_view_fmf, type="arena", output=output_prefix, flash_thresh=av_flash_thresh, reuse=F)

  ## Part 2. Syncing frames and generate frame IDs
  syncing <- sync_frames(fluo_flash=fluo_flash, fly_flash=fly_flash, arena_flash=arena_flash, output=output_prefix, reuse=F)

  ## Part 3. Analyze trajectories
  trj_res <- analyze_trajectories(input=dir, output=output_prefix, fpsfv=syncing$fpsfv, interaction=interaction)

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
  }
  dfstim <- data.frame(flview=stimfr, flyview=syncing$frid[stimfr], arenaview=syncing$frida[stimfr])
  write.table(dfstim, paste0(dir, prefix, "_fridstim.txt"))

  ## Part 5. Detect interaction
  if(interaction==T){
    message("Detecting interaction")
    closefr <- which(trj_res$flydist < dist_thresh)
    closefrid <- sapply(closefr, function(x) which.min(abs(syncing$frida-x)))
    write.table(closefrid, paste0(dir, prefix, "_closefrid.txt"))
  }

  ## Part 6. Load images
  message(sprintf("Reading %s", fluo_view_tif))
  flimg <- readImage(fluo_view_tif)
  flref <- normalize(rotate(flip(flimg[,,fluo_flash$flflashes[1]]), rotate_camera))

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
  flimgrt <- rotate(flip(flimg), rotate_camera)
  # Load fly-view camera images
  fvimgl <- readFMF(fly_view_fmf, frames=frid)
  # Load arena-view camera images
  avimgl <- readFMF(arena_view_fmf, frames=frida)
  writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
  rm(avimgl)

  ## Part 7. Detect window on the head
  fvimgbwbrfh <- detect_window(fvimgl=fvimgl, output=output_prefix, reuse=reuse)

  ## Part 8. Position calibration
  fvref <- readFMF(fly_view_fmf, frames=fly_flash$fvflashes[1])[,,1]/255
  center <- align_cameras(flref=flref, fvref=fvref, output=output_prefix, center=c(0, 0), zoom=0.95, autopos=F)

  ## Part 9. Image registration
  registered_images <- register_images(fvimgl=fvimgl, flimgrt=flimgrt, fvimgbwbrfh=fvimgbwbrfh, angles=trj_res$angles, zoom=zoom, center=center, output=output_prefix, reuse=F)

  ## Part 10. Filtering suggest good frames based on size and position
  goodfr <- find_goodframes(fvimgbwbrfh, output=output_prefix, reuse=reuse)

  ## Part 11. Calculate fluorescence intensity in the brain window
  message("Measuring fluorescence intensity...")
  if(file.exists(paste0(output_prefix, "_intensity_br.RDS"))==T & reuse==T){
    message("Using RDS file")
    intensity_br <- readRDS(paste0(output_prefix, "_intensity_br.RDS"))
  }else{
    intensity_br <- colSums(registered_images$fvimgbwbrfhregimg*registered_images$flimgreg, dims=2)/as.integer(objsize)
    saveRDS(intensity_br, paste0(dir, prefix, "_intensity_br.RDS"))
  }
  intensity <- na.approx(intensity_br)

  # Plot delta F over F0
  message("Creating dF/F0 plots")
  for(fintfr in 1:length(fridstim)){
    Fintfr <- fridstim[fintfr]:(fridstim[fintfr]+stimplotlen-1)
    Fintfr <- Fintfr[which(Fintfr < nframesfl)]
    Fintfrg <- intersect(goodfr, Fintfr)
    Fintfr[!Fintfr%in%goodfr] <- NA
    F0int <- mean(intensity[Fintfrg[1:5]])
    deltaFint <- intensity[Fintfr] - F0int
    dFF0int <- deltaFint/F0int * 100
    dat <- data.frame(x=(1:length(dFF0int)), y=dFF0int, d=flydist[Fintfr])
    p <- ggplot(data=dat, aes(x=x, y=y)) +
      geom_smooth(method="loess", span = 0.4, level=0.95) +
      ylim(-5, 10) +
      geom_line(data=dat, aes(x=x, y=d))
    ggsave(filename = paste0(dir, prefix, "_dFF0int_", fintfr, ".pdf"), width = 8, height = 8)
  }

  # Create delta F over F0 pseudocolor representation only for good frames
  message("Calculating dF/F0 images...")
  for(ffr in 1:length(fridstim)){
    Ffr <- fridstim[ffr]:(fridstim[ffr]+stimplotlen-1)
    Ffr <- Ffr[which(Ffr < nframesfl)]
    oFfr <- Ffr
    Ffr <- intersect(goodfr, Ffr)
    wFfr <- which(oFfr%in%Ffr)
    dFfr <- data.frame(ori=wFfr, Ffr=1:length(Ffr), fl=Ffr, fv=frid[Ffr], av=frida[Ffr])
    write.table(dFfr, file=paste0(dir, prefix, "_dFfr_", ffr, ".txt"), row.names=F)

    Fmean <- rollmeanimg(flimgreg[,,Ffr], 5)
    F0 <- rowMeans(flimgreg[,,Ffr[1:5]], dims=2)
    deltaF <- ssweep(Fmean, F0, op="-")
    dFF0 <- ssweep(deltaF, 1/F0, op="*")
    dFF0[is.na(dFF0)] <- 0
    dFF0[is.infinite(dFF0)] <- 0
    dFF0masked <- fvimgbwbrfhregimg[,,Ffr]*dFF0
    dFF0maskedpos <- dFF0masked * 100 # Convert to %
    dFF0maskedpos[which(dFF0maskedpos < 0)] <- 0

    colmax <- median(apply(dFF0maskedpos, 3, max))
    #colmax <- 300
    dFF0maskedpos <- medianFilter(dFF0maskedpos/colmax, 3) # medianFilter cuts > 1
    dFF0fin <- array(0, dim=c(dim(fvimgbwbrfhregimg)[c(1,2)], 3, length(Ffr)))
    for(cfr in 1:length(Ffr)){
      dFF0fin[,,,cfr] <- pseudoColor(dFF0maskedpos[,,cfr], 64, 256)
    }
    dFF0fin <- Image(dFF0fin, colormode="Color")
    message(sprintf("Pseudocolor range for ffr=%d is 20 to %d", ffr, colmax))

    # Use window size for filtering
    F0size <- mean(objsize[Ffr[1]:(Ffr[1]+4)])
    dFF0size <- which(objsize[Ffr] > (F0size - 50) & objsize[Ffr] < (F0size + 50))

    # Use focus for filtering
    F0focus <- mean(quantcnt[Ffr[1]:(Ffr[1]+4)])
    dFF0focus <- which(quantcnt[Ffr] > (F0focus - 20) & quantcnt[Ffr] < (F0focus + 20))

    # Use both window size and focus for filtering
    dFF0size_focus <- intersect(dFF0size, dFF0focus)
    write.table(dFF0size_focus, file=paste0(dir, prefix, "_dFF0size_focus_", ffr, ".txt"))
    writeImage((100*Fmean)^2, file=paste0(dir, prefix, "_Fmean_", ffr, ".tif"))
    dFF0finmask <- array(0, dim=c(dim(fvimgbwbrfhregimg)[c(1,2)], 3, length(Ffr)))
    dFF0finmask[,,1,] <- fvimgbwbrfhregimg[,,Ffr]
    dFF0finmask[,,2,] <- fvimgbwbrfhregimg[,,Ffr]
    dFF0finmask[,,3,] <- fvimgbwbrfhregimg[,,Ffr]
    dFF0regimg <- array(0, dim=c(dim(fvimgbwbrfhregimg)[c(1,2)], 3, length(Ffr)))
    dFF0regimg[,,1,] <- 255-regimgi[,,Ffr]
    dFF0regimg[,,2,] <- 255-regimgi[,,Ffr]
    dFF0regimg[,,3,] <- 255-regimgi[,,Ffr]
    dFF0finmaskfly <- Image(dFF0fin*dFF0finmask+dFF0regimg/255, colormode="Color")
    dFF0finmask <- Image(dFF0fin*dFF0finmask, colormode="Color")
    writeImage(dFF0finmaskfly[,,,dFF0size_focus], bits.per.sample = 8,
               file=paste0(dir, prefix, "_dFF0finmaskfly_sizefocus_", ffr, ".tif"))
    writeImage(dFF0finmaskfly, bits.per.sample = 8,
               file=paste0(dir, prefix, "_dFF0finmaskfly_", ffr, ".tif"))
    writeImage(dFF0finmask, bits.per.sample = 8,
               file=paste0(dir, prefix, "_dFF0finmask_", ffr, ".tif"))

    # Uncomment for measuring dF/F in ROI
    #      measureROI <- function(img, x, y, w, h){
    #        ROI <- img[x:(x+w-1),y:(y+h-1),]
    #        meanint <- colMeans(ROI, dim=2)
    #        return(meanint)
    #      }
    #      leftROI <- measureROI(Fmean, 80, 96, 10, 10)
    #      rightROI <- measureROI(Fmean, 110, 97, 10, 10)
    #      ROI <- (leftROI + rightROI)/2
    #
    #      ROIF0 <- mean(ROI[1:5])
    #      deltaROIF <- ROI - ROIF0
    #      ROIdFF0 <- deltaROIF/ROIF0 * 100
    #      dat <- data.frame(x=wFfr, y=ROIdFF0, d=flydist[frida[Ffr]])
    #      p <- ggplot(data=dat, aes(x=x, y=y)) +
    #        geom_smooth(method="loess", span = 0.4, level=0.95) +
    #        ylim(-5, 10) +
    #        geom_line(data=dat, aes(x=x, y=d))
    #      ggsave(filename = paste0(dir, prefix, "_dFF0int_ROI_", ffr, ".pdf"), width = 8, height = 8)

    # Save videos of stimulus frames of flyview and arenaview
    writeImage((fvimgl[,,Ffr])/255, file=paste0(dir, prefix, "_fvimgl_stim_", ffr, ".tif"))
    avimglstim <- readFMF(paste0(dir, list.files(dir, pattern="^av.*fmf$")), frames=frida[Ffr])
    writeImage(avimglstim/255, file=paste0(dir, prefix, "_avimglstim_", ffr, ".tif"))

  }

  rm(dFF0)
  rm(dFF0masked)
  rm(dFF0maskedpos)
  rm(dFF0fin)
  rm(dFF0regimg)
  rm(dFF0finmask)
  rm(fvimgbwbrfhregimg)

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
