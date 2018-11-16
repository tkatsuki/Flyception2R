#' Flyception2R main script
#'
#' @param dir path to the directory that contains the data
#' @param autopos logical. Perform camera alignment using FNCC?
#' @param interaction logical. Perform interaction detection? Requires two flies.
#' @param reuse logical. Reuse .RDS files?
#' @param fmf2tif logical. Convert fly-view and arena-view fmf files into tif format?
#' @param zoom numeric. Zoom factor between fly-view and fluo-view cameras.
#' @param FOI a vector of two numbers indicating the first and last frames to be analyzed. If not specified, all frames will be analyzed.
#' @param ROI a vector of four numbers indicating the first and last frames to be analyzed. If not specified, all frames will be analyzed.
#' @param binning integer. Binning of the fluo-view camera.
#' @param fluo_flash_thresh numeric. A threshold for detecting flashes in a fluo-view video.
#' @param fv_flash_thresh integer. A threshold for detecting flashes in a fly-view video.
#' @param av_flash_thresh integer. A threshold for detecting flashes in a arena-view video.
#' @param dist_thresh numeric. A distance threshold for detecting fly-fly interactions.
#' @param rotate_camera integer. Angle of the fluo-view camera.
#' @param window_size a vector of two numbers indicating the size of a window to the brain.
#' @param window_offset a vector of two numbers indicating the position of the window to the brain as an offset from the center of the image.
#' @param flash 1 if the first flash is good, 2 if the first flash is bad and the second flash is good.
#' @export
#' @examples
#' Flyception2R()
#'

Flyception2R <- function(dir, autopos=T, interaction=T, reuse=T,
                         fmf2tif=F, zoom=1.085, FOI=F, ROI=c(391, 7, 240, 240), 
                         binning=1, fluo_flash_thresh=500,
                         fv_flash_thresh=240, av_flash_thresh=100, dist_thresh=4,
                         rotate_camera=-180, window_size=c(68, 28), window_offset=c(-4, 25),
                         colorRange= c(180, 220), flash=1, preprocess=F){
  
  # TO DO
  
  ## Part 0. Initialization
  # Prepare directories and paths
  if(preprocess == T) reuse <- F
  prefix <- strsplit(dir, "/")[[1]][length(strsplit(dir, "/")[[1]])]
  outdir <- paste0(dir, paste0(FOI, collapse="_"), "/")
  dir.create(outdir)
  output_prefix <- paste0(outdir, prefix)
  
  # Start logging 
  loggit::setLogFile(paste0(dir, prefix, "_log.json"))
  
  if(preprocess == T) {
    
    loggit::message(paste0("Preprocessing", prefix, "..."))
    
    # Prepare filenames 
    fluo_view_tif <- paste0(dir, list.files(dir, pattern="Pos0\\.ome\\.tif$"))
    fluo_view_num_vids = length(list.files(dir, pattern="Pos0.*\\.ome\\.tif$"))
    fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
    arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))
    
    # Crop a first channel in fluo_view images using ImageJ
    if(length(list.files(dir, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))==0){
      imageJ_crop_append(dir, ch=1, roi=ROI) # x and y coordinates of the top left corner, width, height
    }
    if(fluo_view_num_vids < 2) {
      fluo_view_tif_ch1 <- paste0(dir, list.files(dir, pattern="ome\\.ch1\\.crop\\.tif$"))
    } else {
      fluo_view_tif_ch1 <- paste0(dir, list.files(dir, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))
    }
   
    flnframe <- dipr::readTIFF2(fluo_view_tif_ch1, getFrames = T)
    
    ## Part 1. Detect flash
    message("Detecting flash in fluo-view")
    fluo_flash <- detect_flash(input=fluo_view_tif_ch1,
                               type="fluo",
                               output=paste0(dir, prefix),
                               flash_thresh=fluo_flash_thresh,
                               reuse=reuse)
    message("Detecting flash in fly-view")
    fly_flash <- detect_flash(input=fly_view_fmf,
                              type="fly",
                              output=paste0(dir, prefix),
                              flash_thresh=fv_flash_thresh,
                              reuse=reuse)
    message("Detecting flash in arena-view")
    arena_flash <- detect_flash(input=arena_view_fmf,
                                type="arena",
                                output=paste0(dir, prefix),
                                flash_thresh=av_flash_thresh,
                                reuse=reuse)
    
    if(flash == 2){
      # fluo-view can miss flashes. If only one flash was detected leave as is.
      if(length(fluo_flash$flflashes)==2){
        fluo_flash$flflashes[1] <- fluo_flash$flflashes[2]
      }
      fly_flash$fvflashes[1] <- fly_flash$fvflashes[2]
      arena_flash$avflashes[1] <- arena_flash$avflashes[2]
    }
    if(flash != 1 && flash !=2){
      if(length(fluo_flash$flflashes) != length(fly_flash$fvflashes)){
        stop("Number of flash detected differ between fluo-view and fly-view.")
      }
    }
    
    ## Part 2. Camera alignment
    # Load fluo-view flash frames for alignment
    fl1ref <- dipr::readTIFF2(fluo_view_tif_ch1, frames=fluo_flash$flflashes[1])
    fl1ref <- EBImage::normalize(fl1ref)
    fl2ref <- dipr::readTIFF2(fluo_view_tif, frames=fluo_flash$flflashes[1])
    fl2ref <- EBImage::normalize(fl2ref)
    fl2refcrop <- fl2ref[1025:2048,1:256] # Split the original image into two halves
    
    # Align two channels of fluo-view
    center <- align_cameras(source=fl2refcrop,
                            template=fl1ref,
                            output=paste0(paste0(dir, prefix), "_fl2fl1"),
                            center=c(0, 0),
                            zoom=1,
                            autopos=T,
                            ROI=c(1, 1, 450, 50))
    
    # Crop a second channel in fluo_view images using ImageJ
    if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
      imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
    }
    if(fluo_view_num_vids < 2) {
      fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.tif$"))
    } else {
      fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
    }
    
    
    
    # Synchronize video frames
    syncing <- sync_frames(dir=dir,
                           fluo_flash=fluo_flash,
                           fly_flash=fly_flash,
                           arena_flash=arena_flash,
                           output=paste0(dir, prefix),
                           reuse=reuse,
                           hypothetical=F)
    
    # Load fly-view flash image
    fvref <- dipr::readFMF(fly_view_fmf, frames=c(fly_flash$fvflashes[1] + 1))[,,1]
    
    # Align fly-view and fluo-view
    center2 <- align_cameras(source=fvref/255,
                             template=flip(fl1ref),
                             output=paste0(paste0(dir, prefix), "_fvfl1"),
                             center=c(0, 0),
                             zoom=1.085,
                             autopos=T,
                             ROI=c(1, 1, 50, 50))
    
    # Check and align template
    ans <- c("N","Y")
    while(!all(stringr::str_to_lower(ans)=="y")){
            
      print(EBImage::display(abind(8*EBImage::resize(fvref/255, dim(fvref)[1]*zoom)[11:250, 11:250]^2,
                                   .75*(EBImage::translate(flip(fl1ref), center2)),
                                   along=3)))
      
      fvfl1ol <- EBImage::resize(fvref/255, dim(fvref)[1]*zoom)[11:250, 11:250] + .75*(EBImage::translate(flip(fl1ref), center2))      
      EBImage::writeImage(normalize(fvfl1ol), file=paste0(output_prefix, "_fvfl1_overlay.tif"))
            
      print(sprintf("Current template center is x=%d y=%d", center2[1], center2[2]))
      ans[1] <- readline("Is template match okay (Y or N)?:")
      if(!stringr::str_to_lower(ans[1])=="y") {
        center2[1] <- as.integer(readline("Enter new x position for center:"))
        center2[2] <- as.integer(readline("Enter new y position for center:"))
      }
    }
    
    savefn <- paste0(dir, prefix,"_prepdata.RData")
    save(arena_view_fmf,center,center2,flnframe,fluo_view_tif_ch1,fluo_view_tif_ch2,fly_view_fmf,savefn,syncing,file=savefn)
    loggit::message("Preprocessing done")
    return()
  } else {
    # Load preprocessed data
    load(paste0(dir, prefix,"_prepdata.RData"))
  }
  
  ## Part 3. Analyze trajectories
  trj_res <- analyze_trajectories(dir=dir,
                                  output=output_prefix,
                                  fpsfv=syncing$fpsfv,
                                  interaction=interaction)
  
  
  ## Part 4. Detect interaction
  if(interaction==T){
    message("Detecting interaction")
    closefr <- which(trj_res$flydist < dist_thresh)
    closefrid <- sapply(closefr, function(x) which.min(abs(syncing$frida-x)))
    write.table(closefrid, paste0(output_prefix, "_closefrid.txt"))
  }
  
  ## Part 5. Image registration
  message(sprintf("Reading %s", fluo_view_tif_ch1))
  
  # Analyze only part of the movie?
  if(FOI!=F && length(FOI)==2){
    flimg1 <- dipr::readTIFF2(fluo_view_tif_ch1, start=FOI[1], end=FOI[2])
    flimg2 <- dipr::readTIFF2(fluo_view_tif_ch2, start=FOI[1], end=FOI[2])
    flimg1 <- flip(flimg1) # flip images to match fly-view
    flimg2 <- flip(flimg2) # flip images to match fly-view
    message(sprintf("Fluo-view frames from %d to %d will be analyzed.", FOI[1], FOI[2]))
    frid <- syncing$frid[FOI[1]:FOI[2]]
    frida <- syncing$frida[FOI[1]:FOI[2]]
  }else{
    message("All frames will be analyzed.")
    frid <- syncing$frid
    frida <- syncing$frida
    FOI <- c(1, flnframe)
  }
  
  # Load fly-view camera images
  fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)
  
  # Apply resize and translation to align with fluo-view
  fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear"), -center2)
  # Match the size of the fvimgl and flimg by cropping (needs a better way though)
  fvimgl <- fvimgl[11:250,11:250,1:dim(fvimgl)[3]]
  EBImage::writeImage(fvimgl/255, file=paste0(output_prefix, "_fvimgl.tif"))
  
  # Detect beads
  fvimglbl <- gblur(fvimgl/255, 2)
  fvimglbw <- thresh(fvimglbl, w=20, h=20, offset=0.2)
  rm(fvimglbl)
  centermask <- drawCircle(matrix(0,dim(fvimglbw)[1],dim(fvimglbw)[2]), dim(fvimglbw)[1]/2,
                           dim(fvimglbw)[2]/2, dim(fvimglbw)[1]*2/5, col=1, fill=1)
  fvimglbw <- ssweep(fvimglbw, centermask, op="*")
  fvimglbwseg <- bwlabel(fvimglbw)
  EBImage::writeImage(fvimglbwseg, file=paste0(output_prefix, "_fvimglbwseg.tif"))
  
  # Calculate head angle
  ang_res <- detect_angle(fvimglbwseg)
  ang <- ang_res[[1]]
  centroid <- ang_res[[2]]
  markernum <- ang_res[[3]]
  
  # Find good frames
  angdiff <- c(0, diff(ang))
  ang_thresh <- 0.02
  goodangfr <- which(angdiff < ang_thresh & angdiff > -ang_thresh)
  goodmarkerfr <- which(markernum == 3)
  
  objdist <- sqrt((centroid[,1]-dim(fvimgl)[1]/2)^2 + (centroid[,2]-dim(fvimgl)[2]/2)^2)
  motion <- c(0, sqrt(diff(centroid[,1])^2 + diff(centroid[,2])^2))
  png(file=paste0(output_prefix, "_motion.png"), width=400, height=400)
  plot(motion)
  dev.off()
  motion_thresh <- 2
  goodmotionfr <- which(motion < motion_thresh)
  
  LoGkern <- round(dipr::LoG(9,9,1.4)*428.5)
  flimg2log <- EBImage::filter2(flimg2, LoGkern)
  centermask <- EBImage::drawCircle(flimg2[,,1]*0, dim(flimg2)[1]/2, dim(flimg2)[2]/2, 100, col=1, fill=T)
  flimg2cntlog <- dipr::ssweep(flimg2log, centermask, op="*")
  quantcnt <- apply(flimg2cntlog, 3, function(x) quantile(x, 0.9))
  png(file=paste0(output_prefix, "_quantcnt.png"), width=400, height=400)
  plot(quantcnt)
  dev.off()
  goodfocusfr <- which(quantcnt > 1000)
  goodfr <- Reduce(intersect, list(goodmarkerfr, goodmotionfr, goodangfr, goodfocusfr))
  
  # Save index of good frames
  if(FOI != F) {
    saveRDS(goodfr + (FOI[1] - 1), paste0(output_prefix, "_gfrid.RDS"))
  } else {
    saveRDS(goodfr, paste0(output_prefix, "_gfrid.RDS"))
  }
  
  # Load arena-view camera images
  avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
  EBImage::writeImage(avimgl/255, file=paste0(output_prefix, "_avimgl.tif"))
  EBImage::writeImage(avimgl[,,goodfr]/255, file=paste0(output_prefix, "_avimgl_goodfr.tif"))
  rm(avimgl)
  
  # Apply rotation compensation
  rot <- fvimgl[,,goodfr]
  for (r in 1:dim(rot)[3]){
    rot[,,r] <- RNiftyReg::rotate(fvimgl[,,goodfr[r]], ang[goodfr[r]], anchor = c("center"))
  }
  
  # Template matching
  centers <- array(0, dim=c(dim(rot)[3],2))
  
  for (c in 1:dim(rot)[3]){
    centers[c,] <- align_cameras(source=rot[,,c],
                                 template=rot[,,1],
                                 output=output_prefix,
                                 center=c(0, 0),
                                 zoom=1,
                                 autopos=T)
  }
  # Apply translation compensation
  rottrans <- fvimgl[,,goodfr]
  for (tr in 1:dim(rottrans)[3]){
    rottrans[,,tr] <- EBImage::translate(rot[,,tr], -centers[tr,])
  }
  EBImage::writeImage(rottrans/255, file=paste0(output_prefix, "_rottrans.tif"))
  rm(fvimgl)
  rm(rot)
  rm(fvimglbw)
  rm(fvimglbwseg)
  rm(flimg2log)
  
  ## Apply transformation functions to fluo-view images
  redrot <- flimg2[,,goodfr]
  for (rr in 1:dim(redrot)[3]){
    redrot[,,rr] <- as.Image(RNiftyReg::rotate(flimg2[,,goodfr[rr]], ang[goodfr[rr]], anchor = c("center")))
  }
  greenrot <- flimg1[,,goodfr]
  for (rg in 1:dim(greenrot)[3]){
    greenrot[,,rg] <- as.Image(RNiftyReg::rotate(flimg1[,,goodfr[rg]], ang[goodfr[rg]], anchor = c("center")))
  }
  
  redrottrans <- redrot
  for (trr in 1:dim(redrottrans)[3]){
    redrottrans[,,trr] <- EBImage::translate(redrot[,,trr], -centers[trr,])
  }
  
  greenrottrans <- greenrot
  for (trg in 1:dim(greenrottrans)[3]){
    greenrottrans[,,trg] <- EBImage::translate(greenrot[,,trg], -centers[trg,])
  }
  rm(redrot)
  rm(greenrot)
  
  ## Part 6. Image segmentation and fluorescence quantification
  # Interactively determine window size and offset to include neurons of interest
  ans <- c("N","N")
  while(!all(stringr::str_to_lower(ans)=="y")){

    redwindow <- redrottrans[(dim(redrottrans)[1]/2 + window_offset[1] - window_size[1]/2):
                               (dim(redrottrans)[1]/2 + window_offset[1] + window_size[1]/2),
                             (dim(redrottrans)[2]/2 + window_offset[2] - window_size[2]/2):
                               (dim(redrottrans)[2]/2 + window_offset[2] + window_size[2]/2),]
    greenwindow <- greenrottrans[(dim(greenrottrans)[1]/2 + window_offset[1] - window_size[1]/2):
                                   (dim(greenrottrans)[1]/2 + window_offset[1] + window_size[1]/2),
                                 (dim(greenrottrans)[2]/2 + window_offset[2] - window_size[2]/2):
                                   (dim(greenrottrans)[2]/2 + window_offset[2] + window_size[2]/2),]
    
    # Show both channels when selecting window/offset
    print(EBImage::display(abind(normalize(redwindow),normalize(greenwindow),along=2)))
    EBImage::writeImage(abind(normalize(redwindow),normalize(greenwindow),along=2), file=paste0(output_prefix, "_redwindow.tif"))
    
    print(sprintf("Current window_size is x=%d y=%d", window_size[1], window_size[2]))
    print(sprintf("Current window_offset is x=%d y=%d", window_offset[1], window_offset[2]))
    ans[1] <- readline("Check redwindow.tif. Is the window size good (Y or N)?:")
    if(ans[1] != "Y" && ans[1] != "y") {
      window_size[1] <- as.integer(readline("Enter new x size:"))
      window_size[2] <- as.integer(readline("Enter new y size:"))
    }
    ans[2] <- readline("Check redwindow.tif. Is the window offset good (Y or N)?:")
    if(ans[2] != "Y" && ans[2] != "y") {
      window_offset[1] <- as.integer(readline("Enter new x offset:"))
      window_offset[2] <- as.integer(readline("Enter new y offset:"))
    }
  }
  
  rm(flimg1)
  rm(flimg2)
  rm(flimg2cntlog)
  
  redwindowmed <- EBImage::medianFilter(redwindow/2^16, size=2)
  greenwindowmed <- EBImage::medianFilter(greenwindow/2^16, size=2)
  redwindowmedth <- EBImage::thresh(redwindowmed, w=10, h=10, offset=0.0003)
  
  # Create F_ratio images  
  redmasked <- redwindowmed*redwindowmedth
  greenmasked <- greenwindowmed*redwindowmedth
  greenperred <- greenmasked/redmasked
  redave <- colMeans(redmasked, dim=2, na.rm=T)
  greenave <- colMeans(greenmasked, dim=2, na.rm=T)
  greenperredave <- colMeans(greenperred, dim=2, na.rm=T)
  goodfrratidx <- !is.na(greenperredave)
  greenperredave <- greenperredave[goodfrratidx]
  goodfrrat <- goodfr[goodfrratidx]
  redave <- redave[goodfrratidx]
  greenave <- greenave[goodfrratidx]

  greenperred[which(is.na(greenperred)==T)] <- 0
  grratiocolor <- array(0, dim=c(dim(greenperred)[c(1,2)], 3, dim(greenperred)[3]))
  for(cfr in 1:dim(greenperred)[3]){
    grratiocolor[,,,cfr] <- dipr::pseudoColor(greenperred[,,cfr], colorRange[1], colorRange[2])
  }
  grratiocolor <- Image(grratiocolor, colormode="Color")
  
  # Overlay fly_view and F_ratio image
  rottransmask <- array(0, dim=c(dim(rottrans)[c(1,2)], dim(rottrans)[3]))
  rottransmask[(dim(rottrans)[1]/2 + window_offset[1] - window_size[1]/2):
                 (dim(rottrans)[1]/2 + window_offset[1] + window_size[1]/2),
               (dim(rottrans)[2]/2 + window_offset[2] - window_size[2]/2):
                 (dim(rottrans)[2]/2 + window_offset[2] + window_size[2]/2),] <- redwindowmedth
  
  rottranscolor <- array(0, dim=c(dim(rottrans)[c(1,2)], 3, dim(rottrans)[3]))
  rottranscolor[,,1,] <- rottrans/255*(1-rottransmask)
  rottranscolor[,,2,] <- rottrans/255*(1-rottransmask)
  rottranscolor[,,3,] <- rottrans/255*(1-rottransmask)
  
  grratiocolorl <- rottranscolor*0
  grratiocolorl[(dim(grratiocolorl)[1]/2 + window_offset[1] - window_size[1]/2):
                  (dim(grratiocolorl)[1]/2 + window_offset[1] + window_size[1]/2),
                (dim(grratiocolorl)[2]/2 + window_offset[2] - window_size[2]/2):
                  (dim(grratiocolorl)[2]/2 + window_offset[2] + window_size[2]/2),,] <- grratiocolor
  flyviewcolor <- rottranscolor + grratiocolorl
  flyviewcolor <- Image(flyviewcolor, colormode="Color")
  #display(flyviewcolor)
  EBImage::writeImage(flyviewcolor, file=paste0(output_prefix, "_flyviewcolor.tif"))
  rm(rottranscolor)
  
  # overlay red channel and F_ratio color image
  redrottranscol <- array(0, dim=c(dim(redrottrans)[c(1,2)], 3, dim(redrottrans)[3]))
  redrottranscol[,,1,] <- redrottrans*(1-rottransmask)
  redrottranscol[,,2,] <- redrottrans*(1-rottransmask)
  redrottranscol[,,3,] <- redrottrans*(1-rottransmask)
  redrottranscol <- normalize(redrottranscol, separate=F, inputRange=c(180, 400))
  redcolor <- redrottranscol + grratiocolorl
  redcolor <- Image(redcolor, colormode="Color")
  rm(redrottranscol)
  
  # Create side-by-side view of fly_view and fluo_view images
  frgcombined <- array(dim=c(dim(rottrans)[1]*4, dim(rottrans)[2], 3, dim(rottrans)[3]))
  frgcombined[1:240,1:240,1,1:dim(rottrans)[3]] <- normalize(rottrans, separate=F)
  frgcombined[1:240,1:240,2,1:dim(rottrans)[3]] <- normalize(rottrans, separate=F)
  frgcombined[1:240,1:240,3,1:dim(rottrans)[3]] <- normalize(rottrans, separate=F)
  frgcombined[241:480,1:240,1,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F, inputRange=c(180, 400))
  frgcombined[241:480,1:240,2,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F, inputRange=c(180, 400))
  frgcombined[241:480,1:240,3,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F, inputRange=c(180, 400))
  frgcombined[481:720,1:240,1,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F, inputRange=c(180, 300))
  frgcombined[481:720,1:240,2,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F, inputRange=c(180, 300))
  frgcombined[481:720,1:240,3,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F, inputRange=c(180, 300))
  frgcombined[721:960,1:240,,1:dim(redrottrans)[3]] <- redcolor
  frgcombined <-  Image(frgcombined, colormode="Color")
  
  EBImage::writeImage(normalize(redrottrans, separate=F, inputRange=c(180, 400)), file=paste0(output_prefix, "_redrottrans.tif"))
  rm(redrottrans)
  EBImage::writeImage(normalize(greenrottrans, separate=F, inputRange=c(180, 300)), file=paste0(output_prefix, "_greenrottrans.tif"))
  rm(greenrottrans)
  EBImage::writeImage(frgcombined, file=paste0(output_prefix, "_frgcombined_goodfr20_normalized.tif"))
  rm(frgcombined)
  
  # Calculate dF/F
  datrawint <- data.frame(x=goodfrrat, y=greenperredave, r=redave, g=greenave)
  
  datloessint <- loess(y ~ x, data=datrawint, span=0.4)
  redloessint <- loess(r ~ x, data=datrawint, span=0.4)
  greenloessint <- loess(g ~ x, data=datrawint, span=0.4)
  datsmoothint <- data.frame(x=goodfrrat, y=predict(datloessint))
  redsmoothint <- data.frame(x=goodfrrat, y=predict(redloessint))
  greensmoothint <- data.frame(x=goodfrrat, y=predict(greenloessint))
  intensity <- zoo::rollmean(greenperredave, 3, align="left")
  datint <- data.frame(x=goodfrrat[1:(length(goodfrrat)-2)], y=intensity)
  
  png(file=paste0(output_prefix, "_datint.png"), width=400, height=400)
  plot(datint)  
  dev.off()
  
  png(file=paste0(output_prefix, "_datsmoothint.png"), width=400, height=400)
  par(mar = c(5,5,2,5))
  plot(datrawint$x, datrawint$y, type="p", pch=16, ylab="F_ratio", xlab="frame", cex=0.5)
  lines(datsmoothint, col="blue")
  par(new=TRUE)
  plot(redsmoothint, col="red", ylim=c(0, 0.0005), axes=F, xlab=NA, ylab=NA)
  par(new=TRUE)
  plot(greensmoothint, col="green", ylim=c(0, 0.0005), axes=F, xlab=NA, ylab=NA)
  axis(side = 4)
  mtext(side = 4, line = 3, 'Fluorescence intensity')
  dev.off()
  
  saveRDS(datrawint, paste0(output_prefix, "_datrawint.RDS"))
  write.table(datrawint, paste0(output_prefix, "_datrawint.csv"), sep = ",", row.names=F)
  saveRDS(datloessint, paste0(output_prefix, "_datloessint.RDS"))
  saveRDS(datsmoothint, paste0(output_prefix, "_datloessint.RDS"))
  saveRDS(datint, paste0(output_prefix, "_datint.RDS"))
  
  F0int <- intensity[1]
  deltaFint <- intensity - F0int
  dFF0int <- deltaFint/F0int * 100
  datdFF0 <- data.frame(x=goodfrrat[1:(length(goodfrrat)-2)], y=dFF0int)
  png(file=paste0(output_prefix, "_datdFF0.png"), width=400, height=400)
  plot(datdFF0)
  dev.off()
  
  p <- ggplot2::ggplot(data=datdFF0, ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_smooth(method="loess", span = 0.4, level=0.95)
  ggplot2::ggsave(filename = paste0(output_prefix, "_dFF0int.pdf"), width = 8, height = 8)
  
  loggit::message(sprintf("window_size was x=%d y=%d", window_size[1], window_size[2]))
  loggit::message(sprintf("window_offset was x=%d y=%d", window_offset[1], window_offset[2]))
  loggit::message(sprintf("FOI was from %d to %d",  FOI[1], FOI[2])) 
  loggit::message(paste0("Max F_ratio intensity in this bout was ", max(intensity)))
  loggit::message(paste0("Max F_ratio smoothed intensity in this bout was ", max(datsmoothint$y)))
  loggit::message(paste0("Number of good frames was ", length(goodfr)))
  
  loggit::message(sprintf("||c(%d, %d) ||c(%d, %d) ||c(%d, %d) ||%d ||%.3f ||%.3f ||", 
                          FOI[1], FOI[2], window_size[1], window_size[2], window_offset[1], window_offset[2], length(goodfr), max(intensity), max(datsmoothint$y)))
  
  ## Part 7. Convert fmf to tif format
  if(fmf2tif==T){
    dipr::fmf2tif(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), skip=10)
    dipr::fmf2tif(paste0(dir, list.files(dir, pattern="^av.*fmf$")), skip=2)
  }
  
  message("Finished processing!")
  gc()
}
