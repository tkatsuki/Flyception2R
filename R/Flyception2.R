#' Flyception2R main script
#'
#' @param dir path to the directory that contains the data
#' @param outdir path to root of directory for analysis outputs. Directory structure of data will be created under root. Default to data directory. 
#' @param autopos logical. Perform camera alignment using FNCC?
#' @param interaction logical. Perform interaction detection? Requires two flies.
#' @param stimulus logical, whether or not stilumation (e.g., odor) is present, only works with data Flyception2 run in stimulation mode
#' @param reuse logical. Reuse .RDS files?
#' @param fmf2tif logical. Convert fly-view and arena-view fmf files into tif format?
#' @param zoom numeric. Zoom factor between fly-view and fluo-view cameras.
#' @param FOI a vector of two integers indicating the first and last frames to be analyzed. If not specified, all frames will be analyzed.
#' @param ROI a vector of four integers indicating the top left corrdinate and width and height of a ROI for cropping the green channel in the fluo-view.
#' @param binning integer. Binning of the fluo-view camera.
#' @param fluo_flash_thresh numeric. A threshold for detecting flashes in a fluo-view video.
#' @param fv_flash_thresh integer. A threshold for detecting flashes in a fly-view video.
#' @param av_flash_thresh integer. A threshold for detecting flashes in a arena-view video.
#' @param dist_thresh numeric. A distance threshold for detecting fly-fly interactions.
#' @param fl1fl2center a vector of two integers indicating the x and y offset between the two fluo-view videos. If not specified, an interactive dialog will pop up.
#' @param flvfl1center a vector of two integers indicating the x and y offset between the fluo-view and fly-view videos. If not specified, an interactive dialog will pop up.
#' @param bgratio float. The ratio of background pixels to total pixels in the ROI. Pixels with intensity values below this percentile of pixels in the ROI are excluded from segmentation mask.
#' @param ratiocutoff float. A percentile below witch ratios are removed from analysis if active region is a subset of the labeled region.
#' @param rotate_camera integer. Angle of the fluo-view camera.
#' @param window_size a list of vectors of two integers indicating the size of a window to the brain. If not specified, an interactive dialog will pop up.
#' @param window_offset a list of vectors of two integers indicating the position of the window to the brain as an offset from the center of the image. If not specified, an interactive dialog will pop up.
#' @param colorRange a vector of two integers indicating the min (blue) and max (red) for pseudocolor representation.
#' @param flash 1 if the first flash is good, 2 if the first flash is bad and the second flash is good.
#' @param preprocess logical. True if flash detection and camera alignment need to be performed.
#' @param size_thresh integer. A threshold for removing object at the segmentation step.
#' @param focus_thresh integer. A threshold for determining out-of-focus frames.
#' @param badfr vector. A list of frames which to be manually excluded from analysis.
#' @param ctr_offset vector. x,y offset from center of markers to center of brain window to compensate for coverslip/bead placement variation.
#' @param baseline 1 or 2 vector frame number to use for df/f f0. If 2 vector f0 is the mean between frames relative to FOI
#' @param input_range_r a vector of two integers setting the contrast range of red channel
#' @param input_range_g a vector of two integers setting the contrast range of green channel
#' @param motion_thresh integer, A threshold for removing frames with motion blur
#' @param stim_pattern a vector of 3 numbers, indicating pre-stimulus, stimulus, and post-stimulus duration
#' @param gen_av_trj_vid bool indicate whether to generate video with arena view tracking indicators
#' @param fly_id zero based index number of the flyview tracked fly corresponding to the arenaview trajectory columns
#' @param process_all logical. True if both preprocess and main process to be executed
#' @export
#' @examples
#' Flyception2R()
#'

Flyception2R <- function(dir, outdir=NA, autopos=T, interaction=T, stimulus=F, reuse=T, fmf2tif=F,
                         zoom=1.085, FOI=F, ROI=c(391, 7, 240, 240), binning=1, 
                         fluo_flash_thresh=500, fv_flash_thresh=240, av_flash_thresh=100, dist_thresh=4,
                         fl1fl2center=NA, flvfl1center=NA,
                         bgratio=0.80, ratiocutoff=0.00,  
                         rotate_camera=-180, window_size=NA, window_offset=NA,
                         colorRange= c(0, 200), flash=NA, preprocess=F,
                         baseline=NA, input_range_r=c(180, 400), input_range_g=c(180, 300),
                         size_thresh=5, focus_thresh=950, badfr=NA, ctr_offset=NA, motion_thresh=10,
                         stim_pattern=c(1,2,10), gen_av_trj_vid=F, fly1_id=0, process_all=F, 
                         loess_span=0.4){
  
  # TO DO
  # - why require restart
  # - fly-fly angle detection needs to be corrected
  # - redwindow contrast too low
  # - Optimize speed
  
  ## Part 0. Initialization
  # Preprocessing ----
  # Prepare directories and paths
  if(preprocess == T) reuse <- F
  prefix <- strsplit(dir, "/")[[1]][length(strsplit(dir, "/")[[1]])]
  # Set output directory if prefix specified
  if(is.na(outdir)) {
    outdirr <- dir
  } else {
    outdirr <- paste0(outdir,
                      substr(dir,nchar(strsplit(dir, "/")[[1]][1]) + 2,nchar(dir)))
  }
  
  # Create ouput directory if it doesn't exist
  dir.create(outdirr, showWarnings=FALSE, recursive=TRUE)
  
  # Start logging 
  loggit::setLogFile(paste0(outdirr, prefix, "_log.json"))
  
  if(preprocess == T | process_all == T | !file.exists(paste0(outdirr, prefix,"_prepdata.RData"))) {
    
    # Log function arguments
    args<-as.list(environment())
    loggit::message(paste(names(args),"->",args,collapse=","))
    
    
    loggit::message(paste0("Preprocessing", prefix, "..."))
    
    # Prepare filenames 
    fluo_view_tif <- paste0(dir, list.files(dir, pattern="Pos0\\.ome\\.tif$"))
    fluo_view_num_vids = length(list.files(dir, pattern="Pos0.*\\.ome\\.tif$"))
    fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
    arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))
    
    # Crop a first channel in fluo_view images using ImageJ
    if(length(list.files(outdirr, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))==0 || preprocess){
      imageJ_crop_append(dir, outdirr, ch=1, roi=ROI) # x and y coordinates of the top left corner, width, height
    }
    if(fluo_view_num_vids < 2) {
      fluo_view_tif_ch1 <- paste0(outdirr, list.files(outdirr, pattern="ome\\.ch1\\.crop\\.tif$"))
    } else {
      fluo_view_tif_ch1 <- paste0(outdirr, list.files(outdirr, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))
    }
    
    flnframe <- dipr::readTIFF2(fluo_view_tif_ch1, getFrames = T)
    
    ## Part 1. Detect flash
    loggit::message("Detecting flash in fluo-view")
    fluo_flash <- detect_flash(input=fluo_view_tif_ch1,
                               type="fluo",
                               output=paste0(outdirr, prefix),
                               flash_thresh=fluo_flash_thresh,
                               reuse=reuse)
    loggit::message("Detecting flash in fly-view")
    fly_flash <- detect_flash(input=fly_view_fmf,
                              type="fly",
                              output=paste0(outdirr, prefix),
                              flash_thresh=fv_flash_thresh,
                              reuse=reuse)
    loggit::message("Detecting flash in arena-view")
    arena_flash <- detect_flash(input=arena_view_fmf,
                                type="arena",
                                output=paste0(outdirr, prefix),
                                flash_thresh=av_flash_thresh,
                                reuse=reuse)
    
    # Assume flyview doesn't miss flashes
    num_flashes <- length(fly_flash$fvflashes)
    # Fluoview didn't miss flash 
    if(length(fluo_flash$flflashes) == num_flashes) {
      # TODO: pose/lighting during flash can be better in a particular flash frame choose best
      if(num_flashes > 1) {
        loggit::message("See display to choose flash frame...")
      } else {
        loggit::message("Only detected one flash...")
      }
      # For now use the first flash as default if Fluoview didn't miss flash
      if(is.na(flash)) flash <- 1
      loggit::message("Using first flash...")
      # Fluoview missed a flash. If second flash is specified, use it.
    } else if(flash >= length(fluo_flash$flflashes)) {
      fly_flash$fvflashes[1] <- fly_flash$fvflashes[flash]
      arena_flash$avflashes[1] <- arena_flash$avflashes[flash]
      loggit::message("Using second flash...")
      # If first flash is specified, use it, otherwise stop.
    } else if(flash != 1) {
      stop("Number of flash detected differ between fluo-view and fly-view.")
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
                            output=paste0(paste0(outdirr, prefix), "_fl2fl1"),
                            center=c(0, 0),
                            zoom=1,
                            autopos=T,
                            ROI=c(1, 1, 450, 50))
    
    # Check and align template between Fluoview cameras
    if(is.na(fl1fl2center[1])) {
      ans <- c("N","Y")
    } else {
      ans <- c("Y","Y")
      center <- fl1fl2center
      
      # Detect out of bounds
      x1 <- ROI[1] + center[1]
      y1 <- ROI[2] + center[2]
      x2 <- ROI[1] + ROI[3] + center[1]
      y2 <- ROI[2] + ROI[4] + center[2]
      
      if(x1 < 0)
        center[1] <- center[1] + x1    
      if(x2 > dim(fl2refcrop)[1])
        center[1] <- center[1] - (x2 - dim(fl2refcrop)[1])
      
      if(y1 < 0)
        center[2] <- center[2] - y1  
      if(y2 > dim(fl2refcrop)[2])
        center[2] <- center[2] - (y2 - dim(fl2refcrop)[2])
      
      x <- (ROI[1] + 1) + center[1]
      y <- (ROI[2] + 1) + center[2]
      
      fl2fl1stack <- abind(1.00*fl2refcrop[x:(x+ROI[3]-1),y:(y+ROI[4]-1)],
                           1.00*(EBImage::translate(fl1ref, center)),
                           along=3)
      
      EBImage::writeImage(normalize(fl2fl1stack),
                          file=paste0(outdirr, prefix,
                                      "_fl2fl1_stack.tif"))
    }
    while(!all(stringr::str_to_lower(ans)=="y")){
      
      # Detect out of bounds
      x1 <- ROI[1] + center[1]
      y1 <- ROI[2] + center[2]
      x2 <- ROI[1] + ROI[3] + center[1]
      y2 <- ROI[2] + ROI[4] + center[2]
      
      if(x1 < 0)
        center[1] <- center[1] + x1    
      if(x2 > dim(fl2refcrop)[1])
        center[1] <- center[1] - (x2 - dim(fl2refcrop)[1])
      
      if(y1 < 0)
        center[2] <- center[2] - y1  
      if(y2 > dim(fl2refcrop)[2])
        center[2] <- center[2] - (y2 - dim(fl2refcrop)[2])
      
      x <- (ROI[1] + 1) + center[1]
      y <- (ROI[2] + 1) + center[2]
      
      fl2fl1stack <- abind(1.00*fl2refcrop[x:(x+ROI[3]-1),y:(y+ROI[4]-1)],
                           1.00*(EBImage::translate(fl1ref, center)),
                           along=3)
      
      print(EBImage::display(fl2fl1stack))
      EBImage::writeImage(normalize(fl2fl1stack),
                          file=paste0(outdirr, prefix,
                                      "_fl2fl1_stack.tif"))
      
      print(sprintf("Current template center is x=%d y=%d", center[1], center[2]))
      ans[1] <- readline("Is template match okay (Y or N)?:")
      if(!stringr::str_to_lower(ans[1])=="y") {
        center[1] <- as.integer(readline("Enter new x position for center:"))
        center[2] <- as.integer(readline("Enter new y position for center:"))
      }
    }
    
    # Crop a second channel in fluo_view images using ImageJ
    if(length(list.files(outdirr, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0 || preprocess){
      imageJ_crop_append(dir, outdirr, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), ROI[3], ROI[4])) # x and y coordinates of the top left corner, width, height
    }
    if(fluo_view_num_vids < 2) {
      fluo_view_tif_ch2 <- paste0(outdirr, list.files(outdirr, pattern="ome\\.ch2\\.crop\\.tif$"))
    } else {
      fluo_view_tif_ch2 <- paste0(outdirr, list.files(outdirr, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
    }
    
    
    
    # Synchronize video frames
    syncing <- sync_frames(dir=dir,
                           fluo_flash=fluo_flash,
                           fly_flash=fly_flash,
                           arena_flash=arena_flash,
                           output=paste0(outdirr, prefix),
                           reuse=reuse,
                           hypothetical=F)
    
    # Load fly-view flash image
    fvref <- dipr::readFMF(fly_view_fmf, frames=c(fly_flash$fvflashes[1] + 1))[,,1]
    
    # Align fly-view and fluo-view
    center2 <- align_cameras(source=fvref/255,
                             template=flip(fl1ref),
                             output=paste0(paste0(outdirr, prefix), "_fvfl1"),
                             center=c(0, 0),
                             zoom=1.085,
                             autopos=T,
                             ROI=c(1, 1, 50, 50))
    
    # Check and align template between Flyview and Fluoview
    if(is.na(flvfl1center[1])) {
      ans <- c("N","Y")
    } else {
      ans <- c("Y","Y")
      center2 <- flvfl1center
      fvzoomcrop <- EBImage::resize(fvref/255, dim(fvref)[1]*zoom)
      x1         <- (dim(fvzoomcrop)[1] - dim(fl1ref)[1])/2
      x2         <- x1 + dim(fl1ref)[1] - 1
      fvzoomcrop <- fvzoomcrop[x1:x2,x1:x2]
      fvfl1stack <- abind(8*fvzoomcrop^2,
                          .75*(EBImage::translate(flip(fl1ref), center2)),
                          along=3)
      EBImage::writeImage(fvfl1stack,
                          file=paste0(outdirr, prefix,
                                      "_fvfl1_stack.tif"))
    }
    
    while(!all(stringr::str_to_lower(ans)=="y")){
      fvzoomcrop <- EBImage::resize(fvref/255, dim(fvref)[1]*zoom)
      x1         <- (dim(fvzoomcrop)[1] - dim(fl1ref)[1])/2
      x2         <- x1 + dim(fl1ref)[1] - 1
      fvzoomcrop <- fvzoomcrop[x1:x2,x1:x2]
      fvfl1stack <- abind(8*fvzoomcrop^2,
                          .75*(EBImage::translate(flip(fl1ref), center2)),
                          along=3)
      
      print(EBImage::display(fvfl1stack))
      EBImage::writeImage(fvfl1stack,
                          file=paste0(outdirr, prefix,
                                      "_fvfl1_stack.tif"))
      
      print(sprintf("Current template center is x=%d y=%d", center2[1], center2[2]))
      ans[1] <- readline("Is template match okay (Y or N)?:")
      if(!stringr::str_to_lower(ans[1])=="y") {
        center2[1] <- as.integer(readline("Enter new x position for center:"))
        center2[2] <- as.integer(readline("Enter new y position for center:"))
      }
    }
    
    closeAllConnections()
    savefn <- paste0(outdirr, prefix,"_prepdata.RData")
    save(arena_view_fmf,center,center2,flnframe,fluo_view_tif_ch1,fluo_view_tif_ch2,
         fly_view_fmf,savefn,syncing,fluo_flash,arena_flash,fly_flash,file=savefn)
    loggit::message("Preprocessing done")
    if(preprocess == T) return()
  } else {
    # Load preprocessed data
    loggit::message(paste0("Loading preprocessed data"))
    load(paste0(outdirr, prefix,"_prepdata.RData"))
  }
  
  ## Detect stimulus frames
  # Find stimulus frame ----
  stimfr <- NA
  if(stimulus==T){
    loggit::message("Detecting stimulus")
    fvtrj <- read.table(paste0(dir, list.files(dir, pattern="fv-traj-")))
    stimtrjfr <- which(fvtrj[,8]==1)
    if(length(stimtrjfr)==0){
      loggit::message(paste0("No stimulus was detected."))
    } else {
      stimfr <- sapply(stimtrjfr, function(x) which.min(abs(syncing$frid-x)))
      loggit::message(paste0("Stimuli were given at the following frames:"))
      loggit::message(stimfr)
      dfstim <- data.frame(flview=stimfr, flyview=syncing$frid[stimfr], arenaview=syncing$frida[stimfr])
      write.table(dfstim, paste0(dir, prefix, "_fridstim.txt"))
    }
  }else{
    loggit::message("Stimulus detection skipped")
  }
  
  if(stimulus==T) event_pattern <- c(rep(0.1, syncing$fpsfl*(stim_pattern[1])), rep(1, syncing$fpsfl*(stim_pattern[2])), rep(1.5, syncing$fpsfl*(stim_pattern[3]) + 1)) # FOI needs to be fixed
  
  
  # Set input paths relative to directory parameters
  arena_view_fmf      <- paste0(dir,tail(strsplit(arena_view_fmf,"/")[[1]],n=1))
  fly_view_fmf        <- paste0(dir,tail(strsplit(fly_view_fmf,"/")[[1]],n=1))
  fluo_view_tif_ch1   <- paste0(outdirr,tail(strsplit(fluo_view_tif_ch1,"/")[[1]],n=1))
  fluo_view_tif_ch2   <- paste0(outdirr,tail(strsplit(fluo_view_tif_ch2,"/")[[1]],n=1))
  
  # Analyze only part of the movie?
  # FOI creation ----
  # If stimulus is given override the FOI
  if(!is.na(stimfr)){
    # Could add multiple stim case
    # 1 sec before stimulus, 2 sec of stimulus, 10 sec after stimulus
    # Could make this interactive or parametric
    
    FOI <- c(stimfr[1] - syncing$fpsfl*stim_pattern[1], stimfr[1] + syncing$fpsfl*(stim_pattern[2] + stim_pattern[3]))
    if(FOI[1]<1) FOI[1] <- 1
    if(FOI[2]>flnframe) FOI[2] <- flnframe
    frid <- syncing$frid[FOI[1]:FOI[2]]
    frida <- syncing$frida[FOI[1]:FOI[2]]
  }else{
    if(FOI[1]!=F && length(FOI)==2){
      frid <- syncing$frid[FOI[1]:FOI[2]]
      frida <- syncing$frida[FOI[1]:FOI[2]]
    }else{
      frid <- syncing$frid
      frida <- syncing$frida
      FOI <- c(1, flnframe)
    }
  }
  loggit::message(sprintf("Fluo-view frames from %d to %d will be analyzed.", FOI[1], FOI[2]))
  
  outdir <- paste0(outdirr, paste0(FOI, collapse="_"), "/")
  dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
  
  output_prefix <- paste0(outdir, prefix)
  if(nchar(output_prefix)>240) stop("Directory name too long")
  
  ## Part 5. 
  # Image registration ----
  loggit::message(sprintf("Reading %s", fluo_view_tif_ch1))
  
  # Load fluo-view camera images
  flimg1 <- dipr::readTIFF2(fluo_view_tif_ch1, start=FOI[1], end=FOI[2])
  flimg2 <- dipr::readTIFF2(fluo_view_tif_ch2, start=FOI[1], end=FOI[2])
  flimg1 <- flip(flimg1) # flip images to match fly-view
  flimg2 <- flip(flimg2) # flip images to match fly-view
  
  # flv exposure 18 ms
  FVOFFSET <- 0 
  frid     <- frid + FVOFFSET
  # Load fly-view camera images
  fvimgl <- dipr::readFMF(fly_view_fmf, frame=(frid)) 
  # Determined by ROI size
  fvsz   <- c(ROI[3],ROI[4])
  # Apply resize and translation to align with fluo-view
  fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear"), -center2)
  x1         <- (dim(fvimgl)[1] - fvsz[1])/2
  x2         <- x1 + fvsz[1] - 1
  # Match the size of the fvimgl and flimg by cropping (Assume size determined by flyview)
  fvimgl <- fvimgl[x1:x2,x1:x2,1:dim(fvimgl)[3]]
  
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
  ang <- ang_res[[1]] # in radianss
  centroid <- ang_res[[2]]
  markernum <- ang_res[[3]]
  
  # Find good frames
  angdiff <- c(0, diff(ang))
  ang_thresh <- 0.1
  goodangfr <- which(angdiff < ang_thresh & angdiff > -ang_thresh)
  goodmarkerfr <- which(markernum == 3)
  
  objdist <- sqrt((centroid[,1]-dim(fvimgl)[1]/2)^2 + (centroid[,2]-dim(fvimgl)[2]/2)^2)
  motion <- c(0, sqrt(diff(centroid[,1])^2 + diff(centroid[,2])^2))
  png(file=paste0(output_prefix, "_motion.png"), width=400, height=400)
  plot(motion)
  dev.off()
  goodmotionfr <- which(motion < motion_thresh)
  
  LoGkern <- round(dipr::LoG(9,9,1.4)*428.5)
  flimg2log <- EBImage::filter2(flimg2, LoGkern)
  centermask <- EBImage::drawCircle(flimg2[,,1]*0, dim(flimg2)[1]/2, dim(flimg2)[2]/2, 100, col=1, fill=T)
  flimg2cntlog <- dipr::ssweep(flimg2log, centermask, op="*")
  quantcnt <- apply(flimg2cntlog, 3, function(x) quantile(x, 0.9))
  png(file=paste0(output_prefix, "_quantcnt.png"), width=400, height=400)
  plot(quantcnt)
  dev.off()
  goodfocusfr <- which(quantcnt > focus_thresh & quantcnt < 8000)
  goodfr <- Reduce(intersect, list(goodmarkerfr, goodmotionfr, goodangfr, goodfocusfr))
  badix  <- badfr - (FOI[1] - 1)
  goodfr <- setdiff(goodfr,badix)
  goodfr <- setdiff(goodfr,(fluo_flash$flflashes - (FOI[1] - 1)))
  loggit::message(paste0("Good frames were ",paste0(goodfr,collapse = " ")))
  
  # Save good frames from each condition
  frix <- array(1,length(quantcnt))
  loggit::message(sprintf("Good Frames:\nmarker: %d\tmotion: %d\tangle: %d\tfocus: %d\n",
                          sum(frix[goodmarkerfr]),sum(frix[goodmotionfr]), sum(frix[goodangfr]), sum(frix[goodfocusfr])))
  
  # Save index of good frames
  if(FOI[1] != F) {
    write.table(cbind(1:length(goodfr),goodfr,goodfr + (FOI[1] - 1)),
                paste0(output_prefix, "_gfrid.csv"), sep = ",", row.names=T)
    saveRDS(goodfr + (FOI[1] - 1), paste0(output_prefix, "_gfrid.RDS"))
  } else {
    write.table(goodfr, paste0(output_prefix, "_gfrid.csv"), sep = ",", row.names=F)
    saveRDS(cbind(1:length(goodfr),goodfr), paste0(output_prefix, "_gfrid.RDS"))
  }
  
  # Load arena-view camera images
  avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
  EBImage::writeImage(avimgl/255, file=paste0(output_prefix, "_avimgl.tif"))
  #EBImage::writeImage(avimgl[,,goodfr]/255, file=paste0(output_prefix, "_avimgl_goodfr.tif"))
  rm(avimgl)
  
  # Apply rotation compensation
  loggit::message(paste0("Applying rotation compensation to the flyview video..."))
  rot <- fvimgl
  for (r in 1:dim(rot)[3]){
    rot[,,r] <- RNiftyReg::rotate(fvimgl[,,r], ang[r], anchor = c("center"))
  }
  
  # Flyview tracking error
  fverr <- read.table(paste0(dir, list.files(dir, pattern="^fv-traj-")), colClasses = "numeric")[(frid),4:5]
  
  # Use flyview trajectory file to find the error
  fverr <- as.matrix(fverr) - 120
  
  # Apply rotation to flyview error vectors
  rtr <- array(0,dim(fverr))
  for(i in 1:length(ang)) {
    rtr[i,] <- fverr[i,] %*% array(c(cos(ang[i]),sin(ang[i]),-sin(ang[i]),cos(ang[i])),c(2,2))
  }
  
  # Round rotated error to nearest pixel
  rtr    <- round(rtr)
  centers <- rtr
  
  # Offset compensation for bead/coverslip offset
  if(!is.na(ctr_offset[1]))
    centers <- t(t(centers) + ctr_offset)
  
  # Apply translation compensation
  loggit::message(paste0("Applying translation compensation to the flyview video..."))
  rottrans <- fvimgl[,,]
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
  loggit::message(paste0("Applying rotation compensation to the fluoview video..."))
  redrot <- flimg2[,,]
  for (rr in 1:dim(redrot)[3]){
    redrot[,,rr] <- as.Image(RNiftyReg::rotate(flimg2[,,rr], ang[rr], anchor = c("center")))
  }
  greenrot <- flimg1[,,]
  for (rg in 1:dim(greenrot)[3]){
    greenrot[,,rg] <- as.Image(RNiftyReg::rotate(flimg1[,,rg], ang[rg], anchor = c("center")))
  }
  
  loggit::message(paste0("Applying translation compensation to the fluoview video..."))
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
  
  ## Part 6. 
  # Image segmentation and fluorescence quantification ----
  # Normalize rotated imgs
  offs   <- as.integer(dim(redrottrans)[1] * (1 - 1/sqrt(2))) 
  redval <- redrottrans[(1+offs):(dim(redrottrans)[2]-offs),(1+offs):(dim(redrottrans)[2]-offs),]
  grnval <- greenrottrans[(1 + offs):(dim(greenrottrans)[2]-offs),(1+offs):(dim(greenrottrans)[2]-offs),]
  EBImage::writeImage(normalize(rottrans[(1+offs):(dim(rottrans)[2]-offs),(1+offs):(dim(rottrans)[2]-offs),goodfr], separate=F), file=paste0(output_prefix, "_rottrans100.tif")) 
  EBImage::writeImage(normalize(redval, separate=F, inputRange=input_range_r), file=paste0(output_prefix, "_redval.tif")) 
  EBImage::writeImage(normalize(grnval, separate=F, inputRange=input_range_g), file=paste0(output_prefix, "_grnval.tif")) 
  redval <- redval[,,goodfr]
  grnval <- grnval[,,goodfr]
  redval <- (redval - min(redval))/(max(redval) - min(redval))
  grnval <- (grnval - min(grnval))/(max(grnval) - min(grnval))
  
  # Dimensions / Number Frames
  wr <- dim(redval)[1]
  hr <- dim(redval)[2]
  fr <- dim(redval)[3]
  
  # ROI creation ----
  # If window size/offsets not passed do dialog
  if(is.na(window_size) || is.na(window_offset)) {
    # Interactively determine window size and offset to include neurons of interest
    num_rois <- as.integer(readline("How many ROIs to be added?: "))
  } else {
    if(length(window_size) != length(window_offset)) 
      stop("Length of window_size and window_offset must be the same")
    num_rois <- length(window_size)
  }
  
  # Initialize ROI masks with zero
  roimasks <- array(rep(FALSE,hr*wr*num_rois),c(hr,wr,num_rois))
  
  # Initialize ROI coords and windows sizes/offsets
  roiix <- array(rep(0,num_rois*4),c(num_rois,4))
  winoffs <- array(rep(0,num_rois*2),c(num_rois,2))
  winsize <- array(rep(0,num_rois*2),c(num_rois,2))
  
  for(i in 1:num_rois) {
    
    if(!(is.na(window_size) && is.na(window_offset))) {
      # Use passed parameters for size/offset
      winoffs[i,] <- window_offset[[i]]
      winsize[i,] <- window_size[[i]]
      ans <- c("Y","Y")
      
      # Check if ROI is valid
      win_valid <- all((winsize[i,]/2 + winoffs[i,]) <= wr/2)
      if(!win_valid)
        stop("Invalid window size/offset")
      
    } else {
      
      # Set ROI size at half crop size no offset
      winoffs[i,] <- c(0,0)
      winsize[i,] <- c(as.integer(wr/2),as.integer(hr/2))
      ans <- c("N","N")
      
    }
    
    # Initialize ROI with default
    roiix[i,] <- c((wr/2 + winoffs[i,1] - winsize[i,1]/2),
                   (wr/2 + winoffs[i,1] + winsize[i,1]/2),
                   (hr/2 + winoffs[i,2] - winsize[i,2]/2),
                   (hr/2 + winoffs[i,2] + winsize[i,2]/2))
    
    # Grab the Roi
    redwindowdisp <- EBImage::normalize(redval, separate=F)
    grnwindowdisp <- EBImage::normalize(grnval, separate=F)
    
    #redwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],]   <- redval[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],]
    #grnwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],] <- grnval[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],]
    
    # Draw box around ROI
    redwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,3],] <- grnwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,3],] <- 1
    redwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,4],] <- grnwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,4],] <- 1
    redwindowdisp[roiix[i,1],roiix[i,3]:roiix[i,4],] <- grnwindowdisp[roiix[i,1],roiix[i,3]:roiix[i,4],] <- 1
    redwindowdisp[roiix[i,2],roiix[i,3]:roiix[i,4],] <- grnwindowdisp[roiix[i,2],roiix[i,3]:roiix[i,4],] <- 1
    
    # Show both channels when selecting window/offset
    print(EBImage::display(abind(redwindowdisp,grnwindowdisp,along=2)))
    EBImage::writeImage(abind(redwindowdisp,grnwindowdisp,along=2), file=paste0(output_prefix, "_redwindow" ,i ,".tif")) 
    
    while(!all(stringr::str_to_lower(ans)=="y")) {
      
      print(sprintf("Current window_size for ROI %d is: x=%d y=%d", i, winsize[i,1], winsize[i,2]))
      print(sprintf("Current window_offset for ROI %d is: x=%d y=%d", i, winoffs[i,1], winoffs[i,2]))
      
      ans[1] <- readline("Check redwindow.tif. Is the window size good (Y or N)?:")
      if(ans[1] != "Y" && ans[1] != "y") {
        winsize[i,1] <- as.integer(readline("Enter new x size:"))
        winsize[i,2] <- as.integer(readline("Enter new y size:"))
      }
      ans[2] <- readline("Check redwindow.tif. Is the window offset good (Y or N)?:")
      if(ans[2] != "Y" && ans[2] != "y") {
        winoffs[i,1] <- as.integer(readline("Enter new x offset (positive to right):"))
        winoffs[i,2] <- as.integer(readline("Enter new y offset (positive to down):"))
      }
      
      # Check if selected ROI is valid else set to default
      win_valid <- all((winsize[i,]/2 + abs(winoffs[i,])) <= wr/2)
      if(!win_valid)
        loggit::message("Invalid window size/offset")
      
      ## Update reference
      if(!all(stringr::str_to_lower(ans)=="y") & win_valid) {
        
        # Update ROI
        roiix[i,] <- c((wr/2 + winoffs[i,1] - winsize[i,1]/2),
                       (wr/2 + winoffs[i,1] + winsize[i,1]/2),
                       (hr/2 + winoffs[i,2] - winsize[i,2]/2),
                       (hr/2 + winoffs[i,2] + winsize[i,2]/2))
        
        # Gamma correction
        redwindowdisp <- EBImage::normalize(redval, separate=F)
        grnwindowdisp <- EBImage::normalize(grnval, separate=F)
        
        # Draw box around ROI
        redwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,3],] <- grnwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,3],] <- 1
        redwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,4],] <- grnwindowdisp[roiix[i,1]:roiix[i,2],roiix[i,4],] <- 1
        redwindowdisp[roiix[i,1],roiix[i,3]:roiix[i,4],] <- grnwindowdisp[roiix[i,1],roiix[i,3]:roiix[i,4],] <- 1
        redwindowdisp[roiix[i,2],roiix[i,3]:roiix[i,4],] <- grnwindowdisp[roiix[i,2],roiix[i,3]:roiix[i,4],] <- 1
        
        print(EBImage::display(abind(redwindowdisp,grnwindowdisp,along=2)))
        EBImage::writeImage(abind(redwindowdisp,grnwindowdisp,along=2), file=paste0(output_prefix, "_redwindow" ,i ,".tif")) 
      }
    } # end selecting ROI i
    
    roimasks[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],i] = 1
  }
  
  # Aggregate all ROI masks
  roimask <- rowSums(roimasks,dims=2)
  loggit::message(paste0("ROIs created."))
  
  # Clean up
  rm(flimg1)
  rm(flimg2)
  rm(flimg2cntlog)
  
  # Mask refinement ----
  # Preallocate Segment Masks and Masked Image
  rroithr  <- greenrois <- redrois <- array(rep(0,wr*hr*fr),c(wr,hr,fr))
  seg_mask <- array(rep(0,wr*hr*fr),c(wr,hr,fr,num_rois)) # Save mask for each ROI
  # Add regions to mask
  for(i in 1:num_rois) {
    
    #[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],]                     # ROI in Valid
    #[roiix[i,1]+offs:roiix[i,2]-offs,roiix[i,3]+offs:roiix[i,4]-offs,] # ROI in Original
    
    rroi     <- redrottrans[(roiix[i,1]+offs):(roiix[i,2]+offs),(roiix[i,3]+offs):(roiix[i,4]+offs),goodfr] 
    groi     <- greenrottrans[(roiix[i,1]+offs):(roiix[i,2]+offs),(roiix[i,3]+offs):(roiix[i,4]+offs),goodfr]
    rroimed  <- EBImage::medianFilter(rroi/2^16, size=3)
    groimed  <- EBImage::medianFilter(groi/2^16, size=3)
    redrois[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],]   <- rroimed
    greenrois[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],] <- groimed
    
    seg_mask_win <- array(0,dim(rroimed))
    toff <- (apply(rroimed,MARGIN=3,max) - apply(rroimed,MARGIN=3,mean))*0.4
    
    for(j in 1:fr) {
      seg_mask_win[,,j] <- EBImage::thresh(rroimed[,,j],
                                           w=as.integer(winsize[i,1]/2),
                                           h=as.integer(winsize[i,2]/2),
                                           offset=toff[j])
    }
    
    # Morphological Operations (Opening + Fill)
    seg_mask_win <- EBImage::erode(seg_mask_win,kern=makeBrush(3,shape="diamond"))    
    seg_mask_win <- EBImage::dilate(seg_mask_win,kern=makeBrush(3,shape="diamond"))    
    seg_mask_win <- EBImage::fillHull(seg_mask_win)
    
    # Add segmented ROI to overall mask
    seg_mask[roiix[i,1]:roiix[i,2],roiix[i,3]:roiix[i,4],,i] <- seg_mask_win
  }
  loggit::message(paste0("ROI masks created."))
  
  seg_maski <- seg_mask # Each ROI's mask frames
  
  # Allocate data frame and means
  rawintroi   <- data.frame(frame=1:length(frid))  
  redave      <- vector()
  grnave      <- vector()
  ratioave    <- vector()
  greenperred <- array(0,dim(redrois))
  seg_masku   <- array(0,dim(redrois)) # Total untion segmentation mask
  
  # For each ROI process separately
  for(r in 1:num_rois) {
    
    seg_mask <- seg_maski[,,,r]      # segmentation mask
    roimask  <- roimasks[,,r]        # ROI mask
    roimask  <- as.vector(roimask)
    redroi   <- redrois*roimask      # pixels in ref channel ROI
    grnroi   <- greenrois*roimask    # pixels in cal channel ROI
    redseg   <- redrois*seg_mask     # ref pixels segmentation mask video
    grnseg   <- greenrois*seg_mask   # cal pixels in segmentation mask video
    
    # Allocate for per frame baseline metrics
    minsred <- array(0,fr)
    meanqred <- array(0,fr)
    quantred <- array(0,fr)
    
    for(i in 1:fr) {
      
      # Per Frame Pixels in entire ROI
      redvalpx      <- redroi[,,i][roimask > 0]
      #redvalpx     <- redseg[,,i][seg_mask[,,i] > 0]
      
      # Per frame quantile (bgratio percentile)
      quantred[i]  <- quantile(redvalpx,bgratio)
      
      # Per frame mean of lowest (1 - bgratio) pixels
      meanqred[i]  <- mean(redvalpx[redvalpx < quantred[i]])
      
    }
    
    # Filter red channel pixels below mean of
    # of the bgratio quantile across frames
    bl <- mean(quantred,na.rm=TRUE)
    seg_mask[redseg <= bl] <- 0
    
    # Morphilogical Operations:
    # Label regions and apply area filtering
    shape_metric <- 1 #s.area s.perimeter s.radius.mean s.radius.sd s.radius.min s.radius.max
    thrsh_map <- apply(seg_mask,3,bwlabel)
    thrsh_map <- array(thrsh_map,dim=dim(seg_mask))
    maskprops <- apply(thrsh_map,3,function(x) list(computeFeatures.shape(x)))
    maskprops <- lapply(maskprops, "[[", 1)
    objrmidx  <- lapply(maskprops,FUN=function(x) which(x[,shape_metric] <= size_thresh))
    seg_mask  <- rmObjects(thrsh_map, objrmidx)
    loggit::message(paste0("Removed objects smaller than ", size_thresh, " pixels in ROI ", r,"."))
    
    # Return filtered regions to binary mask
    seg_mask[seg_mask > 0] <- 1
    
    # Mask pixels around margin of segmented area for background subtraction
    bg_mask  <- (EBImage::dilate(seg_mask,kern=makeBrush(15,shape="diamond")) 
                 - EBImage::dilate(seg_mask,kern=makeBrush(7,shape="diamond")))
    
    # Constrain bg pixels to roi
    bg_mask[!(bg_mask & array(roimask,dim(bg_mask)))]  <- 0
    
    # Get mean background value for each channel per frame
    bggrn    <- grnroi * bg_mask
    bgred    <- redroi * bg_mask
    bgarea   <- apply(bg_mask,MARGIN=3,sum)
    bggrnave <- apply(bggrn,MARGIN=3,sum)/bgarea
    bgredave <- apply(bgred,MARGIN=3,sum)/bgarea
    
    # Background subtraction per frame
    for(i in 1:fr) {
      redseg[,,i]   <- redroi[,,i] - bgredave[i]
      grnseg[,,i]   <- grnroi[,,i] - bggrnave[i]
    }
    
    # Filter values below background
    redseg[redseg < 0] <- 0.0
    grnseg[grnseg < 0] <- 0.0
    
    #Background subtracted Channels/Ratio Image
    redseg      <- redseg*seg_mask
    grnseg      <- grnseg*seg_mask
    ratioseg    <- (grnseg/redseg)
    
    # Filter infinite x/0 
    ratioseg[!is.finite(ratioseg)] <- 0.0
    
    # Filter out lowest ratio values in case of only subset of labeled neurons active
    qfcutoff <- quantile(ratioseg[seg_mask>0],ratiocutoff)
    for(i in 1:fr) {
      # Filter pixel ratios below quantile cutoff
      seg_mask[,,i][ratioseg[,,i] < qfcutoff] <- 0
    }
    
    # Final segmentation mask
    loggit::message(paste0("Segmentation finished for ROI ",r,"."))
    
    # Write mask for this ROI
    EBImage::writeImage(seg_mask, file=paste0(output_prefix, "_segmask_",r,".tif")) 
    
    # Add segmented region to complete mask
    seg_masku[seg_mask > 0] <- 1
    
    # Apply final segmentation mask
    redseg   <- redseg*seg_mask
    grnseg   <- grnseg*seg_mask
    ratioseg <- (grnseg/redseg)
    
    # Filter infinite x/0 
    ratioseg[!is.finite(ratioseg)] <- 0.0
    
    # Ratio video
    greenperred <- greenperred + ratioseg
    
    # Mean of each channel in mask
    segarea   <- apply(seg_mask,MARGIN=3,sum)
    redavei   <- apply(redseg,MARGIN=3,sum)/segarea
    grnavei   <- apply(grnseg,MARGIN=3,sum)/segarea
    ratioavei <- apply(ratioseg,MARGIN=3,sum)/segarea
    
    # Means across both
    redave   <- cbind(redave,redavei)
    grnave   <- cbind(grnave,grnavei)
    ratioave <- cbind(ratioave,ratioavei)
    
    # Write intensities and ratio for this ROI
    redout <- grnout <- ratout <- array(0,length(frid))
    redout[goodfr] <- redavei
    grnout[goodfr] <- grnavei
    ratout[goodfr] <- ratioavei
    
    tmp        <- data.frame(redout,grnout,ratout)
    names(tmp) <- c(paste0("redave",r),paste0("grnave",r),paste0("ratioave",r))
    rawintroi  <- cbind(rawintroi, tmp)
    
  } # end per ROI loop
  
  loggit::message(paste0("Segmentation finished."))
  
  # Write per roi intensity file
  write.table(rawintroi, paste0(output_prefix, "_rawintroi.csv"), sep = ",", row.names=F)
  
  # Mean across ROIs (FOI length)
  redave         <- rowMeans(redave,na.rm=TRUE)
  greenave       <- rowMeans(grnave,na.rm=TRUE)
  greenperredave <- rowMeans(ratioave,na.rm=TRUE)
  
  # Raw intentsities all
  ratrawall <- redrawall <- grnrawall <- array(NA,length(frid))
  ratrawall[goodfr] <- greenperredave
  redrawall[goodfr] <- redave
  grnrawall[goodfr] <- greenave
  
  # Segmentation mask across ROIs
  seg_mask       <- seg_masku
  
  # Check for finite ratio values (reduced lenght no NA)
  goodfrratidx   <- is.finite(greenperredave)
  greenperredave <- greenperredave[goodfrratidx]
  goodfrrat      <- goodfr[goodfrratidx]
  redave         <- redave[goodfrratidx]
  greenave       <- greenave[goodfrratidx]
  greenperred    <- greenperred[,,goodfrratidx]
  seg_mask       <- seg_mask[,,goodfrratidx]
  
  # TODO: Raw result placeholder 
  intensity <- zoo::rollmean(greenperredave, 3, align="left")
  
  # Temp copy of ratio image for heatmap
  gprimage <- greenperred
  gprimage[gprimage >= 1]         <- 0.99 # Threshold ratios > 1 for heatmap
  
  # Full length segmentation mask (0 for frames not analyzed)
  seg_mask_all <- array(0, dim=c(dim(seg_mask)[c(1,2)],dim(rottrans)[3]))
  for(sfr in 1:length(goodfrrat)) {
    seg_mask_all[,,goodfrrat[sfr]] <- seg_mask[,,sfr]
  }
  
  # Create heatmap image
  grratiocolor <- array(0, dim=c(dim(gprimage)[c(1,2)], 3, dim(rottrans)[3]))
  for(cfr in 1:dim(gprimage)[3]){
    gcfr = goodfrrat[cfr]
    grratiocolor[,,,gcfr] <- dipr::pseudoColor(gprimage[,,cfr], colorRange[1], colorRange[2])
  }
  
  grratiocolor <- Image(grratiocolor, colormode="Color")
  EBImage::writeImage(grratiocolor, file=paste0(output_prefix, "_grratiocolor.tif"))
  rm(gprimage)
  
  # Overlay fly_view and F_ratio image
  rottransmask <- array(0, dim=c(dim(rottrans)[c(1,2)], dim(rottrans)[3]))
  rottransmask[(1+offs):(dim(redrottrans)[2]-offs),(1+offs):(dim(redrottrans)[2]-offs),] <- seg_mask_all
  
  rottranscolor <- array(0, dim=c(dim(rottrans)[c(1,2)], 3, dim(rottrans)[3]))
  rottranscolor[,,3,] <-rottranscolor[,,2,] <-rottranscolor[,,1,] <- rottrans/255*(1-rottransmask)
  
  grratiocolorl <- rottranscolor*0
  grratiocolorl[(1+offs):(dim(redrottrans)[2]-offs),(1+offs):(dim(redrottrans)[2]-offs),,] <- grratiocolor
  
  # Create and add marker to denote frames not analyzed in output video
  mrk_img <- array(0,dim(grratiocolorl)[c(1,2)])
  mrk_dim <- 19
  mrk_ctr <- (mrk_dim-1)/2
  mrk_r1  <- 7
  mrk_r2  <- 5
  mrk_sym <- array(0,c(mrk_dim,mrk_dim))
  mrk_sym <- drawCircle(mrk_sym,mrk_ctr,mrk_ctr,mrk_r1,1,fill=T)
  mrk_sym <- drawCircle(mrk_sym,mrk_ctr,mrk_ctr,mrk_r2,0,fill=T) 
  
  xw <- 1
  rx <- as.integer(sqrt(2)/2*mrk_r1 - xw)
  start <- mrk_ctr - rx
  end   <- mrk_ctr + rx
  
  for(ix in start:end)
    mrk_sym[ix,(ix-xw):(ix+xw)] = 1
  
  mrk_sym <- apply(mrk_sym,1,rev)
  mrk_img[(dim(mrk_img)[2]-mrk_dim+1):(dim(mrk_img)[2]),1:mrk_dim]<-mrk_sym
  
  for(gf in setdiff(1:dim(grratiocolorl)[4],goodfrrat)) {
    grratiocolorl[,,1,gf] <- grratiocolorl[,,1,gf] + mrk_img
  }
  
  if(stimulus == T){
    stim_mrk_img <- array(0,dim(grratiocolorl)[c(1,2)])
    stim_mrk_dim <- 19
    stim_mrk_ctr <- (stim_mrk_dim-1)/2
    stim_mrk_r  <- 5
    stim_mrk_sym <- array(0,c(stim_mrk_dim,stim_mrk_dim))
    stim_mrk_sym <- drawCircle(stim_mrk_sym,stim_mrk_ctr,stim_mrk_ctr,stim_mrk_r,1,fill=T)
    
    stim_mrk_img[1:stim_mrk_dim,1:stim_mrk_dim]<-stim_mrk_sym
    
    for(ef in which(event_pattern==1)) {
      grratiocolorl[,,2,ef] <- grratiocolorl[,,2,ef] + stim_mrk_img
    }
  }
  
  rm(rottranscolor)
  
  # overlay red channel and F_ratio color image
  redrottranscol <- array(0, dim=c(dim(redrottrans)[c(1,2)], 3, dim(redrottrans)[3]))
  redrottranscol[,,3,] <-redrottranscol[,,2,] <-redrottranscol[,,1,] <- redrottrans*(1-rottransmask)
  redrottranscol <- normalize(redrottranscol, separate=F, inputRange=input_range_r)
  redcolor <- redrottranscol + grratiocolorl
  redcolor <- Image(redcolor, colormode="Color")
  rm(redrottranscol)
  rm(rottransmask)
  
  
  fr_width  = ROI[3]
  fr_height = ROI[4]
  # Create side-by-side view of fly_view and fluo_view images
  frgcombined <- array(dim=c(dim(rottrans)[1]*4, dim(rottrans)[2], 3, dim(rottrans)[3]))
  frgcombined[1:fr_width,1:fr_height,3,1:dim(rottrans)[3]] <- 
    frgcombined[1:fr_width,1:fr_height,2,1:dim(rottrans)[3]] <- 
    frgcombined[1:fr_width,1:fr_height,1,1:dim(rottrans)[3]] <- 
    normalize(rottrans, separate=F)
  frgcombined[(fr_width+1):(2*fr_width),1:fr_height,3,1:dim(redrottrans)[3]] <-
    frgcombined[(fr_width+1):(2*fr_width),1:fr_height,2,1:dim(redrottrans)[3]] <-
    frgcombined[(fr_width+1):(2*fr_width),1:fr_height,1,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F, inputRange=input_range_r)
  frgcombined[(2*fr_width+1):(3*fr_width),1:fr_height,3,1:dim(greenrottrans)[3]] <-
    frgcombined[(2*fr_width+1):(3*fr_width),1:fr_height,2,1:dim(greenrottrans)[3]] <-
    frgcombined[(2*fr_width+1):(3*fr_width),1:fr_height,1,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F, inputRange=input_range_g)
  frgcombined[(3*fr_width + 1):(4*fr_width),1:fr_height,,1:dim(redrottrans)[3]] <- redcolor
  frgcombined <-  Image(frgcombined, colormode="Color")
  
  redcolor100 <- redcolor[(1 + offs):(dim(greenrottrans)[2]-offs),(1+offs):(dim(greenrottrans)[2]-offs),,]
  EBImage::writeImage(redcolor100, file=paste0(output_prefix, "_redcolor100.tif"))
  EBImage::writeImage(normalize(redrottrans, separate=F, inputRange=input_range_r), file=paste0(output_prefix, "_redrottrans.tif"))
  rm(redrottrans)
  EBImage::writeImage(normalize(greenrottrans, separate=F, inputRange=input_range_g), file=paste0(output_prefix, "_greenrottrans.tif"))
  rm(greenrottrans)
  EBImage::writeImage(frgcombined, file=paste0(output_prefix, "_frgcombined_goodfr20_normalized.tif"))
  rm(frgcombined)
  
  # Quantification of fluorescence intensity ----
  
  # Mean Red/Green/Ratio (Analyzed Frames)
  datrawint <- data.frame(x=goodfrrat, y=greenperredave, r=redave, g=greenave)
  
  # LOESS Model
  datloessint <- loess(y ~ x, data=datrawint, span=loess_span, control=loess.control(surface="direct"))
  redloessint <- loess(r ~ x, data=datrawint, span=loess_span)
  greenloessint <- loess(g ~ x, data=datrawint, span=loess_span)
  # LOESS Predictions
  datsmoothint <- data.frame(x=goodfrrat, y=predict(datloessint))
  redsmoothint <- data.frame(x=goodfrrat, y=predict(redloessint))
  greensmoothint <- data.frame(x=goodfrrat, y=predict(greenloessint))
  
  #LOESS prediction of all time points
  datsmoothintall <- data.frame(x=1:length(frid), y=predict(datloessint, 1:length(frid)))
  redsmoothintall <- predict(redloessint, 1:length(frid))
  grnsmoothintall <- predict(greenloessint, 1:length(frid))
  
  png(file=paste0(output_prefix, "_datsmoothint.png"), width=400, height=400)
  par(mar = c(5,5,2,5))
  plot(datrawint$x, datrawint$y, type="p", pch=16, ylab="F_ratio", xlab="frame", cex=0.5)
  lines(datsmoothint, col="blue")
  par(new=TRUE)
  plot(redsmoothint, col="red", axes=F, xlab=NA, ylab=NA)
  par(new=TRUE)
  plot(greensmoothint, col="green", axes=F, xlab=NA, ylab=NA)
  axis(side = 4)
  mtext(side = 4, line = 3, 'Fluorescence intensity')
  dev.off()
  
  
  # Save a frame map for the sequence
  if(FOI[1] != F) {
    # Write out frame map:
    flframe           <- FOI[1]:FOI[2]
    good              <- array(0,length(frid))
    good[goodfrrat]   <- 1
    fmap              <- cbind(flframe,frid,frida,good)
    colnames(fmap)    <- c("fluo_frame","fly_frame","arena_frame","good_frame")
    write.csv(fmap,paste0(output_prefix,"_frame_map.csv"),row.names = F)
  } else {
    flframe           <- 1:length(frid)
    good              <- array(0,length(frid))
    good[goodfrrat]   <- 1
    fmap              <- cbind(flframe,frid,frida,good)
    colnames(fmap)    <- c("fluo_frame","fly_frame","arena_frame","good_frame")
    write.csv(fmap,paste0(output_prefix,"_frame_map.csv"),row.names = F)
  }
  
  
  saveRDS(datrawint, paste0(output_prefix, "_datrawint.RDS"))
  write.table(datrawint, paste0(output_prefix, "_datrawint.csv"), sep = ",", row.names=F)
  
  # Dataframe LOESS predications
  saveRDS(datsmoothint, paste0(output_prefix, "_datloessint.RDS")) 
  
  # Behavior analysis ----
  ## Analyze trajectories
  trj_res <- analyze_trajectories(dir=dir,
                                  output=output_prefix,
                                  fpsfv=syncing$fpsfv,
                                  interaction=interaction)
  if(interaction==T){
    loggit::message("Detecting interaction")
    closefr <- which(trj_res$flydist[frida] < dist_thresh)
    closefrid <- sapply(closefr, function(x) which.min(abs(syncing$frida-x)))
    write.table(closefrid, paste0(output_prefix, "_closefrid.txt"))
    
    ## Behavior analysis. This part is relevant only when interaction=T
    ## Calculate two lines from the tracked fly: one along the body axis, and the other to the other fly
    ## The angle between the two lines is the view angle of the fly
    # body axis is already given as "ang" in radians but plus pi/2
    # so all we need is the slope of the line between two flies
    
    # Determine which fly in the arena-view is being tracked by fly-view
    
    fly1trjfv <- trj_res$trjfv[frid,]
    fly1trjavmm <- trj_res$trja[frida,1:2]
    fly2trjavmm <- trj_res$trja[frida,3:4]
    
    # Get angles
    # Angle in degrees facing toward Fluoview
    angdf1  <- ang*180/pi - 15 + 90
    
    # Normalize fv angle to -180 - 180
    angdf1n       <- (angdf1 + 360) %% 360
    ixf1          <- angdf1n > 180 & !is.na(angdf1n)
    angdf1n[ixf1] <- angdf1n[ixf1] - 360 # FIX THE OFFSET:
    
    # Vector from male to female
    f1f2vec       <- fly2trjavmm - fly1trjavmm
    
    # Norm of fly1 to fly vector (distance)
    f1f2norm      <- sqrt(rowSums(f1f2vec^2))
    
    # x,y components of unit vector from fly1 to fly2
    f1f2x         <- f1f2vec[,1]/f1f2norm
    f1f2y         <- f1f2vec[,2]/f1f2norm
    
    # Angle vector from fly1 to fly2
    angdf1f2      <- atan2(f1f2y,f1f2x) * 180/pi
    
    # Angle between fly1 heading and fly2 centroid
    theta         <- -angdf1 + angdf1f2
    
    # Normalize fv angle to -180 - 180
    theta         <- (theta + 360) %% 360
    ixt           <- theta > 180 & !is.na(theta)
    theta[ixt]    <- theta[ixt] - 360

    #vecB <- data.frame(Bx=(fly2trja[,1] - fly1trja[,1]), By=(-fly2trja[,2] + fly1trja[,2]))
    # theta ＝ atan2(AxB，A*B) in radian
    # For now left side of the view is positive
    #theta <- 180/pi*atan2((vecA[,1]*vecB[,2] - vecA[,2]*vecB[,1]), (vecA[,1]*vecB[,1] + vecA[,2]*vecB[,2]))
    
  }else{
    if(length(trj_res$trja)==0){
      loggit::message("No valid trajectories found")
      fly1trja  <- data.frame(xr=rep(0, length(ang)), yr=rep(0, length(ang)))
      fly1trjfv <- data.frame(xr=rep(0, length(ang)), yr=rep(0, length(ang)))
      theta     <- array(0,length(ang))
    }else{ # Single Fly
      fly1trja  <- trj_res$trja[frida,1:2]
      fly1trjfv <- trj_res$trjfv[frid,1:2]
      theta     <- array(0,length(ang))
    }
  }
  
  if(gen_av_trj_vid) {
    # Save arenaview tracked video
    avimg   <- dipr::readFMF(arena_view_fmf, frames=frida)/255
    nflies  <- dim(trj_res$trja)[2]/2
    
    # Read arena view trajectory file
    trja   <- read.table(paste0(dir, list.files(dir, pattern="^av-traj-")))
    
    # Allocate first channel and rgb video file
    avimgc  <-array(NA,c(dim(avimg)[c(1,2)],3,dim(avimg)[3]))
    avimgc3 <- avimg
    
    if(nflies == 1) {
      
      fly1_id  <- 0
      fly1_col <- c(2:3) + (fly1_id)*3
      avpix1   <- trja[frida,fly1_col]
      
      for(i in 1:dim(avimg)[3]) {
        
        avimgc3[,,i]   <- drawCircle(avimgc3[,,i],avpix1[i,1],avpix1[i,2],3,col=1,fill=T)
        
      }
      avimgc[,,1,] <- avimg
      avimgc[,,2,] <- avimg
      avimgc[,,3,] <- avimgc3
      avimgc <- Image(avimgc, colormode="Color")
      EBImage::writeImage(avimgc, file=paste0(output_prefix, "_arena_track.tif"))
      
      rm(avimg)
      rm(avimgc)
      rm(avimgc3)
      
    } else if(nflies == 2) {
      #
      fly1_col <- c(2:3) + (fly1_id)*3
      fly2_col <- c(2:3) + (!fly1_id)*3
      
      # Allocate 2nd fly channel
      avimgc1  <- avimgc3
      
      avpix1 <- trja[frida,fly1_col]
      avpix2 <- trja[frida,fly2_col]
      
      for(i in 1:dim(avimg)[3]) {
        
        avimgc3[,,i]   <- drawCircle(avimgc3[,,i],avpix1[i,1],avpix1[i,2],3,col=1,fill=T)
        avimgc1[,,i]   <- drawCircle(avimgc1[,,i],avpix2[i,1],avpix2[i,2],3,col=1,fill=T)
        
      }
      
      avimgc[,,1,] <- avimgc1
      avimgc[,,2,] <- avimg
      avimgc[,,3,] <- avimgc3
      avimgc <- Image(avimgc, colormode="Color")
      EBImage::writeImage(avimgc, file=paste0(output_prefix, "_arena_track.tif"))
      
      rm(avimg)
      rm(avimgc)
      rm(avimgc1)
      rm(avimgc3)
      
    }
  }
  
  
  ## Plotting ----
  if(is.na(baseline)) {
    # Default to using first frame
    F0int    <- intensity[1]
    F0loess  <- datsmoothintall[1,2]
  } else {
    # Specify -1 to use min (Positive Change)
    if(baseline == -1) {
      F0int    <- min(intensity)
      F0loess  <- min(datsmoothintall[,2])
    }else if(baseline == 0) {
      F0int    <- mean(intensity)
      F0loess  <- mean(datsmoothintall[,2])
    } else {
      # Specified a baseline frame or range of frames for mean
      idx      <- intersect(baseline[1]:baseline[length(baseline)],goodfr)
      if(is_empty(idx)) {
        message('Baseline frame is not analyzed... Using first frame.')
        F0int    <- intensity[1]
        F0loess  <- datsmoothintall[1,2]
      } else {
        F0int    <- mean(intensity[baseline[1]:baseline[length(baseline)]])
        F0loess  <- datsmoothintall[1,2]
      }
    }
  }
  
  deltaFint <- intensity - F0int
  dFF0int <- deltaFint/F0int * 100
  
  # Including bad frames
  dFF0intall <- (datsmoothintall[,2]-F0loess)/F0loess*100
  # Output all interesting data to output
  
  goodix            <- array(0,length(frida)) # Generate list of good frames
  goodix[goodfrrat] <- 1
  
  datdFF0all <- data.frame(n=1:length(frida),
                           goodfrix=goodix,
                           flframe=flframe,
                           avframe=frida,
                           fvframe=frid,
                           ratrawall=ratrawall,
                           redrawall=redrawall,
                           grnrawall=grnrawall,
                           fratloess=datsmoothintall[,2],
                           fredloess=redsmoothintall,
                           fgrnloess=grnsmoothintall,
                           dFFloess=dFF0intall,
                           f1f2dist=trj_res$flydist[frida],
                           ang=ang,
                           f1f2angle=theta,
                           dist2=f1f2norm,
                           row.names='n') # TODO: Add trajectories / behavior
  
  write.csv(datdFF0all, paste0(output_prefix, "_datdFF0all.csv"))
  
  datdFF0 <- data.frame(n=goodfrrat[1:(length(goodfrrat)-2)], f=dFF0int, 
                        d=trj_res$flydist[frida[goodfrrat[1:(length(goodfrrat)-2)]]],
                        a=theta[1:(length(goodfrrat)-2)])
  
  p1 <- ggplot2::ggplot(data=datdFF0, ggplot2::aes(x=n, y=f)) +
    ggplot2::geom_smooth(method="loess", span = loess_span, level=0.95) +
    ggplot2::ggsave(filename = paste0(output_prefix, "_dFF0int.pdf"), width = 8, height = 8)
  
  p2 <- ggplot2::ggplot(data=datdFF0, ggplot2::aes(x=n, y=d)) +
    ggplot2::geom_line() +
    ggplot2::ylim(0, 100)
  
  p3 <- ggplot2::ggplot(data=datdFF0, ggplot2::aes(x=n, y=a)) +
    ggplot2::geom_line() +
    ggplot2::ylim(-180, 180)
  
  png(file=paste0(output_prefix, "_dFF0_dist_angle.png"), width=400, height=400)
  Rmisc::multiplot(p1, p2, p3, cols=1)
  dev.off()
  
  # Create trajectory of the flies
  message("Creating trajectory of the flies...")
  
  df1 <- cbind(datdFF0all, fly1trjfv)
  
  if(interaction==T){
    event_pattern <- rep(1, nrow(df1))
    event_pattern[closefr] <- 1.5
  }else{
    event_pattern <- rep(1, nrow(df1))
  }
  
  df1 <- cbind(df1, event_pattern)
  
  if (interaction==T){
    df2 <- cbind(datdFF0all, fly2trjavmm)
    
    p4 <- ggplot2::ggplot(data=df2, ggplot2::aes(x=10*xr, y=10*yr)) + 
      geom_path(linetype=2, lwd = 1, color=1) +
      geom_path(data=df1,  ggplot2::aes(x=10*xr, y=10*yr, color=fratloess), linetype=1, lwd = 1) +
      coord_fixed(ratio = 1) +
      scale_x_continuous(limits=c(-240, 240), expand=c(0,0)) +
      scale_y_reverse(limits=c(220, -220), expand=c(0,0)) +
      scale_colour_gradientn(limits=c(40, max(df1$fratloess)), colours = c("blue", "red")) +
      ggforce::geom_ellipse(aes(x0 = 0, y0 = 0, a = 11.0795*20, b = 10*20, angle = 0)) + # Add an ellipse
      theme(line = element_blank(),
            text = element_blank(),
            title = element_blank(),
            legend.position="none",
            rect= element_blank(),
            plot.margin=unit(c(0,0,-1,-1),"lines"))
  }else{
    
    p4 <- ggplot2::ggplot(data=df1, ggplot2::aes(x=10*xr, y=10*yr, color=fratloess)) + 
      geom_path(linetype=1, lwd = event_pattern, linejoin="round", lineend="round") +
      coord_fixed(ratio = 1) +
      scale_x_continuous(limits=c(-240, 240), expand=c(0,0)) +
      scale_y_reverse(limits=c(220, -220), expand=c(0,0)) +
      scale_colour_gradientn(limits=c(40, max(df1$fratloess)), colours = c("blue", "red"), na.value = "gray70") +
      #scale_alpha(0.5) +
      ggforce::geom_ellipse(aes(x0 = 0, y0 = 0, a = 11.0795*20, b = 10*20, angle = 0), color="black") + # Add an ellipse
      theme(line = element_blank(),
            text = element_blank(),
            title = element_blank(),
            legend.position="none",
            rect= element_blank(),
            plot.margin=unit(c(0,0,-1,-1),"lines"))
    
  }
  
  pdf(file= paste0(output_prefix, "_trackResult_fv.pdf"), width = 4.4, height = 4, bg = "white")
  print(p4)
  dev.off()
  
  ## Output summary ----
  # Format string for multi ROI window sizse/offsets
  wins_str = "list("
  offs_str = "list("
  for(i in 1:num_rois) {
    wins_str = paste0(wins_str,"c(",winsize[i,1],",",winsize[i,2],")")
    offs_str = paste0(offs_str,"c(",winoffs[i,1],",",winoffs[i,2],")")
    if(i < num_rois) {
      wins_str = paste0(wins_str,",")
      offs_str = paste0(offs_str,",")
    } else {
      wins_str = paste0(wins_str,")")
      offs_str = paste0(offs_str,")")
    }
  }
  
  loggit::message(sprintf("Number of ROIs was %d", num_rois))
  loggit::message(sprintf("window size(s): %s", wins_str))
  loggit::message(sprintf("window offset(s): %s", offs_str))
  loggit::message(sprintf("FOI was from %d to %d",  FOI[1], FOI[2])) 
  loggit::message(paste0("Max F_ratio intensity in this bout was ", max(intensity)))
  loggit::message(paste0("Max F_ratio smoothed intensity in this bout was ", max(datsmoothint$y)))
  loggit::message(paste0("Number of good frames was ", length(goodfrratidx)))
  
  out_str = sprintf("||c(%d, %d) ||%s ||%s ||%d ||%.3f/%.3f ||%.3f/%.3f ||", 
                    FOI[1], FOI[2], wins_str, offs_str, length(goodfrratidx),
                    min(intensity), max(intensity), min(datsmoothint$y), max(datsmoothint$y))
  loggit::message(out_str)
  
  ## Convert fmf to tif format ----
  
  if(fmf2tif==T){
    dipr::fmf2tif(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), skip=10)
    dipr::fmf2tif(paste0(dir, list.files(dir, pattern="^av.*fmf$")), skip=2)
  }
  loggit::message("Cleaning up...")
  closeAllConnections()
  rm(list=setdiff(ls(), "out_str"))
  gc()
  loggit::message("Finished processing!")
  return(out_str)
}
