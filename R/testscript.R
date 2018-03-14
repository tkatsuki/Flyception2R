install.packages(c("devtools", "ggplot2", "RNiftyReg"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocInstaller", "EBImage"))
library(devtools)
devtools::install_github("tkatsuki/Flyception2R")
library(Flyception2R)

dir <- "H:/P1_GCaMP6s_tdTomato_02202018/P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy/"  # Don't forget the slash at the end
prefix <- "P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy"       # Will be used as a filename prefix
autopos <- T             # True if you want to align cameras automatically 
reuse <- F               # True if you want to reuse intermediate RDS files
fmf2tif <- T             # True if you want to convert fmf 
zoom <- 1             # Zoom ratio: fluo-view/fly-view
FOI <- c(1000, 1200)                 # A vector specifying start and end frame (e.g. c(10,1000)). False if you want to analyze all frames.
ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
binning <- 1             # Binning of the fluo-view camera
fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
interaction <- T         # True if you want to analyze fly-fly interaction
dist_thresh <- 4         # Threshold for detecting fly-fly interaction based on distance
rotate_camera <- -180    # Rotation angle needed to align fluo-view and fly-view


rlogging::SetLogFile(base.file=paste0(prefix, "_log.txt"), folder=dir)
message(dir)

output_prefix <- paste0(dir, prefix)
fluo_view_tif <- paste0(dir, list.files(dir, pattern="Pos0\\.ome\\.tif$"))
fluo_view_tif_ch1 <- paste0(dir, list.files(dir, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))
fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))


## First crop the fluo-view image with one ROI and find the flash frame
imageJ_crop_append(dir, ch=1, roi=ROI) # x and y coordinates of the top left corner, width, height

message("Detecting flash in fluo-view")
fluo_flash <- detect_flash(input=fluo_view_tif_ch1,
                           type="fluo",
                           output=output_prefix,
                           flash_thresh=fluo_flash_thresh,
                           reuse=reuse)
message("Detecting flash in fly-view")
fly_flash <- detect_flash(input=fly_view_fmf,
                          type="fly",
                          output=output_prefix,
                          flash_thresh=fv_flash_thresh,
                          reuse=reuse)
message("Detecting flash in arena-view")
arena_flash <- detect_flash(input=arena_view_fmf,
                            type="arena",
                            output=output_prefix,
                            flash_thresh=av_flash_thresh,
                            reuse=reuse)

## Find the other channel in fluo-view
fl1ref <- dipr::readTIFF2(fluo_view_tif_ch1, frames=fluo_flash$flflashes[1])
fl1ref <- EBImage::normalize(fl1ref)
fl2ref <- dipr::readTIFF2(fluo_view_tif, frames=fluo_flash$flflashes[1])
fl2ref <- EBImage::normalize(fl2ref)
fl2refcrop <- fl2ref[1025:2048,1:256] # Split the original image into two halves

center <- align_cameras(source=fl2refcrop,
                        template=fl1ref,
                        output=output_prefix,
                        center=c(0, 0),
                        zoom=1,
                        autopos=T)

imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height

# Channel detection


## Part 2. Syncing frames and generate frame IDs
syncing <- sync_frames(dir=dir,
                       fluo_flash=fluo_flash,
                       fly_flash=fly_flash,
                       arena_flash=arena_flash,
                       output=output_prefix,
                       reuse=reuse)

## Part 3. Analyze trajectories
trj_res <- analyze_trajectories(dir=dir,
                                output=output_prefix,
                                fpsfv=syncing$fpsfv,
                                interaction=interaction)

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
message(sprintf("Reading %s", fluo_view_tif_ch1))

# Analyze only part of the movie?
if(FOI!=F && length(FOI)==2){
  flimg1 <- dipr::readTIFF2(fluo_view_tif_ch1, start=FOI[1], end=FOI[2])
  flimg2 <- dipr::readTIFF2(fluo_view_tif_ch2, start=FOI[1], end=FOI[2])
  message(sprintf("Fluo-view frames from %d to %d will be analyzed.", FOI[1], FOI[2]))
  frid <- syncing$frid[FOI[1]:FOI[2]]
  frida <- syncing$frida[FOI[1]:FOI[2]]
}else{
  message("All frames will be analyzed.")
  frid <- syncing$frid
  frida <- syncing$frida
}

# Green or red channel?
flimg1int <- colMeans(flimg1, dim=2)
if(flimg1int[1] < mean(flimg1int)){
  green <- flimg1[,,seq(2, dim(flimg1)[3], 2)]
  red <- flimg2[,,seq(1, dim(flimg2)[3], 2)]
} else {
  green <- flimg1[,,seq(1, dim(flimg1)[3], 2)]
  red <- flimg2[,,seq(2, dim(flimg2)[3], 2)]
}

# Calculate dF/F


# Crop fluo-view movie for speed
if(dim(flimg1)[1] > 130){
  flimg1 <- flimg1[round((dim(flimg1)[1] - 128)/2):(round((dim(flimg1)[1]/2+128/2))-1),
                 round((dim(flimg1)[2] - 128)/2):(round((dim(flimg1)[2]/2+128/2))-1),]
  flimg2 <- flimg2[round((dim(flimg2)[1] - 128)/2):(round((dim(flimg2)[1]/2+128/2))-1),
                   round((dim(flimg2)[2] - 128)/2):(round((dim(flimg2)[2]/2+128/2))-1),]
  
}
flimg1rt <- EBImage::rotate(EBImage::flip(flimg1), rotate_camera)
flimg2rt <- EBImage::rotate(EBImage::flip(flimg2), rotate_camera)

# Load fly-view camera images
fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)

# Load arena-view camera images
avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
EBImage::writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
rm(avimgl)

## Part 7. Detect window on the head
fvimgbwbrfh <- detect_window(fvimgl=fvimgl, output=output_prefix, reuse=reuse)

## Part 8. Position calibration
fl1ref <- dipr::readTIFF2(fluo_view_tif_ch1, frames=fluo_flash$flflashes[1])
fl1ref <- EBImage::normalize(fl1ref)
fl2ref <- dipr::readTIFF2(fluo_view_tif_ch2, frames=fluo_flash$flflashes[1])
fl2ref <- EBImage::normalize(fl2ref)


center <- align_cameras(flref=fl1ref,
                        fvref=fl2ref,
                        output=output_prefix,
                        center=c(0, 0),
                        zoom=1,
                        autopos=T)

## Part 9. Image registration
regresi <- list()
for(rg in 1:dim(fvimgl)[3]){
  regresi[[rg]] <- niftyreg(img[,,rg], img[,,1],
                            scope="rigid", symmetric=F)
}

test <- niftyreg(fvimgl, fvimgl[,,1],
         scope="rigid", symmetric=F, sequentialInit=T)

regimgi <- test
regimgi <- array(sapply(regresi, function(x) x$image), dim=dim(red))
regimgi[which(is.na(regimgi)==T)] <- 0
display(normalize(regimgi))

registered_images <- register_images(fvimgl=fvimgl,
                                     flimgrt=flimgrt,
                                     fvimgbwbrfh=fvimgbwbrfh,
                                     angles=trj_res$angles,
                                     zoom=zoom,
                                     center=center,
                                     output=output_prefix,
                                     reuse=reuse)

## Part 10. Look for good frames based on size, position, motion, and focus
goodfr <- find_goodframes(window_mask=fvimgbwbrfh,
                          fvimgl=fvimgl,
                          output=output_prefix,
                          reuse=reuse)
