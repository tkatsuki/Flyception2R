install.packages(c("devtools", "ggplot2", "RNiftyReg"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocInstaller", "EBImage"))
library(devtools)
devtools::install_github("tkatsuki/FlyceptionR")
library(FlyceptionR)

#dir <- "H:/P1_GCaMP6s_tdTomato_02202018/P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy/"  # Don't forget the slash at the end
dir <- "C:/Users/tkatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_4/"  # Don't forget the slash at the end
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
fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))

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
  flimg1 <- flip(flimg1) # flip images to match fly-view
  flimg2 <- flip(flimg2) # flip images to match fly-view
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

# Load fly-view camera images
fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)

# Load arena-view camera images
avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
EBImage::writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
rm(avimgl)

# Calculate head angles
fvimglbl <- gblur(fvimgl/255, 2)
fvimglbw <- thresh(fvimglbl, w=20, h=20, offset=0.1)
display(fvimglbw)
fvimglbwseg <- bwlabel(fvimglbw)

ang <- c()

#for (i in 1:dim(fvimglbwseg)[3]){
for (i in 1:100){
  m <- computeFeatures.moment(fvimglbwseg[,,i])
  distmat <- dist(m[1:3,1:2])
  maxpair <- which(distmat == max(distmat))
  
  if(maxpair == 1){ # pair 1-2
    angle <- atan((m[2,2] - m[1,2])/(m[2,1] - m[1,1]))
    if (angle < 0){
      dleft <- (m[1,1] - 1 - m[1,1])*(m[2,2] - m[1,2]) - (m[1,2] - m[1,2])*(m[2,1] - m[1,1])
      d <- (m[3,1] - m[1,1])*(m[2,2] - m[1,2]) - (m[3,2] - m[1,2])*(m[2,1] - m[1,1])
      if (dleft * d > 0){ # Triangle facing left
        angle <- angle + pi
      }
    }
    if (angle > 0){
      dleft <- (m[1,1] - 1 - m[1,1])*(m[2,2] - m[1,2]) - (m[1,2] - m[1,2])*(m[2,1] - m[1,1])
      d <- (m[3,1] - m[1,1])*(m[2,2] - m[1,2]) - (m[3,2] - m[1,2])*(m[2,1] - m[1,1])
      if (dleft * d < 0){ # Triangle facing left
        angle <- angle + pi
      }
    }
  }
  
  if(maxpair == 2){ # pair 1-3
    angle <- atan((m[3,2] - m[1,2])/(m[3,1] - m[1,1]))
    if (angle < 0){
      dleft <- (m[1,1] - 1 - m[1,1])*(m[3,2] - m[1,2]) - (m[1,2] - m[1,2])*(m[3,1] - m[1,1])
      d <- (m[2,1] - m[1,1])*(m[3,2] - m[1,2]) - (m[2,2] - m[1,2])*(m[3,1] - m[1,1])
      if (dleft * d > 0){ # Triangle facing left
        angle <- angle + pi
      }
    }
    if (angle > 0){
      dleft <- (m[1,1] - 1 - m[1,1])*(m[3,2] - m[1,2]) - (m[1,2] - m[1,2])*(m[3,1] - m[1,1])
      d <- (m[2,1] - m[1,1])*(m[3,2] - m[1,2]) - (m[2,2] - m[1,2])*(m[3,1] - m[1,1])
      if (dleft * d < 0){ # Triangle facing left
        angle <- angle + pi
      }
    }
  }
  
  if(maxpair == 3){ # pair 2-3
    angle <- atan((m[2,2] - m[3,2])/(m[2,1] - m[3,1]))
    if (angle < 0){
      dleft <- (m[3,1] - 1 - m[3,1])*(m[2,2] - m[3,2]) - (m[3,2] - m[3,2])*(m[2,1] - m[3,1])
      d <- (m[1,1] - m[3,1])*(m[2,2] - m[3,2]) - (m[1,2] - m[3,2])*(m[2,1] - m[3,1])
      if (dleft * d > 0){ # Triangle facing left
        angle <- angle + pi
      }
    }
    if (angle > 0){
      dleft <- (m[3,1] - 1 - m[3,1])*(m[2,2] - m[3,2]) - (m[3,2] - m[3,2])*(m[2,1] - m[3,1])
      d <- (m[1,1] - m[3,1])*(m[2,2] - m[3,2]) - (m[1,2] - m[3,2])*(m[2,1] - m[3,1])
      if (dleft * d < 0){ # Triangle facing left
        angle <- angle + pi
      }
    }
  }
  ang[i] <- angle
}

# Build affine matrix for rotation
aff <- list()
#for(a in 1:dim(img)[3]){
for(a in 1:100){
    aff[[a]] <- buildAffine(angles=c(0,0, ang[a]), source=fvimgl[,,1], anchor="center")
}

# Apply rotation compensation
rot <- fvimgl
#for (r in 1:dim(img)[3]){
  for (r in 1:100){
  rot[,,r] <- as.Image(rotate(fvimgl[,,r], ang[r], anchor = c("center")))
  }

# Template matching
centers <- array(0, dim=c(100,2))

for (c in 1:100){
  centers[c,] <- align_cameras(source=rot[,,c],
                           template=rot[,,1],
                           output=output_prefix,
                           center=c(0, 0),
                           zoom=1,
                           autopos=T)
}

# Apply translation compensation
rottrans <- fvimgl
for (tr in 1:100){
  rottrans[,,tr] <- EBImage::translate(rot[,,tr], -centers[tr,])
}

EBImage::writeImage(rot/255, file=paste0(dir, prefix, "_rot.tif"))

# Run image registration using the initial angles
regresi <- list()
if(cores==1){
  for(rg in 1:dim(img)[3]){
    regresi[[rg]] <- niftyreg(rotc[,,rg], rotc[,,1],
                              scope="rigid", symmetric=F)
  }
}else{
  regresi <- foreach(rg = 1:dim(fvimgli)[3]) %dopar% niftyreg(fvimgli[,,rg], fvimgrt1sti, init=aff[[rg]], scope="rigid", symmetric=F, internal=FALSE)
}

## Part 9. Image registration
regresi <- list()
for(rg in 1:dim(fvimgl)[3]){
  regresi[[rg]] <- niftyreg(img[,,rg], img[,,1],
                            scope="rigid", symmetric=F)
}

test <- niftyreg(rot, rot[,,1],
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
