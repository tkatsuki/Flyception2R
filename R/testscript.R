install.packages(c("devtools", "ggplot2", "RNiftyReg"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocInstaller", "EBImage"))
library(devtools)
devtools::install_github("tkatsuki/FlyceptionR")
library(FlyceptionR)
library(zoo)

#dir <- "H:/P1_GCaMP6s_tdTomato_02202018/P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_4/"  # Don't forget the slash at the end
dir <- "C:/Users/tkatsuki/Desktop/P1/"  # Don't forget the slash at the end
prefix <- "P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy"       # Will be used as a filename prefix
autopos <- T             # True if you want to align cameras automatically 
reuse <- F               # True if you want to reuse intermediate RDS files
fmf2tif <- T             # True if you want to convert fmf 
zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
FOI <- c(250, 350)                 # A vector specifying start and end frame (e.g. c(10,1000)). False if you want to analyze all frames.
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
fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))


## First crop the fluo-view image with one ROI and find the flash frame
imageJ_crop_append(dir, ch=1, roi=ROI) # x and y coordinates of the top left corner, width, height
fluo_view_tif_ch1 <- paste0(dir, list.files(dir, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))

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


# Load fly-view camera images
fvref <- dipr::readFMF(fly_view_fmf, frames=c(fly_flash$fvflashes[1] + 1))[,,1]

# Align fly-view and fluo-view
center2 <- align_cameras(source=fvref,
                        template=flip(fl1ref),
                        output=output_prefix,
                        center=c(0, 0),
                        zoom=1.085,
                        autopos=T)

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
  greenfr <- seq(2, dim(flimg1)[3], 2)
  redfr <- seq(1, dim(flimg2)[3], 2)
} else {
  greenfr <- seq(1, dim(flimg1)[3], 2)
  redfr <- seq(2, dim(flimg2)[3], 2)
}

green <- flimg1[,,greenfr]
red <- flimg2[,,redfr]

EBImage::writeImage(green/1000, file=paste0(dir, prefix, "green.tif"))
EBImage::writeImage(red/1000, file=paste0(dir, prefix, "red.tif"))


# Load fly-view camera images
#fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)
fvimgl <- dipr::readFMF(fly_view_fmf, frames=seq(frid[1], frid[length(frid)], by=1))

# Apply resize and translation to align with fluo-view
#fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear", output.dim=dim(red)[1:2]), (-center2 - 10))
fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear"), c(0, -4))
fvimgl <- fvimgl[11:250,11:250,1:dim(fvimgl)[3]]

# Load arena-view camera images
avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
EBImage::writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
rm(avimgl)

# Calculate head angles
fvimglbl <- gblur(fvimgl/255, 2)
fvimglbw <- thresh(fvimglbl, w=20, h=20, offset=0.1)
centermask <- drawCircle(matrix(0,dim(fvimglbl)[1],dim(fvimglbl)[2]), dim(fvimglbl)[1]/2,
                         dim(fvimglbl)[2]/2, dim(fvimglbl)[1]*2/5, col=1, fill=1)
fvimglbw <- ssweep(fvimglbw, centermask, op="*")
#display(fvimglbw)
fvimglbwseg <- bwlabel(fvimglbw)

ang <- c()
centroid <- array(0, dim=c(dim(fvimgl)[3],2))
  
for (i in 1:dim(fvimglbwseg)[3]){
  m <- computeFeatures.moment(fvimglbwseg[,,i])
  if(nrow(m)==3){
    distmat <- dist(m[1:3,1:2])
    maxpair <- which(distmat == max(distmat))
    centroid[i,] <- colMeans(m[,1:2])
      
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
  # if(){
  #   
  # }

}

# Apply rotation compensation
rot <- fvimgl
for (r in 1:dim(fvimgl)[3]){
#  for (r in 1:220){
  rot[,,r] <- as.Image(RNiftyReg::rotate(fvimgl[,,r], ang[r], anchor = c("center")))
  }

# Template matching
centers <- array(0, dim=c(dim(fvimgl)[3],2))

for (c in 1:dim(fvimgl)[3]){
  centers[c,] <- align_cameras(source=rot[,,c],
                           template=rot[,,1],
                           output=output_prefix,
                           center=c(0, 0),
                           zoom=1,
                           autopos=T)
}


objdist <- sqrt((centroid[,1]-dim(fvimgl)[1]/2)^2 + (centroid[,2]-dim(fvimgl)[2]/2)^2)
motion <- c(0, sqrt(diff(centroid[,1])^2 + diff(centroid[,2])^2))
motionsum <- zoo::rollsumr(motion, 20)
motion_thresh <- 20
goodmotionfr <- which(motionsum < motion_thresh)
fvimglfr20 <- seq(1, dim(fvimgl)[3], by=19)
goodmotionfr20 <- which(fvimglfr20 %in% goodmotionfr)
message(sprintf("The following frames have too large motion: %s", paste((1:length(motion))[-goodmotionfr], collapse=" ")))

matplot(centers, type="l")

# Apply translation compensation
rottrans <- fvimgl
for (tr in 1:dim(fvimgl)[3]){
  rottrans[,,tr] <- EBImage::translate(rot[,,tr], -centers[tr,])
}

EBImage::writeImage(rottrans/255, file=paste0(dir, prefix, "_rottrans.tif"))
display(normalize(rottrans))

## Apply transformation functions to fluo-view images
redrot <- flimg2
for (rr in 1:dim(fvimgl)[3]){
  #for (rr in 2:110){
  redrot[,,rr] <- as.Image(RNiftyReg::rotate(flimg2[,,rr], ang[rr], anchor = c("center")))
}
greenrot <- flimg1
for (rg in 1:dim(fvimgl)[3]){
  #for (rg in 2:110){
  greenrot[,,rg] <- as.Image(RNiftyReg::rotate(flimg1[,,rg], ang[rg], anchor = c("center")))
}
# 
# redrot <- red
# for (rr in 2:dim(fvimgl)[3]){
# #for (rr in 2:110){
#   redrot[,,rr] <- as.Image(RNiftyReg::rotate(red[,,rr], ang[redfr[rr-1]], anchor = c("center")))
# }
# greenrot <- green
# for (rg in 2:dim(fvimgl)[3]){
# #for (rg in 2:110){
#   greenrot[,,rg] <- as.Image(RNiftyReg::rotate(green[,,rg], ang[greenfr[rg-1]], anchor = c("center")))
# }

display(normalize(redrot))
display(normalize(greenrot))

centerr <- array(0, dim=c(dim(fvimgl)[3],2))
for (cr in 1:dim(fvimgl)[3]){
  centerr[cr,] <- align_cameras(source=redrot[,,cr],
                                template=redrot[,,1],
                                output=output_prefix,
                                center=c(0, 0),
                                zoom=1,
                                autopos=T)
}

redrottrans <- flimg2
for (trr in 1:dim(fvimgl)[3]){
  redrottrans[,,trr] <- EBImage::translate(redrot[,,trr], -centers[trr,])
}
display(normalize(redrottrans))

greenrottrans <- flimg1
for (trg in 1:dim(fvimgl)[3]){
  greenrottrans[,,trg] <- EBImage::translate(greenrot[,,trg], -centers[trg,])
}

frgcombined <- array(dim=c(dim(fvimgl)[1]*3, dim(fvimgl)[2], dim(fvimgl)[3]))
frgcombined[1:240,1:240,1:dim(fvimgl)[3]] <- rottrans
frgcombined[241:480,1:240,1:dim(fvimgl)[3]] <- greenrottrans
frgcombined[481:720,1:240,1:dim(fvimgl)[3]] <- redrottrans
display(normalize(frgcombined))
  
display(normalize(redrottrans))
display(normalize(greenrottrans))
EBImage::writeImage(redrottrans/2^16, file=paste0(dir, prefix, "_redrottrans.tif"))
EBImage::writeImage(greenrottrans/2^16, file=paste0(dir, prefix, "_greenrottrans.tif"))
EBImage::writeImage(frgcombined/2^16, file=paste0(dir, prefix, "_frgcombined.tif"))


centermask <- drawCircle(matrix(0,dim(ratioimg)[1],dim(ratioimg)[2]), dim(ratioimg)[1]/2,
                         dim(ratioimg)[2]/2, dim(ratioimg)[1]*2/5, col=1, fill=1)
redregres <- list()
for(redrg in 1:dim(ratioimg)[3]){
  redregres[[redrg]] <- niftyreg(ratioimg[,,redrg], ratioimg[,,1],
                                 scope="rigid", symmetric=F)
}

redregimg <- array(sapply(redregres, function(x) x$image), dim=dim(ratioimg))
redregimg[which(is.na(redregimg)==T)] <- 0
display(normalize(redregimg))
writeImage((255-regimgi)/255, file=paste0(output, "_regimgi.tif"))


## Part 10. Look for good frames based on size, position, motion, and focus
goodfr <- find_goodframes(window_mask=fvimgbwbrfh,
                          fvimgl=fvimgl,
                          output=output_prefix,
                          reuse=reuse)
