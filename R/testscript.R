#install.packages(c("devtools", "ggplot2", "RNiftyReg"))
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("BiocInstaller", "EBImage"))
#library(devtools)
#devtools::install_github("tkatsuki/FlyceptionR")
library(FlyceptionR)
library(zoo)
source("~/Flyception2R/R/align_cameras.R")
source("~/Flyception2R/R/imageJ_crop_append.R")
source("~/Flyception2R/R/sync_frames.R")

#dir <- "H:/P1_GCaMP6s_tdTomato_02202018/P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_4/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_7/"
dir <- "/Users/takeokatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_5/"
prefix <- "P1-Gal4_UAS-GCaMP6s_tdTomato_5"       # Will be used as a filename prefix
autopos <- T             # True if you want to align cameras automatically 
reuse <- F               # True if you want to reuse intermediate RDS files
fmf2tif <- T             # True if you want to convert fmf 
zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
FOI <- c(2438, 2793)                 # A vector specifying start and end frame (e.g. c(10,1000)). False if you want to analyze all frames.
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
#imageJ_crop_append(dir, ch=1, roi=ROI) # x and y coordinates of the top left corner, width, height
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
                        output=paste0(output_prefix, "_fl2fl1"),
                        center=c(0, 0),
                        zoom=1,
                        autopos=T)

#imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))


# Load fly-view camera images
fvref <- dipr::readFMF(fly_view_fmf, frames=c(fly_flash$fvflashes[1] + 1))[,,1]

# Align fly-view and fluo-view
center2 <- align_cameras(source=fvref/255,
                         template=flip(fl1ref),
                         output=paste0(output_prefix, "_fvfl1"),
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
flimgtest <- dipr::readTIFF2(fluo_view_tif_ch2, start=1, end=1)
fvimgltest <- dipr::readFMF(fly_view_fmf, frames=syncing$frid[1])

display(fvimgltest/255)
display(normalize(flimgtest)*5)

flimgtest2 <- EBImage::flip(dipr::readTIFF2(fluo_view_tif_ch2, start=517, end=527))
fvimgltest2 <- dipr::readFMF(fly_view_fmf, frames=syncing$frid[517:527])

display(fvimgltest2/255)
display(normalize(flimgtest2)*5)

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

# # Green or red channel?
# flimg1int <- colMeans(flimg1, dim=2)
# if(flimg1int[1] < mean(flimg1int)){
#   greenfr <- seq(2, dim(flimg1)[3], 2)
#   redfr <- seq(1, dim(flimg2)[3], 2)
# } else {
#   greenfr <- seq(1, dim(flimg1)[3], 2)
#   redfr <- seq(2, dim(flimg2)[3], 2)
# }
# 
# green <- flimg1[,,greenfr]
# red <- flimg2[,,redfr]
# 
# EBImage::writeImage(green/1000, file=paste0(dir, prefix, "green.tif"))
# EBImage::writeImage(red/1000, file=paste0(dir, prefix, "red.tif"))


# Load fly-view camera images
#fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)
fvimgl <- dipr::readFMF(fly_view_fmf, frames=seq(frid[1], frid[length(frid)]+20, by=1))

# Apply resize and translation to align with fluo-view
#fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear", output.dim=dim(red)[1:2]), (-center2 - 10))
fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear"), c(0, -4))
fvimgl <- fvimgl[11:250,11:250,1:dim(fvimgl)[3]]
#EBImage::writeImage(fvimgl/255, file=paste0(dir, prefix, "_fvimgl.tif"))

# Load arena-view camera images
#avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
#EBImage::writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
#rm(avimgl)

# Calculate head angles
fvimglbl <- gblur(fvimgl/255, 2)
fvimglbw <- thresh(fvimglbl, w=20, h=20, offset=0.2)
rm(fvimglbl)
centermask <- drawCircle(matrix(0,dim(fvimglbw)[1],dim(fvimglbw)[2]), dim(fvimglbw)[1]/2,
                         dim(fvimglbw)[2]/2, dim(fvimglbw)[1]*2/5, col=1, fill=1)
fvimglbw <- ssweep(fvimglbw, centermask, op="*")
fvimglbwseg <- bwlabel(fvimglbw)

ftrs <- list()
for (i in 1:dim(fvimglbwseg)[3]){
  ftrs[[i]] <- computeFeatures.moment(fvimglbwseg[,,i])
}

ang <- c()
centroid <- array(0, dim=c(dim(fvimglbwseg)[3],2))
markernum <- c()
for (im in 1:dim(fvimglbwseg)[3]){
  m <- ftrs[[im]]
  markernum[im] <- nrow(m)
  
  if(markernum[im]==3){
    distmat <- dist(m[1:3,1:2])
    maxpair <- which(distmat == max(distmat))
    centroid[im,] <- colMeans(m[,1:2])
    
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
    ang[im] <- angle 
  }
  
  if(markernum[im]!=3){
    print(paste0("Less or more than 3 beads detected in the ", im, "th frame."))
    centroid[im,] <- centroid[im-1,]
    ang[im] <- ang[im-1]
  }
}

angdiff <- c(0, diff(ang))
angsum <- zoo::rollsumr(angdiff, 20)
ang_thresh <- 0.02
goodangfr <- which(angsum < ang_thresh & angsum > -ang_thresh)
fvimglfr20 <- seq(1, dim(fvimgl)[3], by=20)
goodangfr20 <- which(fvimglfr20 %in% goodangfr)
goodmarkerfr <- which(markernum == 3) 
goodmarkerfr20 <- which(fvimglfr20 %in% goodmarkerfr)

objdist <- sqrt((centroid[,1]-dim(fvimgl)[1]/2)^2 + (centroid[,2]-dim(fvimgl)[2]/2)^2)
motion <- c(0, sqrt(diff(centroid[,1])^2 + diff(centroid[,2])^2))
motionsum <- zoo::rollsumr(motion, 20)
motion_thresh <- 10
goodmotionfr <- which(motionsum < motion_thresh)
goodmotionfr20 <- which(fvimglfr20 %in% goodmotionfr)

LoGkern <- round(dipr::LoG(9,9,1.4)*428.5)
flimg2log <- EBImage::filter2(flimg2, LoGkern)
centermask <- EBImage::drawCircle(flimg2[,,1]*0, dim(flimg2)[1]/2, dim(flimg2)[2]/2, 100, col=1, fill=T)
flimg2cntlog <- dipr::ssweep(flimg2log, centermask, op="*")
quantcnt <- apply(flimg2cntlog, 3, function(x) quantile(x, 0.9))
goodfocusfr <- which(quantcnt > 1100)
goodfr20 <-  Reduce(intersect, list(goodmarkerfr20, goodmotionfr20, goodangfr20, goodfocusfr))
goodfr20 <-  Reduce(intersect, list(goodmarkerfr20, goodmotionfr20, goodangfr20))


# Apply rotation compensation
rot <- fvimgl[,,fvimglfr20[goodfr20]]
for (r in 1:dim(rot)[3]){
  rot[,,r] <- RNiftyReg::rotate(fvimgl[,,fvimglfr20[goodfr20[r]]], ang[fvimglfr20[goodfr20[r]]], anchor = c("center"))
  
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
rottrans <- fvimgl[,,fvimglfr20[goodfr20]]
for (tr in 1:dim(rottrans)[3]){
  rottrans[,,tr] <- EBImage::translate(rot[,,tr], -centers[tr,])
}

EBImage::writeImage(rottrans/255, file=paste0(dir, prefix, "_rottrans.tif"))
display(normalize(rottrans))

## Apply transformation functions to fluo-view images
redrot <- flimg2[,,goodfr20]
for (rr in 1:dim(redrot)[3]){
  redrot[,,rr] <- as.Image(RNiftyReg::rotate(flimg2[,,goodfr20[rr]], ang[fvimglfr20[goodfr20[rr]]], anchor = c("center")))
}
greenrot <- flimg1[,,goodfr20]
for (rg in 1:dim(greenrot)[3]){
  greenrot[,,rg] <- as.Image(RNiftyReg::rotate(flimg1[,,goodfr20[rg]], ang[fvimglfr20[goodfr20[rg]]], anchor = c("center")))
}

display(normalize(redrot))
display(normalize(greenrot))

# Run image registration using the initial angles
# centermask <- drawCircle(matrix(0,dim(redrot)[1],dim(redrot)[2]), dim(redrot)[1]/2,
#                          dim(redrot)[2]/2, dim(redrot)[1]/5*2, col=1, fill=1)
# regresi <- list()
# for(rg in 1:dim(redrot)[3]){
#     regresi[[rg]] <- niftyreg(redrot[,,rg], redrot[,,1], scope="rigid", targetMask = centermask, symmetric=F)
# }
# 
# regimgi <- array(sapply(regresi, function(x) x$image), dim=dim(redrot))
# regimgi[which(is.na(regimgi)==T)] <- 0
# display(normalize(regimgi))


centerr <- array(0, dim=c(dim(redrot)[3],2))
for (cr in 1:dim(redrot)[3]){
  centerr[cr,] <- align_cameras(source=redrot[,,cr],
                                template=redrot[,,1],
                                output=output_prefix,
                                center=c(0, 0),
                                zoom=1,
                                autopos=T)
}

redrottrans <- redrot
for (trr in 1:dim(redrottrans)[3]){
  redrottrans[,,trr] <- EBImage::translate(redrot[,,trr], -centers[trr,])
}
display(normalize(redrottrans))

greenrottrans <- greenrot
for (trg in 1:dim(greenrottrans)[3]){
  greenrottrans[,,trg] <- EBImage::translate(greenrot[,,trg], -centers[trg,])
}

# Segment neurons
redwindow <- redrottrans[(dim(redrottrans)[1]/2-35):(dim(redrottrans)[1]/2+35),
                         (dim(redrottrans)[2]/2-20):(dim(redrottrans)[2]/2+15),]
greenwindow <- greenrottrans[(dim(greenrottrans)[1]/2-35):(dim(greenrottrans)[1]/2+35),
                             (dim(greenrottrans)[2]/2-20):(dim(greenrottrans)[2]/2+15),]

display(normalize(redwindow))
redwindowmed <- EBImage::medianFilter(redwindow/2^16, size=2)
greenwindowmed <- EBImage::medianFilter(greenwindow/2^16, size=2)
display(normalize(redwindowmed))
redwindowmedth <- thresh(redwindowmed, w=10, h=10, offset=0.0003)
display(redwindowmedth)

# F_ratio image
redmasked <- redwindowmed*redwindowmedth
greenmasked <- greenwindowmed*redwindowmedth
greenperred <- greenmasked/redmasked
greenperredave <- colMeans(greenperred, dim=2, na.rm=T)
plot(greenperredave)
greenperred[which(is.na(greenperred)==T)] <- 0
grratiocolor <- array(0, dim=c(dim(greenperred)[c(1,2)], 3, dim(greenperred)[3]))
for(cfr in 1:dim(greenperred)[3]){
  grratiocolor[,,,cfr] <- dipr::pseudoColor(greenperred[,,cfr], 160, 250)
}
grratiocolor <- Image(grratiocolor, colormode="Color")
display(grratiocolor)

# Overlay fly_view and F_ratio image
rottransmask <- array(0, dim=c(dim(rottrans)[c(1,2)], dim(rottrans)[3]))
rottransmask[(dim(rottrans)[1]/2-35):(dim(rottrans)[1]/2+35),
             (dim(rottrans)[2]/2-20):(dim(rottrans)[2]/2+15),] <- redwindowmedth

rottranscolor <- array(0, dim=c(dim(rottrans)[c(1,2)], 3, dim(rottrans)[3]))
rottranscolor[,,1,] <- rottrans/255*(1-rottransmask)
rottranscolor[,,2,] <- rottrans/255*(1-rottransmask)
rottranscolor[,,3,] <- rottrans/255*(1-rottransmask)

grratiocolorl <- rottranscolor*0
grratiocolorl[(dim(grratiocolorl)[1]/2-35):(dim(grratiocolorl)[1]/2+35),
              (dim(grratiocolorl)[2]/2-20):(dim(grratiocolorl)[2]/2+15),,] <- grratiocolor
flyviewcolor <- rottranscolor + grratiocolorl
flyviewcolor <- Image(flyviewcolor, colormode="Color")
display(flyviewcolor)

# overlay green channel and F_ratio color image
greenrottranscol <- array(0, dim=c(dim(greenrottrans)[c(1,2)], 3, dim(greenrottrans)[3]))
greenrottranscol[,,1,] <- greenrottrans/2^16*(1-rottransmask)
greenrottranscol[,,2,] <- greenrottrans/2^16*(1-rottransmask)
greenrottranscol[,,3,] <- greenrottrans/2^16*(1-rottransmask)
greenrottranscol <- normalize(greenrottranscol, separate=F)
greencolor <- greenrottranscol + grratiocolorl
greencolor <- Image(greencolor, colormode="Color")
display(greencolor)

# Create side-by-side view of fly_view and fluo_view images
frgcombined <- array(dim=c(dim(rottrans)[1]*4, dim(rottrans)[2], 3, dim(rottrans)[3]))
frgcombined[1:240,1:240,1,1:dim(rottrans)[3]] <- normalize(rottrans, separate=F)
frgcombined[1:240,1:240,2,1:dim(rottrans)[3]] <- normalize(rottrans, separate=F)
frgcombined[1:240,1:240,3,1:dim(rottrans)[3]] <- normalize(rottrans, separate=F)
frgcombined[241:480,1:240,1,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F)
frgcombined[241:480,1:240,2,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F)
frgcombined[241:480,1:240,3,1:dim(redrottrans)[3]] <- normalize(redrottrans, separate=F)
frgcombined[481:720,1:240,1,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F)
frgcombined[481:720,1:240,2,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F)
frgcombined[481:720,1:240,3,1:dim(greenrottrans)[3]] <- normalize(greenrottrans, separate=F)
frgcombined[721:960,1:240,,1:dim(greenrottrans)[3]] <- greencolor
frgcombined <-  Image(frgcombined, colormode="Color")

display(frgcombined)
EBImage::writeImage(redrottrans/2^16, file=paste0(dir, prefix, "_redrottrans.tif"))
EBImage::writeImage(greenrottrans/2^16, file=paste0(dir, prefix, "_greenrottrans.tif"))
EBImage::writeImage(frgcombined, file=paste0(dir, prefix, "_frgcombined_goodfr20_normalized.tif"))
