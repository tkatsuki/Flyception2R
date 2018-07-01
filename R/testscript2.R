#install.packages(c("devtools", "ggplot2", "RNiftyReg"))
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("BiocInstaller", "EBImage"))
#library(devtools)
#devtools::install_github("tkatsuki/FlyceptionR")
library(FlyceptionR)
library(zoo)
library(magick)
library(loggit)
source("~/Flyception2R/R/align_cameras.R")
source("~/Flyception2R/R/imageJ_crop_append.R")
source("~/Flyception2R/R/sync_frames.R")
source("~/Flyception2R/R/detect_flash.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/align_cameras.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/imageJ_crop_append.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/sync_frames.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/detect_flash.R")

#To do
# Should stop when the number of flash detected do not match between flou-view and fly-view
# Save flash data

#dir <- "H:/P1_GCaMP6s_tdTomato_02202018/P1-Gal4_UAS-GCaMP6s_tdTomato_4Copy/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_4/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1/"  # Don't forget the slash at the end
#dir <- "C:/Users/tkatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_7/"
#dir <- "/Users/takeokatsuki/Desktop/P1-Gal4_UAS-GCaMP6s_tdTomato_5/"
dir <- "/Volumes/LaCie/P1_GCaMP6s_tdTomato_06212018_CW_Dual_Laser/P1-Gal4_UAS-GCaMP6s_tdTomato_7/"
#dir <- "/Volumes/LaCie/P1_GCaMP6s_tdTomato_06212018_CW_Dual_Laser/P1-Gal4_UAS-GCaMP6s_tdTomato_6/"
prefix <- paste0("P1-Gal4_UAS-GCaMP6s_tdTomato_7")       # Will be used as a filename prefix
autopos <- T             # True if you want to align cameras automatically 
reuse <- T               # True if you want to reuse intermediate RDS files
fmf2tif <- T             # True if you want to convert fmf 
zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
FOI <-  c(4200, 4460)                 # A vector specifying start and end frame (e.g. c(10,1000)). False if you want to analyze all frames.
ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
binning <- 1             # Binning of the fluo-view camera
fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
interaction <- T         # True if you want to analyze fly-fly interaction
dist_thresh <- 4         # Threshold for detecting fly-fly interaction based on distance
rotate_camera <- -180    # Rotation angle needed to align fluo-view and fly-view
window_size <- c(68, 28) # Size of a rectangle window on the head for segmentation. Choose even numbers.
window_offset <- c(-4, 25)     # Offset of the window from the center of the image. Positive x moves right
outdir <- paste0(dir, paste0(FOI, collapse="_"), "/")

dir.create(outdir)
output_prefix <- paste0(outdir, prefix)

loggit::setLogFile(paste0(output_prefix, "_log.json"))
loggit::message(outdir)

fluo_view_tif <- paste0(dir, list.files(dir, pattern="Pos0\\.ome\\.tif$"))
fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))

## First crop the fluo-view image with one ROI and find the flash frame

if(length(list.files(dir, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))==0){
  imageJ_crop_append(dir, ch=1, roi=ROI) # x and y coordinates of the top left corner, width, height
}
fluo_view_tif_ch1 <- paste0(dir, list.files(dir, pattern="ome\\.ch1\\.crop\\.concat\\.tif$"))
flnframe <- dipr::readTIFF2(fluo_view_tif_ch1, getFrames = T)
  
message("Detecting flash in fluo-view")

# if(file.exists(paste0(output_prefix, "_fluo_flash.RDS"))==T &
#    file.exists(paste0(output_prefix, "_fly_flash.RDS"))==T &
#    file.exists(paste0(output_prefix, "_arena_flash.RDS"))==T &
#    reuse==T){
#   message("Loading RDS file")
#   fluo_flash <- readRDS(paste0(output_prefix, "_fluo_flash.RDS"))
#   fly_flash <- readRDS(paste0(output_prefix, "_fly_flash.RDS"))
#   arena_flash <- readRDS(paste0(output_prefix, "_arena_flash.RDS"))
#   
# } elase {
  
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
  # if(saveRDS==T){
  #   saveRDS(regimgi, paste0(output, "_regimgi.RDS"))
  #   saveRDS(regresi, paste0(output, "_regresi.RDS"))
  # }




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
                        autopos=T,
                        ROI=c(1, 1, 450, 50))


if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
  imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
}
fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))


# Load fly-view camera images
fvref <- dipr::readFMF(fly_view_fmf, frames=c(fly_flash$fvflashes[1] + 1))[,,1]

# Align fly-view and fluo-view
center2 <- align_cameras(source=fvref/255,
                         template=flip(fl1ref),
                         output=paste0(output_prefix, "_fvfl1"),
                         center=c(0, 0),
                         zoom=1.085,
                         autopos=T,
                         ROI=F)

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
  FOI <- c(1, flnframe)
}

# Load fly-view camera images
fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)

# Apply resize and translation to align with fluo-view
fvimgl <- EBImage::translate(EBImage::resize(fvimgl, dim(fvimgl)[1]*1.085, filter="bilinear"), -center2)
fvimgl <- fvimgl[11:250,11:250,1:dim(fvimgl)[3]]
EBImage::writeImage(fvimgl/255, file=paste0(output_prefix, "_fvimgl.tif"))

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
ang_thresh <- 0.02
goodangfr <- which(angdiff < ang_thresh & angdiff > -ang_thresh)
goodmarkerfr <- which(markernum == 3) 

objdist <- sqrt((centroid[,1]-dim(fvimgl)[1]/2)^2 + (centroid[,2]-dim(fvimgl)[2]/2)^2)
motion <- c(0, sqrt(diff(centroid[,1])^2 + diff(centroid[,2])^2))
motion_thresh <- 2
goodmotionfr <- which(motion < motion_thresh)

LoGkern <- round(dipr::LoG(9,9,1.4)*428.5)
flimg2log <- EBImage::filter2(flimg2, LoGkern)
centermask <- EBImage::drawCircle(flimg2[,,1]*0, dim(flimg2)[1]/2, dim(flimg2)[2]/2, 100, col=1, fill=T)
flimg2cntlog <- dipr::ssweep(flimg2log, centermask, op="*")
quantcnt <- apply(flimg2cntlog, 3, function(x) quantile(x, 0.9))
goodfocusfr <- which(quantcnt > 1000)
goodfr <- Reduce(intersect, list(goodmarkerfr, goodmotionfr, goodangfr, goodfocusfr))

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
display(normalize(rottrans))

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

# Segment neurons
redwindow <- redrottrans[(dim(redrottrans)[1]/2 + window_offset[1] - window_size[1]/2):
                           (dim(redrottrans)[1]/2 + window_offset[1] + window_size[1]/2),
                         (dim(redrottrans)[2]/2 + window_offset[2] - window_size[2]/2):
                           (dim(redrottrans)[2]/2 + window_offset[2] + window_size[2]/2),]
greenwindow <- greenrottrans[(dim(greenrottrans)[1]/2 + window_offset[1] - window_size[1]/2):
                             (dim(greenrottrans)[1]/2 + window_offset[1] + window_size[1]/2),
                           (dim(greenrottrans)[2]/2 + window_offset[2] - window_size[2]/2):
                             (dim(greenrottrans)[2]/2 + window_offset[2] + window_size[2]/2),]

display(normalize(redwindow))
redwindowmed <- EBImage::medianFilter(redwindow/2^16, size=2)
greenwindowmed <- EBImage::medianFilter(greenwindow/2^16, size=2)
display(normalize(redwindowmed))
redwindowmedth <- EBImage::thresh(redwindowmed, w=10, h=10, offset=0.0003)
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
  grratiocolor[,,,cfr] <- dipr::pseudoColor(greenperred[,,cfr], 180, 220)
}
grratiocolor <- Image(grratiocolor, colormode="Color")
display(grratiocolor)

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
display(flyviewcolor)
EBImage::writeImage(flyviewcolor, file=paste0(output_prefix, "_flyviewcolor.tif"))

# overlay red channel and F_ratio color image
redrottranscol <- array(0, dim=c(dim(redrottrans)[c(1,2)], 3, dim(redrottrans)[3]))
redrottranscol[,,1,] <- redrottrans*(1-rottransmask)
redrottranscol[,,2,] <- redrottrans*(1-rottransmask)
redrottranscol[,,3,] <- redrottrans*(1-rottransmask)
redrottranscol <- normalize(redrottranscol, separate=F, inputRange=c(180, 400))
redcolor <- redrottranscol + grratiocolorl
redcolor <- Image(redcolor, colormode="Color")
display(redcolor)

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

display(frgcombined)
EBImage::writeImage(normalize(redrottrans, separate=F, inputRange=c(180, 400)), file=paste0(output_prefix, "_redrottrans.tif"))
EBImage::writeImage(normalize(greenrottrans, separate=F, inputRange=c(180, 300)), file=paste0(output_prefix, "_greenrottrans.tif"))
EBImage::writeImage(frgcombined, file=paste0(output_prefix, "_frgcombined_goodfr20_normalized.tif"))

# Calculate dF/F
intensity <- zoo::rollmean(greenperredave, 3, align="left")
datint <- data.frame(x=goodfr[1:(length(goodfr)-2)], y=intensity)
png(file=paste0(output_prefix, "_datint.png"), width=400, height=400)
plot(datint)  
dev.off()

F0int <- intensity[1]
deltaFint <- intensity - F0int
dFF0int <- deltaFint/F0int * 100
datdFF0 <- data.frame(x=goodfr[1:(length(goodfr)-2)], y=dFF0int)
png(file=paste0(output_prefix, "_datdFF0.png"), width=400, height=400)
plot(datdFF0)
dev.off()

loggit::message(sprintf("window_size was x=%d y=%d", window_size[1], window_size[2]))
loggit::message(sprintf("window_offset was x=%d y=%d", window_offset[1], window_offset[2]))
loggit::message(sprintf("FOI was from %d to %d",  FOI[1], FOI[2])) 
loggit::message(paste0("Max F_ratio intensity in this bout was ", max(intensity)))
loggit::message(paste0("Number of good frames was ", length(goodfr)))

loggit::message(sprintf("||c(%d, %d) ||c(%d, %d) ||c(%d, %d) ||%d ||%.3f ||", 
                        FOI[1], FOI[2], window_size[1], window_size[2], window_offset[1], window_offset[2], length(goodfr), max(intensity)))

