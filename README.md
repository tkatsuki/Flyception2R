# FlyceptionR
R scripts and utilities for analyzing Flyception data

## Requirement
Windows: R and Rtools

Mac: R and Xcode

ImageJ https://imagej.nih.gov/ij/

Set ImageJ (and Fiji if present) to save Tiff in Intel byte order from Edit -> Options -> Input/Output... 

## Installation
The following commands will automatically install packages necessary for running FlyceptionR.

```
install.packages(c("devtools", "ggplot2", "RNiftyReg"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocInstaller", "EBImage"))
library(devtools)
devtools::install_github("tkatsuki/Flyception2R")
library(Flyception2R)
```

## Usage example

### Part 0. Initialization

#### Set variables
```
dir <- "/example/data/"  # Don't forget the slash at the end
prefix <- "data_1"       # Will be used as a filename prefix
autopos <- T             # True if you want to align cameras automatically 
reuse <- F               # True if you want to reuse intermediate RDS files
fmf2tif <- T             # True if you want to convert fmf 
zoom <- 1.13             # Zoom ratio: fluo-view/fly-view
FOI <- F                 # A vector specifying start and end frame (e.g. c(10,1000)). False if you want to analyze all frames.
binning <- 1             # Binning of the fluo-view camera
fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
fv_flash_thresh <- 135   # Threshold for detecting flash in fly-view
av_flash_thresh <- 50    # Threshold for detecting flash in arena-view
interaction <- T         # True if you want to analyze fly-fly interaction
dist_thresh <- 4         # Threshold for detecting fly-fly interaction based on distance
rotate_camera <- -180    # Rotation angle needed to align fluo-view and fly-view
```

#### Start logging
```
rlogging::SetLogFile(base.file=paste0(prefix, "_log.txt"), folder=dir)
message(dir)
```

#### Prepare filenames
```
output_prefix <- paste0(dir, prefix)
fluo_view_tif <- paste0(dir, list.files(dir, pattern="ome\\.tif$"))
fly_view_fmf <- paste0(dir, list.files(dir, pattern="^fv.*fmf$"))
arena_view_fmf <- paste0(dir, list.files(dir, pattern="^av.*fmf$"))
```

### Part 1. Detect flash
Required step.

```
message("Detecting flash in fluo-view")
fluo_flash <- detect_flash(input=fluo_view_tif,
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
```

### Part 2. Synchronize frames and generate frame IDs
Required step. Identifies which frame of fly-view and arena-view videos correspond to the fluo-view frames.

```
syncing <- sync_frames(dir=dir,
                       fluo_flash=fluo_flash,
                       fly_flash=fly_flash,
                       arena_flash=arena_flash,
                       output=output_prefix,
                       reuse=reuse)
```

### Part 3. Analyze trajectories
Required step.

```
trj_res <- analyze_trajectories(dir=dir,
                                output=output_prefix,
                                fpsfv=syncing$fpsfv,
                                interaction=interaction)
```

### Part 4. Detect stimulus
Optional. Perform this step if you want to detect stimulus in the fv-traj file.

```
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
```

### Part 5. Detect interaction
Optional. Perform this step if you want to detect fly-fly interaction based on distance.

```
if(interaction==T){
  message("Detecting interaction")
  closefr <- which(trj_res$flydist < dist_thresh)
  closefrid <- sapply(closefr, function(x) which.min(abs(syncing$frida-x)))
  write.table(closefrid, paste0(dir, prefix, "_closefrid.txt"))
}
```

### Part 6. Load images
Required step.

```
message(sprintf("Reading %s", fluo_view_tif))

# Analyze only part of the movie?
if(FOI!=F && length(FOI)==2){
  flimg <- dipr::readTIFF2(fluo_view_tif, start=FOI[1], end=FOI[2])
  message(sprintf("Fluo-view frames from %d to %d will be analyzed.", FOI[1], FOI[2]))
  frid <- syncing$frid[FOI[1]:FOI[2]]
  frida <- syncing$frida[FOI[1]:FOI[2]]
}else{
  message("All frames will be analyzed.")
  frid <- syncing$frid
  frida <- syncing$frida
}
# Crop fluo-view movie for faster processing
if(dim(flimg)[1] > 130){
  flimg <- flimg[round((dim(flimg)[1] - 128)/2):(round((dim(flimg)[1]/2+128/2))-1),
                 round((dim(flimg)[2] - 128)/2):(round((dim(flimg)[2]/2+128/2))-1),]
}
flimgrt <- rotate(flip(flimg), rotate_camera)
# Load fly-view camera images
fvimgl <- dipr::readFMF(fly_view_fmf, frames=frid)
# Load arena-view camera images
avimgl <- dipr::readFMF(arena_view_fmf, frames=frida)
writeImage(avimgl/255, file=paste0(dir, prefix, "_avimgl_fr_", frida[1], "-", tail(frida, n=1), ".tif"))
rm(avimgl)
```


### Part 7. Detect window on the head
Required step.

```
fvimgbwbrfh <- detect_window(fvimgl=fvimgl, output=output_prefix, reuse=reuse)
```

### Part 8. Position calibration
Required step.

```
flref <- dipr::readTIFF2(fluo_view_tif, frames=fluo_flash$flflashes[1])
flref <- EBImage::normalize(EBImage::rotate(EBImage::flip(flref), rotate_camera))
fvref <- dipr::readFMF(filepath=fly_view_fmf,
                       frames=fly_flash$fvflashes[1])[,,1]/255
center <- align_cameras(flref=flref,
                        fvref=fvref,
                        output=output_prefix,
                        zoom=zoom,
                        autopos=autopos)
```

### Part 9. Image registration
Required step. This generates the following 3 files, and also returns the corresponding objects. One can take advantage of parallel processing by specifying the number of cores to be used.

* regimgi.tif: registered fly-view image
* flimgreg.tif: registered fluo-view image
* fvfvwindflregimg.tif: registered window image

```  
registered_images <- register_images(fvimgl=fvimgl,
                                     flimgrt=flimgrt,
                                     fvimgbwbrfh=fvimgbwbrfh,
                                     angles=trj_res$angles,
                                     zoom=zoom,
                                     center=center,
                                     output=output_prefix,
                                     cores=1,
                                     saveRDS=F,
                                     reuse=reuse)
```

### Part 10. Look for good frames based on size, position, motion, and focus
Optional.

```
goodfr <- find_goodframes(window_mask=fvimgbwbrfh,
                          fvimgl=fvimgl,
                          output=output_prefix,
                          reuse=reuse)
```

### Part 11. Calculate fluorescence intensity in the brain window
Required step.
```
message("Measuring fluorescence intensity...")
if(file.exists(paste0(output_prefix, "_intensity_br.RDS"))==T & reuse==T){
  message("Using RDS file")
  intensity_br <- readRDS(paste0(output_prefix, "_intensity_br.RDS"))
}else{
  intensity_br <- colSums(registered_images$fvimgbwbrfhregimg*registered_images$flimgreg, dims=2)/as.integer(goodfr$objsize)
  saveRDS(intensity_br, paste0(dir, prefix, "_intensity_br.RDS"))
}
intensity <- zoo::na.approx(intensity_br)
```

### Part 12. Plot delta F over F0
```
message("Creating dF/F0 plots")
F0int <- mean(intensity[1:5])
deltaFint <- intensity - F0int
dFF0int <- deltaFint/F0int * 100
dat <- data.frame(x=(1:length(dFF0int)), y=dFF0int, d=trj_res$flydist[frida])
p <- ggplot(data=dat, aes(x=x, y=y)) +
  geom_smooth(method="loess", span = 0.4, level=0.95) +
  ylim(-5, 10) +
  geom_line(data=dat, aes(x=x, y=d))
ggsave(filename = paste0(output_prefix, "_dFF0int.pdf"), width = 8, height = 8)
```

### Part 13. Create delta F over F0 pseudocolor images
This step produces the following 4 output files: F0, Fmean, dF over F0 grayscale image, dF over F0 pseudocolor image

```
message("Calculating dF/F0 images...")
dF_F0_image(flimgreg=registered_images$flimgreg,
            fvimgbwbrfhregimg=registered_images$fvimgbwbrfhregimg,
            regimgi=registered_images$regimgi,
            colmax=100, cmin=30, cmax=200,
            goodfr=which(goodfr$goodfr==T),
            output=output_prefix)
```

### Part 14. ROI measurement

```
# Create ROI mask
# Rectangle ROI example
ROI_mask <- array(0, dim=dim(registered_images$flimgreg)[1:2])
ROI_mask[120:(120+10-1),120:(120+10-1)] <- 1
# Circular ROI example
# ROI_mask <- array(0, dim=dim(registered_images$flimgreg)[1:2])
# drawCircle(img=ROI_mask, x=120, y=120, radius=5, col=1, fill=T)

ROI_dFF0 <- measureROI(img=registered_images$flimgreg,
                       mask=ROI_mask,
                       output=output_prefix)
```

### Part 15. Create trajectory of the flies

```
message("Creating trajectory of the flies...")
pdf(file= paste0(output_prefix, "_trackResult.pdf"), width = 4.4, height = 4, bg = "white")
par(plt = c(0, 1, 0, 1), xaxs = "i", yaxs = "i")
plot(trj_res$trja[frida,1]*10, -trj_res$trja[frida,2]*10,
     type = "l", lty = 1, col="red",
     axes = F, xlim = c(-240, 240), ylim = c(-220, 220))
par(new=T)
plotrix::draw.ellipse(0,0,11.0795*20,10*20)
dev.off()
```

### Part 16. Convert fmf to tif format
Optional.

```
if(fmf2tif==T){
  dipr::fmf2tif(paste0(dir, list.files(dir, pattern="^fv.*fmf$")), skip=10)
  dipr::fmf2tif(paste0(dir, list.files(dir, pattern="^av.*fmf$")), skip=2)
}
```
