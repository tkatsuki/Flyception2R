library(FlyceptionR)
library(zoo)
library(magick)
library(loggit)
source("~/Flyception2R/R/align_cameras.R")
source("~/Flyception2R/R/imageJ_crop_append.R")
source("~/Flyception2R/R/sync_frames.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/align_cameras.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/imageJ_crop_append.R")
source("C:/Users/tkatsuki/Documents/GitHub/Flyception2R/R/sync_frames.R")

dirlist <- list.dirs("/Volumes/LaCie/P1_GCaMP6s_tdTomato_06182018_CW_Dual_Laser", recursive=F)

for(dl in dirlist){
  dir <- paste0(dl, "//")
  prefix <- strsplit(dir, "/")[[1]][5]       # Will be used as a filename prefix
  
  autopos <- T             # True if you want to align cameras automatically 
  zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
  ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
  fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
  fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
  av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
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
                          autopos=T,
                          ROI=c(1, 1, 450, 50))
  
  
  if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
    imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
  }
  fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
}

dirlist <- list.dirs("/Volumes/LaCie/P1_GCaMP6s_tdTomato_06202018_CW_Dual_Laser", recursive=F)

for(dl in dirlist){
  dir <- paste0(dl, "//")
  prefix <- strsplit(dir, "/")[[1]][5]       # Will be used as a filename prefix
  
  autopos <- T             # True if you want to align cameras automatically 
  zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
  ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
  fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
  fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
  av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
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
                          autopos=T,
                          ROI=c(1, 1, 450, 50))
  
  
  if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
    imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
  }
  fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
}

dirlist <- list.dirs("/Volumes/LaCie/P1_GCaMP6s_tdTomato_06112018_CW_Dual_Laser", recursive=F)

for(dl in dirlist){
  dir <- paste0(dl, "//")
  prefix <- strsplit(dir, "/")[[1]][5]       # Will be used as a filename prefix
  
  autopos <- T             # True if you want to align cameras automatically 
  zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
  ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
  fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
  fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
  av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
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
                          autopos=T,
                          ROI=c(1, 1, 450, 50))
  
  
  if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
    imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
  }
  fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
}

dirlist <- list.dirs("/Volumes/LaCie/P1_GCaMP6s_tdTomato_06082018_CW_Dual_Laser", recursive=F)

for(dl in dirlist){
  dir <- paste0(dl, "//")
  prefix <- strsplit(dir, "/")[[1]][5]       # Will be used as a filename prefix
  
  autopos <- T             # True if you want to align cameras automatically 
  zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
  ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
  fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
  fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
  av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
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
                          autopos=T,
                          ROI=c(1, 1, 450, 50))
  
  
  if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
    imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
  }
  fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
}

dirlist <- list.dirs("/Volumes/LaCie/P1_GCaMP6s_tdTomato_06022018_CW_Dual_Laser", recursive=F)

for(dl in dirlist){
  dir <- paste0(dl, "//")
  prefix <- strsplit(dir, "/")[[1]][5]       # Will be used as a filename prefix
  
  autopos <- T             # True if you want to align cameras automatically 
  zoom <- 1.085             # Zoom ratio: fluo-view/fly-view. Measure this using a resolution target.
  ROI <- c(391, 7, 240, 240) # Top left corner is (0, 0)
  fluo_flash_thresh <- 500 # Threshold for detecting flash in fluo-view
  fv_flash_thresh <- 240   # Threshold for detecting flash in fly-view
  av_flash_thresh <- 100    # Threshold for detecting flash in arena-view
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
                          autopos=T,
                          ROI=c(1, 1, 450, 50))
  
  
  if(length(list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))==0){
    imageJ_crop_append(dir, ch=2, roi=c((1024 + ROI[1] + center[1]), (ROI[2] + center[2]), 240, 240)) # x and y coordinates of the top left corner, width, height
  }
  fluo_view_tif_ch2 <- paste0(dir, list.files(dir, pattern="ome\\.ch2\\.crop\\.concat\\.tif$"))
}
