#' Calculate window size, position, motion, and blurriness
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' find_goodframes()
#'

find_goodframes <- function(window_mask, fvimgl, output, motion_thresh=10, dist_thresh=20, sharp_thresh=150, reuse=F){
  # This function calculate 4 variables to assess suitability of each frame

  # Compute features
  ftrs <- sfeatures(window_mask, window_mask)
  # Detect a biggest object in each frame
  maxobj <- lapply(ftrs, function(x) x[which(x[, 'm.pxs'] == max(x[, 'm.pxs'])),])
  # NA if no object was present
  maxobj[sapply(maxobj, length)==0]<-NA
  # Reshape into a vector of object size
  objsize <- unlist(lapply(maxobj, function(x) x[c('m.pxs')]))
  # Compute the distance of the object from the center of the image
  objdist <- unlist(lapply(maxobj, function(x) sqrt((x['m.x']-dim(window_mask)[1]/2)^2 + (x['m.y']-dim(window_mask)[2]/2)^2)))
  # Compute the motion from the previous frame
  objx <- unlist(lapply(maxobj, function(x) x[c('m.x')]))
  objy <- unlist(lapply(maxobj, function(x) x[c('m.y')]))
  motion <- sqrt(diff(objx)^2 + diff(objy)^2)

  objsizemedian <- median(objsize, na.rm=T)
  message(sprintf("Median window size is %.1f", objsizemedian))
  objsizemad <- mad(objsize, na.rm=T)
  goodsizefr <- which(objsize > objsizemedian - 4*objsizemad & objsize < objsizemedian + 4*objsizemad)
  message(sprintf("The following frames have too big or too small window: %s", paste((1:length(objsize))[-goodsizefr], collapse=" ")))
  goodposfr <- which(objdist < dist_thresh)
  message(sprintf("The following frames have the window too far from the center: %s", paste((1:length(objdist))[-goodposfr], collapse=" ")))
  goodmotionfr <- which(motion < motion_thresh)
  message(sprintf("The following frames have too large motion: %s", paste((1:length(motion))[-goodmotionfr], collapse=" ")))

  # Detect blurriness
  message("Detecting blurriness...")
  if(file.exists(paste0(output, "_quantcnt.RDS"))==T & reuse==T){
    message("Loading RDS file")
    quantcnt <- readRDS(paste0(output, "_quantcnt.RDS"))
  }else{
    LoGkern <- round(LoG(9,9,1.4)*428.5)
    fvimgllog <- filter2(fvimgl, LoGkern)
    centermask <- drawCircle(fvimgl[,,1]*0, dim(fvimgl)[1]/2, dim(fvimgl)[2]/2, 100, col=1, fill=T)
    fvimgcntlog <- ssweep(fvimgllog, centermask, op="*")
    quantcnt <- apply(fvimgcntlog, 3, function(x) quantile(x, 0.9))
    saveRDS(quantcnt, paste0(output, "_quantcnt.RDS"))
  }
  sharpfr <- which(quantcnt > sharp_thresh)
  message(sprintf("The following frames are blurry: %s", paste((1:length(quantcnt))[-sharpfr], collapse=" ")))

  # Intersect results
  goodfr <- Reduce(intersect, list(goodsizefr, goodposfr, goodmotionfr, sharpfr))
  message(sprintf("The following frames have passed size, centering, motion, and sharpness filters: %s", paste(goodfr, collapse=" ")))

  # Summarize results in a dataframe
  df <- data.frame(objsize=objsize, objdist=objdist, motion=motion, blurriness=quantcnt, goodfr=1:length(motion)[goodfr])

  # Plot results
  png(file=paste0(output, "_motion.png"), width=400, height=400)
  plot(motion)
  dev.off()
  png(file=paste0(output, "_objsize.png"), width=400, height=400)
  plot(objsize)
  dev.off()
  png(file=paste0(output, "_objdist.png"), width=400, height=400)
  plot(objdist)
  dev.off()
  png(file=paste0(output, "_cntLoG_quant.png"), width=800, height=800)
  plot(quantcnt, ylim=c(0, max(quantcnt)))
  dev.off()

  # Output results to a file
  write.table(df, file=paste0(output, "_goodfr_df.txt"), row.names=F)
  return(df)
}
