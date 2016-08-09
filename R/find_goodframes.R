find_goodframes <- function(window_mask, output, reuse=F){
  ftrs <- sfeatures(window_mask, window_mask)
  maxobj <- lapply(ftrs, function(x) x[which(x[, 'm.pxs'] == max(x[, 'm.pxs'])),])
  maxobj[sapply(maxobj, length)==0]<-NA
  objsize <- unlist(lapply(maxobj, function(x) x[c('m.pxs')]))
  objdist <- unlist(lapply(maxobj, function(x) sqrt((x['m.x']-dim(window_mask)[1]/2)^2 + (x['m.y']-dim(window_mask)[2]/2)^2)))
  objx <- unlist(lapply(maxobj, function(x) x[c('m.x')]))
  objy <- unlist(lapply(maxobj, function(x) x[c('m.y')]))
  motion <- sqrt(diff(objx)^2 + diff(objy)^2)
  nomotion <- which(motion < 1)

  png(file=paste0(output, "_motion.png"), width=400, height=400)
  plot(motion)
  dev.off()  png(file=paste0(output, "_objsize.png"), width=400, height=400)
  plot(objsize)
  dev.off()
  png(file=paste0(output, "_objdist.png"), width=400, height=400)
  plot(objdist)
  dev.off()

  objsizemedian <- median(objsize, na.rm=T)
  message(sprintf("Median window size is %.1f", objsizemedian))
  objsizemad <- mad(objsize, na.rm=T)
  goodsizefr <- which(objsize > objsizemedian - 4*objsizemad & objsize < objsizemedian + 4*objsizemad)
  message(sprintf("Good size frame: %d", length(goodsizefr)))
  message(sprintf("The following frames have too big or too small window: %s", paste((1:length(objsize))[-goodsizefr], collapse=" ")))
  goodposfr <- which(objdist < 20)
  message(sprintf("The following frames have too big motion: %s", paste((1:length(motion))[-nomotion], collapse=" ")))
  goodfr <- intersect(nomotion, goodfr)
  message(sprintf("The following frames have too large motion: %s", paste((1:length(objsize))[-gooderrorfr], collapse=" ")))
  goodfr <- intersect(gooderrorfr, goodsizefr)
  message(sprintf("Number of good size and good error frame: %d", length(goodfr)))
  goodfr <- intersect(goodposfr, goodfr)
  message(sprintf("Number of good size, good error, good position frame: %d", length(goodfr)))
  message(sprintf("Good size, good error, good position, good focus, good registration frame: %d", length(goodfr)))
  message(sprintf("The following frames have passed size, centering, motion, and sharpness filters: %s", paste(goodfr, collapse=" ")))
  
  # Detect blurriness
  message("Detecting blurriness...")
  if(file.exists(paste0(output, "_quantcnt.RDS"))==T & reuse==T){
    message("Loading RDS file")
    quantcnt <- readRDS(paste0(output, "_quantcnt.RDS"))
  }else{
    
    LoG <- function(x,y,s){
      fn <- function(x, y, s) 1/(pi*s^4)*((x^2 + y^2)/(2 * s^2)-1)*exp(-(x^2 + y^2)/(2*s^2))
      x <- seq(-floor(x/2), floor(x/2), len = x)
      y <- seq(-floor(y/2), floor(y/2), len = y)
      w <- outer(x, y, fn, s)
      w
    }
    LoGkern <- round(LoG(9,9,1.4)*428.5)
    fvimgllog <- filter2(fvimgl, LoGkern)
    centermask <- drawCircle(fvimgl[,,1]*0, dim(fvimgl)[1]/2, dim(fvimgl)[2]/2, 100, col=1, fill=T)
    fvimgcntlog <- sweep(fvimgllog, 1:2, centermask, FUN="*")
    rm(fvimgllog)
    quantcnt <- apply(fvimgcntlog, 3, function(x) quantile(x, 0.9))
    png(file=paste0(output, "_cntLoG_quant.png"), width=800, height=800)
    plot(quantcnt, ylim=c(0, max(quantcnt)))
    dev.off()
    rm(fvimgcntlog)
    saveRDS(quantcnt, paste0(output, "_quantcnt.RDS"))
    
  }
  
  sharpfr <- which(quantcnt > 150)
  blurfr <- which(quantcnt <= 150)
  message(sprintf("The following frames are blurry: %s", paste(blurfr, collapse=" ")))
  goodfr <- intersect(goodfr, sharpfr)
  message(sprintf("Good size, good error, good position, good focus frame: %d", length(goodfr)))
  return(goodfr)
}