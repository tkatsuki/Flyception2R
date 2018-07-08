#' Create F_ratio images
#'
#'
#' @param redwindowmed Red channel image.
#' @param redwindowmedth A binary red channel image.
#' @param greenwindowmed Green channel image.
#' @export
#' @examples
#' F_ratio_image()
#'

F_ratio_image <- function(redwindowmed, redwindowmedth, greenwindowmed){
  
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
  
}