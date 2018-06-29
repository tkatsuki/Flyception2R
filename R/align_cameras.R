#' Align camera axis
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' align_cameras()
#'

align_cameras <- function(source, template, output, center=c(0, 0), zoom=1, autopos=T, ROI=F){
  # template image must always be smaller than source image
  # Manual position calibration with fly contour during flash
  if(autopos==F){
    fvrefrs <- EBImage::resize(template, dim(template)[1]*zoom)
    flrefpad <- fvrefrs*0
    if(dim(fvrefrs)[1] > dim(source)[1]){
      flrefpad[round((dim(fvrefrs)[1]-dim(source)[1])/2):
                 (round((dim(fvrefrs)[1]-dim(source)[1])/2)+dim(source)[1]-1),
               round((dim(fvrefrs)[2]-dim(source)[2])/2):
                 (round((dim(fvrefrs)[2]-dim(source)[2])/2)+dim(source)[2]-1)] <- source
    } else{
      flrefpad <- source[round((dim(source)[1]-dim(fvrefrs)[1])/2):
                           (round((dim(source)[1]-dim(fvrefrs)[1])/2)+dim(fvrefrs)[1]-1),
                         round((dim(source)[2]-dim(fvrefrs)[2])/2):
                           (round((dim(source)[2]-dim(fvrefrs)[2])/2)+dim(fvrefrs)[2]-1)]
    }
    flrefpadmv <- EBImage::translate(flrefpad, center)
    display(normalize(flrefpadmv))
    display(normalize(fvrefrs + flrefpadmv))
    writeImage(normalize(fvrefrs + flrefpadmv), file=paste0(output, "_aligned.png"))
    
  } else {
    # Automated position calibration using template matching
    message("Automatically aligning two cameras...")
    writeImage(source, file=paste0(output, "_source.png"))
    writeImage(template, file=paste0(output, "_template.png"))
    source_rs <- EBImage::resize(source, dim(source)[1]*zoom)
    
    if(dim(source_rs)[1] > dim(template)[1] && dim(source_rs)[2] > dim(template)[2]){
      fncc <- dipr::FNCC(source_rs, template)
      if(length(ROI)!=1){
        fncc <- fncc[ROI[1]:(ROI[1]+ROI[3]-1),
                     ROI[2]:(ROI[2]+ROI[4]-1)]
      }
      maxpeak <- which(fncc==max(fncc), arr.ind=TRUE)
      centerx <- (maxpeak[1] + round(nrow(template)/2)) - round(dim(source_rs)[1]/2)
      centery <- (maxpeak[2] + round(ncol(template)/2)) - round(dim(source_rs)[2]/2)
      center <- c(centerx, centery)
      template_pad <- source_rs*0
      template_pad[(round((dim(source_rs)[1]-dim(template)[1])/2) + 1):
                     (round((dim(source_rs)[1]-dim(template)[1])/2)+dim(template)[1]),
                   (round((dim(source_rs)[2]-dim(template)[2])/2) + 1):
                     (round((dim(source_rs)[2]-dim(template)[2])/2)+dim(template)[2])] <- template
    }else if (dim(source_rs)[1] < dim(template)[1]){
      stop(paste("Use smaller template"))
    } else{
      # Template needs to be smaller than the reference therefore crop one of the images
      x1 <- round(nrow(template)/2)-round(nrow(template)/4)
      x2 <- round(nrow(template)/2)+round(nrow(template)/4)
      y1 <- round(ncol(template)/2)-round(ncol(template)/4)
      y2 <- round(ncol(template)/2)+round(ncol(template)/4)
      templatecrop <- template[x1:x2, y1:y2]
      fncc <- dipr::FNCC(source_rs, templatecrop)
      if(length(ROI)!=1){
        fncc <- fncc[ROI[1]:(ROI[1]+ROI[3]-1),
                     ROI[2]:(ROI[2]+ROI[4]-1)]
      }
      maxpeak <- which(fncc==max(fncc), arr.ind=TRUE)
      centerx <- (maxpeak[1] + round(nrow(templatecrop)/2)) - round(dim(source_rs)[1]/2)
      centery <- (maxpeak[2] + round(ncol(templatecrop)/2)) - round(dim(source_rs)[2]/2)
      center <- c(centerx, centery)
      template_pad <- source_rs*0
      template_pad[round((dim(source_rs)[1]-dim(templatecrop)[1])/2):
                     (round((dim(source_rs)[1]-dim(templatecrop)[1])/2)+dim(templatecrop)[1]-1),
                   round((dim(source_rs)[2]-dim(templatecrop)[2])/2):
                     (round((dim(source_rs)[2]-dim(templatecrop)[2])/2)+dim(templatecrop)[2]-1)] <- templatecrop
    }
    
    template_mv <- EBImage::translate(template_pad, center)
    EBImage::writeImage(normalize(source_rs), file=paste0(output, "_source_aligned.png"))
    EBImage::writeImage(normalize(template_mv), file=paste0(output, "_template_aligned.png"))
    EBImage::writeImage(normalize(fncc), file=paste0(output, "_fncc.png"))
    
  }
  message(sprintf("FNCC max peak was: x=%d, y=%d", maxpeak[1], maxpeak[2]))
  message(sprintf("Center offset: x=%d, y=%d", center[1], center[2]))
  return(center)
}
