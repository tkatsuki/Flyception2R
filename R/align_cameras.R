#' Align camera axis
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' align_cameras()
#'

align_cameras <- function(flref, fvref, output, center=c(0, 0), zoom=1, autopos=T){
  # Manual position calibration with fly contour during flash
  if(autopos==F){
    EBImage::display(flref)
    EBImage::display(fvref)
    fvrefrs <- EBImage::resize(fvref, dim(fvref)[1]*zoom)
    flrefpad <- fvrefrs*0
    if(dim(fvrefrs)[1] > dim(flref)[1]){
      print(1)
      flrefpad[round((dim(fvrefrs)[1]-dim(flref)[1])/2):
                 (round((dim(fvrefrs)[1]-dim(flref)[1])/2)+dim(flref)[1]-1),
               round((dim(fvrefrs)[2]-dim(flref)[2])/2):
                 (round((dim(fvrefrs)[2]-dim(flref)[2])/2)+dim(flref)[2]-1)] <- flref
    } else{
      flrefpad <- flref[round((dim(flref)[1]-dim(fvrefrs)[1])/2):
                          (round((dim(flref)[1]-dim(fvrefrs)[1])/2)+dim(fvrefrs)[1]-1),
                        round((dim(flref)[2]-dim(fvrefrs)[2])/2):
                          (round((dim(flref)[2]-dim(fvrefrs)[2])/2)+dim(fvrefrs)[2]-1)]
    }
    flrefpadmv <- EBImage::translate(flrefpad, center)
    EBImage::display(flrefpadmv)
    EBImage::display(EBImage::normalize(fvrefrs + flrefpadmv))
    EBImage::writeImage(normalize(fvrefrs + flrefpadmv), file=paste0(output, "_aligned.png"))

  } else {
    # Automated position calibration using template matching
    message("Automatically aligning two cameras...")
    EBImage::writeImage(flref, file=paste0(output, "_flref.png"))
    EBImage::writeImage(fvref, file=paste0(output, "_fvref.png"))
    fvrefrs <- EBImage::resize(fvref, dim(fvref)[1]*zoom)
    if(dim(fvrefrs)[1] > dim(flref)[1]){
      fncc <- dipr::FNCC(fvrefrs, flref)
      maxpeak <- which(fncc==max(fncc), arr.ind=TRUE)
      centerx <- (maxpeak[1] + round(nrow(flref)/2)) - round(dim(fvrefrs)[1]/2)
      centery <- (maxpeak[2] + round(ncol(flref)/2)) - round(dim(fvrefrs)[2]/2)
      center <- c(centerx, centery)
      flrefpad <- fvrefrs*0
      flrefpad[round((dim(fvrefrs)[1]-dim(flref)[1])/2):
                 (round((dim(fvrefrs)[1]-dim(flref)[1])/2)+dim(flref)[1]-1),
               round((dim(fvrefrs)[2]-dim(flref)[2])/2):
                 (round((dim(fvrefrs)[2]-dim(flref)[2])/2)+dim(flref)[2]-1)] <- flref
    }else{
      fncc <- dipr::FNCC(flref, fvrefrs)
      maxpeak <- which(fncc==max(fncc), arr.ind=TRUE)
      centerx <- (maxpeak[1] + round(nrow(flref)/2)) - round(dim(fvrefrs)[1]/2)
      centery <- (maxpeak[2] + round(ncol(flref)/2)) - round(dim(fvrefrs)[2]/2)
      center <- c(centerx, centery)
      flrefpad <- fvrefrs*0
      flrefpad <- flref[round((dim(flref)[1]-dim(fvrefrs)[1])/2):
                          (round((dim(flref)[1]-dim(fvrefrs)[1])/2)+dim(fvrefrs)[1]-1),
                        round((dim(flref)[2]-dim(fvrefrs)[2])/2):
                          (round((dim(flref)[2]-dim(fvrefrs)[2])/2)+dim(fvrefrs)[2]-1)]
    }

    flrefpadmv <- EBImage::translate(flrefpad, center)
    EBImage::writeImage(normalize(fvrefrs + flrefpadmv), file=paste0(output, "_aligned.png"))

  }
  message(sprintf("Center offset: x=%d, y=%d", center[1], center[2]))
  return(center)
}
