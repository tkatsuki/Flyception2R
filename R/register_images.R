#' Perform image registration for fly-view, window, and fluo-view
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' register_images()
#'

register_images <- function(fvimgl, flimgrt, fvimgbwbrfh, angles, output, zoom, center, cores=1, saveRDS=F, reuse=F){
  message("Performing image registration...")
  if(file.exists(paste0(output, "_regimgi.RDS"))==T &
     file.exists(paste0(output, "_regresi.RDS"))==T & reuse==T){
    message("Loading RDS file")
    regimgi <- readRDS(paste0(output, "_regimgi.RDS"))
    regresi <- readRDS(paste0(output, "_regresi.RDS"))
  }else{
    # Prepare fvimg for registration
    fvimgli <- resize(255 - fvimgl, dim(fvimgl)[1]*zoom)
    centermask <- drawCircle(matrix(0,dim(fvimgli)[1],dim(fvimgli)[2]), dim(fvimgli)[1]/2,
                             dim(fvimgli)[2]/2, dim(fvimgli)[1]/2-1, col=1, fill=1)
    # Create first image, which will be the target in registration
    fvimgrt1sti <- rotate(fvimgli[,,1], angles[frid[1]]*180/pi, output.dim=dim(fvimgli)[1:2])
    # Prepare cores
    if(cores!=1){
      library(doParallel)
      registerDoParallel(cores=cores)
      message(print0(getDoParWorkers(), " cores are used..."))
    }
    # Build affine matrix for rotation
    aff <- list()
    for(a in 1:dim(fvimgli)[3]){
      aff[[a]] <- buildAffine(angles=c(0,0, -angles[frid[a]]), source=fvimgli[,,1], anchor="center")
    }
    # Run image registration using the initial angles
    regresi <- list()
    if(cores==1){
      for(rg in 1:dim(fvimgli)[3]){
        regresi[[rg]] <- niftyreg(fvimgli[,,rg], fvimgrt1sti,
                                  init=aff[[rg]], scope="rigid", symmetric=F)
      }
    }else{
      regresi <- foreach(rg = 1:dim(fvimgli)[3]) %dopar% niftyreg(fvimgli[,,rg], fvimgrt1sti, init=aff[[rg]], scope="rigid", symmetric=F, internal=FALSE)
    }
    regimgi <- array(sapply(regresi, function(x) x$image), dim=dim(fvimgli))
    regimgi[which(is.na(regimgi)==T)] <- 0
    if(saveRDS==T){
      saveRDS(regimgi, paste0(output, "_regimgi.RDS"))
      saveRDS(regresi, paste0(output, "_regresi.RDS"))
    }
    writeImage((255-regimgi)/255, file=paste0(output, "_regimgi.tif"))
    rm(fvimgli)
  }

  # Apply affine transformation to fluo-view
  message("Transforming fluo-view...")
  if(file.exists(paste0(output, "_flimgreg.RDS"))==T & reuse==T){
    message("Loading RDS file")
    flimgreg <- readRDS(paste0(output, "_flimgreg.RDS"))

  }else{

    # Pad flimgrt to match the size of fvimg
    if(dim(regimgi)[1] > dim(flimgrt)[1]){
      flimgpad <- regimgi*0
      flimgpad[round((dim(regimgi)[1]-dim(flimgrt)[1])/2):
                 (round((dim(regimgi)[1]-dim(flimgrt)[1])/2)+dim(flimgrt)[1]-1),
               round((dim(regimgi)[2]-dim(flimgrt)[2])/2):
                 (round((dim(regimgi)[2]-dim(flimgrt)[2])/2)+dim(flimgrt)[2]-1),
               1:dim(flimgrt)[3]] <- flimgrt

    }else{
      flimgpad <- flimgrt[round((dim(flimgrt)[1]-dim(regimgi)[1])/2):
                            (round((dim(flimgrt)[1]-dim(regimgi)[1])/2)+dim(regimgi)[1]-1),
                          round((dim(flimgrt)[2]-dim(regimgi)[2])/2):
                            (round((dim(flimgrt)[2]-dim(regimgi)[2])/2)+dim(regimgi)[2]-1),]
    }
    flimgpadmv <- translate(flimgpad, center)
    flimgrgres <- list()
    if(cores==1){
      for(app in 1:dim(fvimgl)[3]){
        flimgrgres[[app]] <- applyTransform(forward(regresi[[app]]), flimgpadmv[,,app])
      }
    }else{
      flimgrgres <- foreach(app = 1:dim(fvimgl)[3]) %dopar%  applyTransform(forward(regresi[[app]]), flimgpadmv[,,app])
    }
    flimgreg <- array(unlist(flimgrgres), dim=dim(flimgpadmv))
    writeImage(flimgreg, file=paste0(output, "_flimgreg.tif"))
    rm(flimgrgres)
    rm(flimgpad)
    rm(flimgpadmv)
    if(saveRDS==T){
      saveRDS(flimgreg, paste0(output, "_flimgreg.RDS"))
    }
  }

  # Apply affine transformation to window mask
  message("Transforming window images...")
  if(file.exists(paste0(output, "_fvimgbwbrfhregimg.RDS"))==T & reuse==T){
    message("Loading RDS file")
    fvimgbwbrfhregimg <- readRDS(paste0(output, "_fvimgbwbrfhregimg.RDS"))
  }else{

    fvimgbwbrfhrs <- resize(fvimgbwbrfh, dim(fvimgl)[1]*zoom)
    rm(fvimgbwbrfh)
    fvimgbwbrfhreg <- list()
    if(cores==1){
      for(app in 1:dim(fvimgl)[3]){
        fvimgbwbrfhreg[[app]] <- applyTransform(forward(regresi[[app]]), fvimgbwbrfhrs[,,app])
      }
    }else{
      fvimgbwbrfhreg <- foreach(app = 1:dim(fvimgl)[3]) %dopar% applyTransform(forward(regresi[[app]]), fvimgbwbrfhrs[,,app])
    }
    fvimgbwbrfhregimg <- array(unlist(fvimgbwbrfhreg), dim=dim(fvimgbwbrfhrs))
    fvimgbwbrfhregimg <- fvimgbwbrfhregimg >= 0.5
    writeImage((255-regimgi)/255*(1-fvimgbwbrfhregimg) + normalize(flimgreg),
               file=paste0(output, "_fvfvwindflregimg.tif"))
    rm(fvimgbwbrfhrs)
    rm(fvimgbwbrfhreg)
    rm(regresi)
    if(saveRDS==T){
      saveRDS(fvimgbwbrfhregimg, paste0(output, "_fvimgbwbrfhregimg.RDS"))
    }
  }
  return(list("regimgi"=regimgi, "flimgreg"=flimgreg, "fvimgbwbrfhregimg"=fvimgbwbrfhregimg))
}
