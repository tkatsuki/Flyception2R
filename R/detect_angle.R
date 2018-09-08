#' Detect angle of the beads
#'
#'
#' @param img A binary labelled image or array.
#' @export
#' @examples
#' detect_angle()
#'

detect_angle <- function(img){
  
  ftrs <- list()
  for (i in 1:dim(img)[3]){
    tmp <- computeFeatures.moment(img[,,i])
      if(is.null(tmp)){
        ftrs[[i]] <- 0
      }else{
        ftrs[[i]] <- tmp
      }
  }

  ang <- c()
  centroid <- array(0, dim=c(dim(img)[3],2))
  markernum <- c()
  
  for (im in 1:dim(img)[3]){
    m <- ftrs[[im]]
    if(length(m)==1){
      markernum[im] <- 0
    }else{
      markernum[im] <- nrow(m)
    }
    
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
  
  result <- list(ang, centroid, markernum)
  return(result)
}