#' Crop and concatenate fluo_view image with up to 2 channels using ImageJ
#'
#' @param dir path to the directory that contains the data
#' @param ch channel name used for output files
#' @param roi roi coordinates (x, y, w, h), note origin is (0,0)
#' @export
#' @examples
#' imageJ_crop_append()

imageJ_crop_append <- function(dir, ch=1, roi=c(383, 0, 256, 256)){
  
  os <-.Platform$OS.type
  
  # Crop a ROI for each file
  fluo_view_files <- list.files(dir, pattern="ome\\.tif$", full.names=T) # The first file is created without a numeric extension
  file_order <- c(length(fluo_view_files), 1:(length(fluo_view_files)-1)) # Therefore we need to bring it front
  for(s in file_order){
    macro <- paste0('open("',fluo_view_files[s],'");\n makeRectangle(', paste(roi, collapse=","), ');\n run("Crop");\n saveAs("tiff", "',tools::file_path_sans_ext(fluo_view_files[s]),'.ch',ch,'.crop.tif");\n run("Quit");\n')
    write(macro, file=paste0(dir,"macro.txt"))
    if (os == "windows"){
      bat <- paste0('pushd "C:\\Program Files\\ImageJ"', '\n jre\\bin\\java -jar -Xmx1024m ij.jar  -batch "H:\\P1_GCaMP6s_tdTomato_041117\\P1-Gal4_UAS-GCaMP6s_tdTomato_test_8\\macro.txt" ', dir, '\n pause\n exit')
      tempbat <- paste(tempfile('bat'),".bat",sep="")
      write(bat, file=tempbat)
      shell(tempbat,wait=T)   
    }else{
      system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar -ijpath /Applications/ImageJ -batch ", dir, "macro.txt"), wait=T) 
    }
  }
  
  # Load cropped videos
  fluo_view_cropped_files <- list.files(dir, pattern=paste0("ome\\.ch", ch, "\\.crop\\.tif$"))
  cropped_file_order <- c(length(fluo_view_cropped_files), 1:(length(fluo_view_cropped_files)-1))
  
  for(cr in 1:length(fluo_view_cropped_files)){
    if(cr == 1){
      write(paste0('open("',fluo_view_cropped_files[cropped_file_order[cr]],'");\n'), file="macro.txt")
    }else{
      write(paste0('open("',fluo_view_cropped_files[cropped_file_order[cr]],'");\n'), file="macro.txt", append=T)
    }
  }
  
  # Concatenate cropped videos into one
  strs <- c()
  for(st in 1:length(fluo_view_cropped_files)){
    strs[st] <- paste0('image',st,'=',fluo_view_cropped_files[cropped_file_order[st]])
  }
  write(paste('run("Concatenate...", "  title=[Concatenated Stacks]', paste(strs, collapse=" "), '");\n'), file="macro.txt", append=T)
  write(paste0('saveAs("tiff", "', tools::file_path_sans_ext(fluo_view_cropped_files[length(fluo_view_cropped_files)]), '.concat.tif");\n run("Quit");\n'), file="macro.txt", append=T)
  
  # Execute the macro
  system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar -ijpath /Applications/ImageJ -batch macro.txt"), wait=T) 
  
}

