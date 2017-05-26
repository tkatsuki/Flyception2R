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
  windir <- gsub("/", "\\", dir, fixed=TRUE)
    
  # Configure ImageJ so that TIFF is saved in little endian
  macro <- paste0('run("Input/Output...", "jpeg=85 gif=-1 file=.xls use_file save copy_row save_column save_row");\n')
  write(macro, file=paste0(dir,"macro.txt"))
  if (os == "windows"){
    bat <- paste0('pushd "C:\\Program Files\\ImageJ"', '\n jre\\bin\\java -jar -Xmx8g ij.jar  -batch "', windir, 'macro.txt" ', dir, '\n pause\n exit')
    tempbat <- paste(tempfile('bat'),".bat",sep="")
    write(bat, file=tempbat)
    shell(tempbat, translate=T, wait=T)   
  }else{
    system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar -ijpath /Applications/ImageJ -batch ", dir, "macro.txt"), wait=T) 
  }

  # Crop a ROI for each file
  fluo_view_files <- list.files(dir, pattern="ome\\.tif$", full.names=T) # The first file is created without a numeric extension
  file_order <- c(length(fluo_view_files), 1:(length(fluo_view_files)-1)) # Therefore we need to bring it front
  for(s in file_order){
    macro <- paste0('open("',fluo_view_files[s],'");\n makeRectangle(', paste(roi, collapse=","), ');\n run("Crop");\n saveAs("tiff", "',tools::file_path_sans_ext(fluo_view_files[s]),'.ch',ch,'.crop.tif");\n run("Quit");\n')
    write(macro, file=paste0(dir,"macro.txt"))
    if (os == "windows"){
      bat <- paste0('pushd "C:\\Program Files\\ImageJ"', '\n jre\\bin\\java -jar -Xmx8g ij.jar  -batch "', windir, 'macro.txt" ', dir, '\n pause\n exit')
      tempbat <- paste(tempfile('bat'),".bat",sep="")
      write(bat, file=tempbat)
      shell(tempbat,wait=T)   
    }else{
      system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar -ijpath /Applications/ImageJ -batch ", dir, "macro.txt"), wait=T) 
    }
  }
  
  # Load cropped videos
  fluo_view_cropped_files_full <- list.files(dir, pattern=paste0("ome\\.ch", ch, "\\.crop\\.tif$"),  full.names=T)
  fluo_view_cropped_files <- list.files(dir, pattern=paste0("ome\\.ch", ch, "\\.crop\\.tif$"))
  cropped_file_order <- c(length(fluo_view_cropped_files), 1:(length(fluo_view_cropped_files)-1))
  
  for(cr in 1:length(fluo_view_cropped_files_full)){
    if(cr == 1){
      write(paste0('open("',fluo_view_cropped_files_full[cropped_file_order[cr]],'");\n'), file=paste0(dir,"macro.txt"))
    }else{
      write(paste0('open("',fluo_view_cropped_files_full[cropped_file_order[cr]],'");\n'), file=paste0(dir,"macro.txt"), append=T)
    }
  }
  
  # Concatenate cropped videos into one
  strs <- c()
  for(st in 1:length(fluo_view_cropped_files)){
    strs[st] <- paste0('image',st,'=',fluo_view_cropped_files[cropped_file_order[st]])
  }
  write(paste('run("Concatenate...", "  title=[Concatenated Stacks]', paste(strs, collapse=" "), '");\n'), file=paste0(dir,"macro.txt"), append=T)
  write(paste0('saveAs("tiff", "', tools::file_path_sans_ext(fluo_view_cropped_files_full[length(fluo_view_cropped_files_full)]), '.concat.tif");\n run("Quit");\n'), file=paste0(dir,"macro.txt"), append=T)
  
  # Execute the macro
  
  if (os == "windows"){
    bat <- paste0('pushd "C:\\Program Files\\ImageJ"', '\n jre\\bin\\java -jar -Xmx8g ij.jar  -batch "', windir, 'macro.txt" ', dir, '\n pause\n exit')
    tempbat <- paste(tempfile('bat'),".bat",sep="")
    write(bat, file=tempbat)
    shell(tempbat,wait=T)   
  }else{

  system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar -ijpath /Applications/ImageJ -batch macro.txt"), wait=T) 
  }
}

