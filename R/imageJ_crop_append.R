#' Crop and concatenate fluo_view image with up to 2 channels using ImageJ
#'
#' @param dir path to the directory that contains the data
#' @param outdir path to directory in where to save outputs (defaults to dir)
#' @param ch channel name used for output files
#' @param roi roi coordinates (x, y, w, h), note origin is (0,0)
#' @export
#' @examples
#' imageJ_crop_append()

imageJ_crop_append <- function(dir, outdir=NA, ch=1, roi=c(383, 0, 256, 256)){
  
  os <-.Platform$OS.type
  
  if(is.na(outdir))
    outdir <- dir
  
  windir <- gsub("/", "\\", outdir, fixed=TRUE)

    # Crop a ROI for each file
  fluo_view_files <- list.files(dir, pattern="ome\\.tif$", full.names=T) # The first file is created without a numeric extension
  if(length(fluo_view_files)==1){
    file_order <- 1
  }else{
    file_order <- c(length(fluo_view_files), 1:(length(fluo_view_files)-1)) # Therefore we need to bring it front
  }
  for(s in file_order){
    macro <- paste0('run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file save copy_row save_column save_row");',
                    '\nopen("', fluo_view_files[s],'");\n',
                    'makeRectangle(', paste(roi, collapse=","),');\n',
                    'run("Crop");\n saveAs("tiff", "',
                    paste0(outdir,basename(tools::file_path_sans_ext(fluo_view_files[s]))),
                    '.ch',ch,'.crop.tif");\n',
                    'run("Quit");\n')
    
    write(macro, file=paste0(outdir,"macro1.txt"))
    
    if (os == "windows"){
      bat <- paste0('pushd "C:\\Program Files\\ImageJ"',
                    '\n jre\\bin\\java -jar -Xmx8g ij.jar  -batch "',
                    windir, 'macro1.txt" ', dir, '\n pause\n exit')
      tempbat <- paste(tempfile('bat'),".bat",sep="")
      write(bat, file=tempbat)
      shell(tempbat,wait=T)   
    } else if (os == "unix") {
      system(paste0("~/Fiji.app/ImageJ-linux64 --headless -batch ", outdir, "macro1.txt"), wait=T) 
    }else{
      system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar",
                    " -ijpath /Applications/ImageJ ","--headless -batch ", outdir, "macro1.txt"), wait=T) 
    }
  } 
  
  # Load cropped videos
  fluo_view_cropped_files_full <- list.files(outdir, pattern=paste0("ome\\.ch", ch, "\\.crop\\.tif$"),  full.names=T)
  fluo_view_cropped_files <- list.files(outdir, pattern=paste0("ome\\.ch", ch, "\\.crop\\.tif$"))
  if(length(fluo_view_cropped_files)==1){
    return()
  }else{
    cropped_file_order <- c(length(fluo_view_cropped_files), 1:(length(fluo_view_cropped_files)-1))
  }
  
  for(cr in 1:length(fluo_view_cropped_files_full)){
    if(cr == 1){
      write(paste0('open("',fluo_view_cropped_files_full[cr],'");'), file=paste0(outdir,"macro2.txt"))
    }else{
      write(paste0('open("',fluo_view_cropped_files_full[cr],'");'), file=paste0(outdir,"macro2.txt"), append=T)
    }
  }
  # Concatenate cropped videos into one
  if (os == "windows"){
    strs <- c()
    for(st in 1:length(fluo_view_cropped_files)){
      strs[st] <- paste0('image',st,'=',fluo_view_cropped_files[st])
    }
  } else {
    strs <- c()
    for(st in 1:length(fluo_view_cropped_files)){
      strs[st] <- paste0('image',st,'=',fluo_view_cropped_files[cropped_file_order[st]])
    }
  }
  write(paste('run("Concatenate...", "  title=[Concatenated Stacks]',
              paste(strs, collapse=" "), '");'),
        file=paste0(outdir,"macro2.txt"), append=T)
  write('run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file save copy_row save_column save_row");',
        file=paste0(outdir,"macro2.txt"), append=T)
  write(paste0('saveAs("tiff", "',
               paste0(outdir,basename(tools::file_path_sans_ext(fluo_view_cropped_files_full[1]))),
               '.concat.tif");\nrun("Quit");'),
        file=paste0(outdir,"macro2.txt"), append=T)
  
  # Execute the macro
  if (os == "windows"){
    bat <- paste0('pushd "C:\\Program Files\\ImageJ"',
                  '\n jre\\bin\\java -jar -Xmx8g ij.jar  -batch "',
                  windir, 'macro2.txt" ', outdir, '\n pause\n exit')
    tempbat <- paste(tempfile('bat'),".bat",sep="")
    write(bat, file=tempbat)
    shell(tempbat,wait=T)
    
  } else if (os == "unix") {
    system(paste0("~/Fiji.app/ImageJ-linux64 --headless -batch ", outdir, "macro2.txt"), wait=T) 
  } else {
    print(os)
    system(paste0("java -Xmx8g -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar ",
                  "-ijpath /Applications/ImageJ --headless -batch ", outdir, "macro2.txt"), wait=T) 
  }
  
  # Remove intermedite files
  for(i in 1:length(fluo_view_cropped_files_full))
    file.remove(fluo_view_cropped_files_full[i])
  
}