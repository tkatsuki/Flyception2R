# Flyception2R
R scripts and utilities for analyzing Flyception2 data

## Requirement
Windows: R and Rtools

Mac: R and Xcode

ImageJ https://imagej.nih.gov/ij/

Set ImageJ (and Fiji if present) to save Tiff in Intel byte order from Edit -> Options -> Input/Output... 

## Installation
The following commands will install packages necessary for running Flyception2R.

```
install.packages(c("devtools", "ggplot2", "RNiftyReg", "zoo", "loggit"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EBImage")
library(devtools)
install.Rtools() # only for Windows
devtools::install_github("tkatsuki/dipr")
devtools::install_github("tkatsuki/Flyception2R")
library(Flyception2R)
```

## Usage example
The function Flyception2R processes and analyzes the data automatically or interactively depending on the inputs to the command. Normally, you just need to specify which frames of a video file you want to analyze. Then the results will be saved in a folder named with the frame range.

To run the following example, first download the data (~13GB) from the link below. The data consists of 100sec recordings from a P1-GCaMP6s,tdTomato male fly.

https://www.dropbox.com/sh/okey9zkuhpw44zb/AABmV_mT09rPtAWEU7lgcd-Ua?dl=0

```
dir <- "/PATH/TO/DATA/P1_GCaMP6s_tdTomato_06182018_CW_Dual_Laser/P1-Gal4_UAS-GCaMP6s_tdTomato_12/"
Flyception2R(dir=dir, FOI=c(4242, 4556), interaction=F, flash=1)
```

Because the marker position, hence the center of the image, relative to the fly brain varies from fly to fly, the position of the window (ROI) for segmenting neurons needs to be manually adjusted. Flyception2R will ask you if the window_size and window_offset are acceptable. Check the file ```_redwindow.tif``` and adjust the window size or offset so that the neurons of your interest are in the window (see example below). When you want to move the window up, decrease the y value. When you want to move the window left, decrease the x value. You can create multiple ROIs if you wish.

![redwindow](https://github.com/tkatsuki/Flyception2R/blob/master/redwindow.png)

```
[1] "Current window_size is x=68 y=28"
[1] "Current window_offset is x=-4 y=12"
Check redwindow.tif. Is the window size good (Y or N)?:Y
Check redwindow.tif. Is the window offset good (Y or N)?:N
Enter new x offset:0
Enter new y offset:0
```

The function spits out results like so:

```
window_size was x=68 y=28
window_offset was x=-4 y=12
FOI was from 4242 to 4556
Max F_ratio intensity in this bout was 0.719705128074784
Number of good frames was 197
||c(4242, 4556) ||c(68, 28) ||c(-4, 12) ||197 ||0.720 ||
```

Check ```_frgcombined_goodfr20_normalized.tif``` that the neurons are correctly segmented.

![frgcombined](https://github.com/tkatsuki/Flyception2R/blob/master/frgcombined.png)

