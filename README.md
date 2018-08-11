# FlyceptionR
R scripts and utilities for analyzing Flyception data

## Requirement
Windows: R and Rtools

Mac: R and Xcode

ImageJ https://imagej.nih.gov/ij/

Set ImageJ (and Fiji if present) to save Tiff in Intel byte order from Edit -> Options -> Input/Output... 

## Installation
The following commands will automatically install packages necessary for running FlyceptionR.

```
install.packages(c("devtools", "ggplot2", "RNiftyReg", "zoo", "loggit"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocInstaller", "EBImage"))
library(devtools)
devtools::install_github("tkatsuki/Flyception2R")
library(Flyception2R)
```

## Usage example
Normally, you just need to specify which frames of a video file you want to analyze.

```
dir <- "/Volumes/LaCie/P1_GCaMP6s_tdTomato_06182018_CW_Dual_Laser/P1-Gal4_UAS-GCaMP6s_tdTomato_12/"
Flyception2R(dir=dir, FOI=c(4242, 4556), interaction=F, flash=1)
```

While the code is running you will be asked if the window_size and window_offset are acceptable. Check the file _redwindow.tif with ImageJ, for example, and adjust the window size or offset so that the neurons of your interest are in the window.

```
[1] "Current window_size is x=68 y=28"
[1] "Current window_offset is x=-4 y=12"
Check redwindow.tif. Is the window size good (Y or N)?:Y
Check redwindow.tif. Is the window offset good (Y or N)?:N
Enter new x offset:0
Enter new y offset:0
```

The function spits out results like so;

```
window_size was x=68 y=28
window_offset was x=-4 y=12
FOI was from 4242 to 4556
Max F_ratio intensity in this bout was 0.719705128074784
Number of good frames was 197
||c(4242, 4556) ||c(68, 28) ||c(-4, 12) ||197 ||0.720 ||
```