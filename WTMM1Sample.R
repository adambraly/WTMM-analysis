# REGRESSION
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
library(wmtsa)
Working.Dir <- "C:/Users/delucia.TTU277174/Desktop/Adam"   
setwd     (Working.Dir)
  source  ("C:/Users/delucia.TTU277174/Desktop/Adam/functions.R")
  source  ("C:/Users/delucia.TTU277174/Desktop/Adam/mfdfa.R")
  
 
# Create File
  file.base   <- "Sampledata"                 # file name here,*.csv
    sub.file  <- "f.name"     # subfile you will be working on
    ts.header <- "RT"                # header for the timeseries in the csv 
    startBlk  <- 1
    endBlk    <- 1
  data.dir <- "/data/"                     # Insert data directory here
  output.dir  <- "/output/"                # Insert output directory here
 now <-format(Sys.time(), "%b%d%H%M%S")  # Setting a timestamp

    input.dir  <- paste(Working.Dir,output.dir,file.base, "$", 
                        sub.file, "--",now, sep="")
       dir.create(input.dir)
# Creating Filenames
  input.file  <- paste(Working.Dir,data.dir,file.base,".csv", sep="")
  output.file <- paste(input.dir,"/",file.base, "-",now, sep="")

# Load Data
require(car, quietly = TRUE)
  f.name      <- read.csv(input.file, sep=",", na.strings="3996") 
    f.name$CorrectAnswer<- recode(f.name$Correct,"FALSE=NA")
  f.name      <- na.omit(f.name)            # Omit NA
      # View (f.name)
      # str(f.name)                         # Display $frame Titles
      # summary(f.name)
# Subset Creation
  d.302  <- subset(f.name, ID=="302" ) 
 

# Plot by individual block
sub.file.p   <- eval (parse(text = sub.file))
# View(sub.file.p)

Blk <- startBlk
# Loop Start
while (Blk <= endBlk){
  sub.file.blk <- subset(sub.file.p, Block==Blk)
  par(resetPar())  
# creating timeseries
  sub.file.blk.dp <- deparse(substitute(sub.file.blk))
  sub.ts.blk <- paste(sub.file.blk.dp,"$",ts.header, sep= "")
  sub.ts.blk.p   <- eval (parse(text = sub.ts.blk))
  strLen   <-  nchar(sub.ts.blk)
  adjTitle <-  substr(sub.file, 3, strLen-4)
  finTitle <-  paste(adjTitle,"_Block", as.character(Blk), sep= "")


# Creating x
require(outliers, quietly = TRUE)
sub.ts.blk.p <- sub.ts.blk.p*1000 #converting from seconds to miliseconds
x <- (sub.ts.blk.p)
y <- rm.outlier(x)
y <- na.omit(y)
bk.sub.blk.l <- length(y)
x.mean <- mean(x)
x.sd   <- sd(x)

# MFDFA Control Panel
  x <- y
  scmin = 4
  scmax = bk.sub.blk.l
  ressc = 19
  qmin = -5
  qmax = 5
  qres = 1147

  par(mfrow = c(1, 2))
  boxplot(sub.ts.blk.p, main = "preTransform")
  boxplot(y, main = "post.outlier.Transform")
    x <- (y)
# View (x)
# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "BoxPlot", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# CWT/WTMM Creation
require(tseries, quietly = TRUE)
require(wmtsa, quietly = TRUE)
timeseries <-ts(x)  
  bk.sub.l   <- length(x)
  W        <- wavCWT(timeseries, wavelet="gaussian2")
  scales <- attr(W, "scale")  
  W.Tree   <- wavCWTTree(W, n.octave.min=1)

# Hurst Stuff 
require(fractal, quietly = TRUE)
  h.x= (timeseries)
  set.seed(bk.sub.l)
  walk <- cumsum(h.x)

# calculate the Hurst coefficient
  methods <- c("standard","smoothed")
  z.hurst <- lapply(methods, function(method, walk){
              hurstSpec(walk, method=method, sdf.method="multitaper")
              },walk=walk )
  names(z.hurst) <- methods

# MFDFA Function call      
  x.mfdfa <- mfdfa(x,scmin,scmax,ressc,qmin,qmax,qres)

# Making Graphs
  numdata <- length(W)
  pal <- colorRampPalette(c("black", "red"))
  colors <- pal(numdata)

par(mfrow = c(1, 2))
plot (W, col=colors, xlab="Time", ylab="Scale", main = "CWT")
  plot(x.mfdfa[[7]],x.mfdfa[[8]],xlim=c(-1,4.5), ylim=c(0.1,1), 
       xlab="Holder exponent", ylab="Fractal Dimension", 
       main = "Multifractal Spectrum")

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "-CWT_DynSpect-", sep="")
out.file <- paste(output.file.pre, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# WTMM Tree Plots
par(resetPar())  
par (mfrow = c(1, 1))
plot (W, col=colors, xlab="Time", ylab="Scale", main = "CWT")
  if (is.R()) plot(W.Tree, extrema=TRUE, add=TRUE)

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "-CWTWTMMoverlay-", sep="")
out.file <- paste(output.file.pre, ".png", sep="")
dev.copy(png, out.file)
dev.off()

plot (W.Tree)
# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "-WTMM-", sep="")
out.file <- paste(output.file.pre, ".png", sep="")
dev.copy(png, out.file)
dev.off()

layout(matrix(c(1,1,1,1), 2, 4, byrow = TRUE))
plot(W.Tree[range=c(0, 200)], fit=TRUE)

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "-WTMMPlots-", sep="")
out.file <- paste(output.file.pre, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# Hurst/Hist Plots
par(mfrow = c(1, 2))
  plot(z.hurst[[1]])

hist(x,breaks=5,xlim=c(0,400),freq=FALSE)
curve(dnorm(x, mean = x.mean, sd = x.sd), add = T, lwd=2)

# Output to File
  output.file.pre <- paste(output.file, "-", finTitle, "-hist_hurst-", sep="")
  out.file <- paste(output.file.pre,now, ".png", sep="")
  dev.copy(png, out.file)
  dev.off()
  
  Blk = Blk + 1
  par(resetPar())  
}
# creating timeseries
sub.ts <- paste(sub.file,"$",ts.header, sep="")
r.ID   <- eval (parse(text = sub.ts))
ts     <- r.ID
strLen <- nchar(sub.ts)
adjTitle <- substr(sub.ts, 3, strLen-4)
  timeseries <-ts(ts)
      bk.tot.l   <- length(ts)
      W      <- wavCWT(timeseries, wavelet="gaussian2")
      W.Tree <- wavCWTTree(W, n.octave.min=1)

# Making CWT Graph
par(mfrow = c(1, 1))
  numdata <- length(W)
  pal <- colorRampPalette(c("black", "green"))
  colors <- pal(numdata)
  plot (W, col=colors, main= adjTitle)
  
# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "Overall_CWT", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# CWT With Tree
plot (W, col=colors, main= adjTitle)
if (is.R()) plot(W.Tree, extrema=TRUE, add=TRUE)

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "Overall_CWTTree", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# CWT With Tree
plot (W.Tree,  main= adjTitle)

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "Overall_W.Tree", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# W.tree Log
layout(matrix(c(1,1,1,1), 2, 4, byrow = TRUE))
plot(W.Tree[range=c(0, 200)], fit=TRUE)

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "Overall_Log", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# Creating x
require(outliers, quietly = TRUE)
ts <- ts*1000 #converting from seconds to miliseconds
x <-  sqrt(ts)
y <- rm.outlier(x)
y <- na.omit(y)
tot.l <- length(y)
x.mean <- mean(x)
x.sd   <- sd(x)

# Control Panel 
x<-y
scmin = 4
scmax = tot.l
ressc = 19
qmin = -5
qmax = 5
qres = 1147

# Hurst Stuff
require(fractal, quietly = TRUE)
h.x= (x)
set.seed(tot.l)
walk <- cumsum(h.x)

# calculate the Hurst coefficient
methods <- c("standard","smoothed")
z.hurst <- lapply(methods, function(method, walk){
  hurstSpec(walk, method=method, sdf.method="multitaper")
},walk=walk )
names(z.hurst) <- methods

# Hurst/Hist Plots
par(mfrow = c(1, 2))
plot(z.hurst[[1]])

hist(x,breaks=5,xlim=c(1,600),freq=FALSE)
curve(dnorm(x, mean = x.mean, sd = x.sd), add = T, lwd=2)

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "hurst_hist", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()

# Function Call
par(mfrow = c(1, 1))
  x.mfdfa <- mfdfa(x,scmin,scmax,ressc,qmin,qmax,qres)
  plot(x.mfdfa[[7]],x.mfdfa[[8]],xlim=c(-1,4.5), ylim=c(0.1,1),main=adjTitle, 
    xlab="Holder exponent", ylab="Fractal Dimension")

# Output to File
output.file.pre <- paste(output.file, "-", finTitle, "Overall_MFSpect", sep="")
out.file <- paste(output.file.pre,now, ".png", sep="")
dev.copy(png, out.file)
dev.off()



