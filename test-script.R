library(devtools)
library(flowPloidyData)
load_all()


batch1 <-batchFlowHist(files = flowPloidyFiles, channel = "FL3.INT.LIN")
batch1b <- browseFlowHist(batch1)

thresh = 0.0002
par(mfrow = c(2,2))
deb <- batch1[["240+S.LMD"]]
debm <- batch1[["240+S.LMD"]]
debm2 <- batch1[["240+S.LMD"]]
plot(deb)
sscfl <- exprs(fhRaw(deb))[, "SS.INT.LIN"] /
  exprs(fhRaw(deb))[, "FL3.INT.LIN"]
fl <- exprs(fhRaw(deb))[, "FL3.INT.LIN"]
plot(sscfl ~ fl, col = "#44444411", pch = 16, ylim = c(0, 0.1))
dat <- as.data.frame(exprs(fhRaw(deb)))
dat$ssfl <- sscfl
abline(h = thresh, col = 2)
dat1 <- dat[dat$ssfl > thresh, ]
dat2 <- dat[dat$ssfl <= thresh, ]
debm@raw@exprs <- as.matrix(dat1)
debm2@raw@exprs <- as.matrix(dat2)
debm <- setBins(debm)
debm2 <- setBins(debm2)
plot(debm)
plot(debm2)




#######################
## Linearity testing ##
#######################

linFiles <- list.files("~/research/flow/linearity_data/", "LMD",
                       full.names = TRUE)
lintest <- batchFlowHist(linFiles, channel = "FL3.INT.LIN")
lintest <- batchFlowHist(linFiles, channel = "FL3.INT.LIN",
                         linearity = "variable")

lintestv <- browseFlowHist(lintest)

fh1m <-FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN",
                 analyze = TRUE)
fh1s <- updateFlowHist(fh1m, opts = list("SC"), analyze = TRUE)
fh1sl<- updateFlowHist(fh1s, linearity = "variable", analyze = TRUE)
fh1ml<- updateFlowHist(fh1mc, linearity = "variable", analyze = TRUE)

plot(fh1m)
dev.new()
##dev.set()
plot(fh1s)
dev.new()
plot(fh1ml)
dev.new()
##dev.set()
plot(fh1sl)


fh1mcl <-FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN",
                  linearity = "variable")
fh1sc <-FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN",
                 opts = list("SC"))

plotFH(fh1mcl)
plot(fh1mcl, init = TRUE)

fh1sc <- fhAnalyze(fh1sc)
plot(fh1sc)

fh1mc <- fhAnalyze(fh1mc)
plot(fh1mc)

fh1mcl <- fhAnalyze(fh1mcl)
plot(fh1mcl)

batch1 <-batchFlowHist(files = flowPloidyFiles, channel = "FL3.INT.LIN")
batch1b <- browseFlowHist(batch1)

batch2 <-batchFlowHist(files = flowPloidyFiles, channel = "FL3.INT.LIN")

tabulateFlowHist(batch1)

tmp <- browseFlowHist(batch1)

fh4S4 <-FlowHist(file = flowPloidyFiles[4], channel = "FL3.INT.LIN")
plot(fh4S4, init = TRUE)
fh4S4 <- pickInit(fh4S4)

tmp <- list()

for(i in fh1S4@comps){
  tmp <- c(tmp, i@initParams(fh1S4))
}
  
fA1@initParams(fh1S4)

fh1S4b <- setBins(fh1S4, 512)

library(devtools)
load_all()
chan = "FL3.INT.LIN"
files <- list.files(system.file("extdata", package = "flowPloidy"),
                    full.names = TRUE)

i <- 1
#filei <- system.file("extdata", files[i], package = "flowPloidy")
filei <- files[i]
fhi <- flowHist(FILE = filei, channel = chan)
plot(fhi, init = TRUE)
fhi <- pickInit(fhi)
fhi <- fhAnalyze(fhi)
plot(fhi)
fhi

my.files <- list.files(system.file("extdata/", package = "flowPloidy"),
                       pattern = "*.LMD", full.names = TRUE)

batch1 <- histBatch(my.files, channel = "FL3.INT.LIN", window = 20,
                    smooth = 20)

batch1 <- flowShiny(batch1)


parOld <- par(ask = TRUE)
lapply(batch1, FUN = plot)
## press enter to scroll through your files!
par(parOld)

library(pander)
myReport <- Pandoc$new("Tyler Smith", "flowPloidy Test")
myReport$format <- "html"
res <- list()
for(i in seq_along(files)){
  message("processing ", files[i])
  list[[files[i]]] <- NULL
  filei <- system.file("extdata", files[i], package = "flowPloidy")
  fhi <- flowHist(FILE = filei, CHANNEL = chan)
  try(list[[files[[files[i]]] <- fhAnalyze(fhi))
  myReport$add.paragraph(paste("#", files[i]))
  myReport$add(exportFlowHist(fhi))
  myReport$add(plot(fhi))
}              

histReport <- function(hb, author, title, reportFile = NULL,
                       verbose = TRUE){
  if(is.null(reportFile))
    reportFile <- tempfile(tmpdir = getwd())
  myReport <- Pandoc$new(author, title)
  for(i in seq_along(hb)){
    myReport$add(plot(hb[[i]], init = TRUE))
  }

  myReport$format <- "html"
  myReport$export(reportFile, options = " ")

}
  
histProcess <- function(files, chan, author, title,
                        dataFile = NULL, reportFile = NULL, bins = 256,
                        verbose = TRUE){ 
  if(is.null(reportFile))
    reportFile <- tempfile(tmpdir = getwd())
  myReport <- Pandoc$new(author, title)
  res <- list()

  for(i in seq_along(files)){
    if(verbose) message("processing ", files[i])
    #filei <- system.file("extdata", files[i], package = "flowPloidy")
    res[[files[i]]] <- flowHist(FILE = files[i], CHANNEL = chan)
    tryVal <- try(res[[files[i]]] <- fhAnalyze(res[[files[i]]]))
    if(verbose && inherits(tryVal, "try-error")) message("-- analysis failed")
    myReport$add(plot(res[[files[i]]], init = TRUE))
  }              

  exportFlowHist(res, file = dataFile)
  myReport$format <- "html"
  myReport$export(reportFile, options = " ")

  return(res)
}
  
## 12: 734.LMD requires manually setting the init values
## 8: "337.LMD" requires manually setting the init values
## "240+S.LMD" requires manually setting the init values
## "248+r.LMD" check linearity

file1 <- system.file("extdata", "188-15.LMD", package = "flowPloidy")
fh1 <- flowHist(FILE = file1, CHANNEL = chan)
plot(fh1, init = TRUE)
fh1 <- fhAnalyze(fh1)
plot(fh1)
fh1

file2 <- system.file("extdata", "SM239.LMD", package = "flowPloidy")
fh2c <- flowHist(FILE = file2, channel = chan)
plot(fh2, init = TRUE)
fh2 <- fhAnalyze(fh2)
fh2b <- fhAnalyze(fh2b)
plot(fh2)
plot(fh2b)
fh2

file3 <- system.file("extdata", "226.LMD", package = "flowPloidy")
fh3 <- flowHist(FILE = file3, CHANNEL = chan)
plot(fh3, init = TRUE)
fh3 <- fhAnalyze(fh3)
plot(fh3)
fh3


222.LMD

