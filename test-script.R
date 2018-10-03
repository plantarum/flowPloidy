library(devtools)
library(flowPloidyData)
load_all()

batch1 <- batchFlowHist(flowPloidyFiles, channel="FL3.INT.LIN",
                        standards = 2.5)

batch1 <- batchFlowHist(flowPloidyFiles, channel="FL3.INT.LIN")


endo <- FlowHist(file = "~/research/flow/problems/endopolyploidy.LMD",
                channel = "FL3.INT.LIN", g2 = FALSE)

etFiles <- list.files("~/research/flow/problems/etienne/",
                      pattern = ".fcs", full.names = TRUE) 

viewFlowChannels(etFiles[1])

etienne <- batchFlowHist(etFiles, channel = "FL3.A", standard = 2.5)

et  <- browseFlowHist(etienne)

etienne <- FlowHist(file = "~/research/flow/problems/etienne/",
                channel = "FL3.INT.LIN", g2 = FALSE)

b1 <- batch1[[1]]

b1Br <- browseFlowHist(batch1)

b3 <- b1Br[[3]]

nuria <- FlowHist(file = "~/research/flow/nuria.fcs", channel = "DAPI.A",
                  g2 = FALSE, debrisLimit = 10, samples = 5)


fhWMean <- function(fh, range){
  sum(fhHistData(fh)$intensity[range] * range) /
    sum(fhHistData(fh)$intensity[range])
}

fhVisual <- function(fh, range){
  dat <- rep(range, times = fhHistData(fh)$intensity[range])
  message("Mean: ", mean(dat))
  message("CV: ", sd(dat)/ mean(dat))
  message("counts: ", length(dat))
}

svg(filename = "visual_gate.svg", width = 25/2.54, height = 25*0.66/2.54,
    family = "lato")
par(mar = c(5,5, 1, 1))
plot(b1Br[[11]], nls = FALSE, init = FALSE, comps= FALSE, main = "",
     xlim = c(0, 175))
rcol = "#FF000060"
rect(xleft = 42, xright = 47, ybottom = 0, ytop=300, col = rcol,
     border = rcol)
rect(xleft = 55, xright = 60, ybottom = 0, ytop=300, col = rcol,
     border = rcol)
rect(xleft = 125, xright = 133, ybottom = 0, ytop=300, col = rcol,
     border = rcol)
rect(xleft = 144, xright = 152, ybottom = 0, ytop=300, col = rcol,
     border = rcol)
dev.off()


svg(filename = "nls_fit.svg", width = 25/2.54, height = 25*0.66/2.54,
    family = "lato")
#dev.new(width = 25/2.54, height = 25 * 0.66/2.54)
par(mar = c(5,5, 1, 1))
plot(b1Br[[11]], nls = TRUE, init = FALSE, comps= TRUE, main = "",
     xlim = c(0, 175))
dev.off()


svg(filename = "endo.svg", width = 25/2.54, height = 25*0.66/2.54,
    family = "lato")
#dev.new(width = 25/2.54, height = 25 * 0.66/2.54)
par(mar = c(5,5, 1, 1))
plot(endo, nls = TRUE, init = FALSE, comps= TRUE, main = "",
     xlim = c(0, 125))
dev.off()


abline(v = 42)
abline(v = 47)

abline(v = 55)
abline(v = 60)

fhVisual(b1Br[[11]], 42:60)
fhVisual(b1Br[[11]], 47:55)
fhVisual(b1Br[[11]], 42:55)
fhVisual(b1Br[[11]], 47:60)
## mean range from 49.8 to 51.5
## CV range from 7.3% to 4.0%
## counts range from 1883 to 1475

abline(v = 125)
abline(v = 133)

abline(v = 144)
abline(v = 152)

fhVisual(b1Br[[11]], 125:152)
fhVisual(b1Br[[11]], 133:144)
fhVisual(b1Br[[11]], 125:144)
fhVisual(b1Br[[11]], 133:152)
## mean range from 137.4 to 139.4
## CV range from 2.1% to 3.6%
## counts range from 1664 to 2126
## ratio ranges from 0.357 to 0.37
## with a standard size of 3 pg, this translates to a range of 0.05 pg,
## introduce 5% or more variation due to subjectivity




browseFlowHist(nuria)


chh <- FlowHist(file = "~/research/flow/problems/valentin/chh101.fcs",
                channel = "FL4.A", standards = 1.99)

chc <- FlowHist(file = "~/research/flow/problems/valentin/chc203.fcs",
                channel = "FL4.A", standards = 1.99, debrisLimit = 0)


viewFlowChannels("~/research/flow/sample.001.fcs")

brian <- FlowHist("~/research/flow/problems/brian/sample.001.fcs",
                  channel = "FL2.A", samples = 1)

library(digest)
digest(brian@histData$intensity, "md5")
# [1] "f0724874dd1beeb6bb5e20fac8e3de5e"
brian@raw@description$`$DATE`
# [1] "17-Jan-18"
brian@raw@description$`$BTIM`
# [1] "17:11:47"
brian@raw@description$`FILE GUID`
# [1] "8E581367-FBEC-11E7-AE75-001124882388"

lowPeaks <- list.files("/home/tws/research/flow/paul/low-peaks",
                       full.names = TRUE)

low <- batchFlowHist(lowPeaks, channel= "FL2.A", debrisLimit = 10)

low <- browseFlowHist(low)

badFiles <- list.files("/home/tws/research/flow/paul/badFiles", full.names
                       = TRUE)
fpBad <- badFiles[8]
bad <- batchFlowHist(badFiles, channel = "FL2.A")

bad <- browseFlowHist(bad)


etienne <- batchFlowHist(list.files("/home/tws/research/flow/etienne/",
                                    full.names = TRUE),
                         channel = "FL3.INT.LIN", standards = 2.5) 

etBr <- browseFlowHist(etienne)

e1 <- FlowHist("Calli+laxa-moved gate 2017-02-07 553.LMD",
               channel = "FL3.INT.LIN")

e1 <- FlowHist("Calliciprus+laxa 2017-02-07 550.LMD",
               channel = "FL3.INT.LIN")

e1 <- FlowHist("call1.LMD", channel = "FL3.INT.LIN")


viewFlowChannels("~/research/flow/160930_HCT116_Ca-new.fcs")
viewFlowChannels("~/research/flow/160930_HCT116_Ca.fcs")

hct <- FlowHist("~/research/flow/160930_HCT116_Ca-new.fcs",
                channel = "Propidium.Iodide.A")

hctB <- browseFlowHist(hct)

hctNoDebris <- dropComponents(hct, c("SC", "AG"))
hctNoDebris <- fhAnalyze(hctNoDebris)

plot(hctNoDebris, ylim = c(0, 800))

browseFlowHist(hct)
browseFlowHist(hctNoDebris)
plot(hct, ylim = c(0, 800))

hct <- FlowHist("~/research/flow/160930_HCT116_Ca-new.fcs",
                channel = "Propidium.Iodide.A", bins = 512)


batch1 <- batchFlowHist(flowPloidyFiles, channel="FL3.INT.LIN",
                        standards = 2.5)
b1 <- batch1[[1]]

b1Br <- browseFlowHist(batch1)

b1 <- batch1[[1]]

do.call(integrate, c(substitute(fhModel(b1)), as.list(coef(fhNLS(b1))),
                     lower = 0, upper = 256, subdivisions = 512,
                     SCvals = substitute(fhHistData(b1)["SCvals"]),
                     DBvals = substitute(fhHistData(b1)["DBvals"]),
                     TRvals = substitute(fhHistData(b1)["TRvals"]),
                     QDvals = substitute(fhHistData(b1)["QDvals"]))) 

integrate(fhModel(b1), 
          a1 = coef(fhNLS(b1))["a1"],
          Ma = coef(fhNLS(b1))["Ma"],
          Sa = coef(fhNLS(b1))["Sa"],
          a2 = coef(fhNLS(b1))["a2"],
          d = coef(fhNLS(b1))["d"],
          BRA = coef(fhNLS(b1))["BRA"],
          b1 = coef(fhNLS(b1))["b1"],
          Mb = coef(fhNLS(b1))["Mb"],
          Sb = coef(fhNLS(b1))["Sb"],
          SCa = coef(fhNLS(b1))["SCa"],
          SCvals = fhHistData(b1)["SCvals"], 
          Ap = coef(fhNLS(b1))["Ap"],
          DBvals = fhHistData(b1)["DBvals"],
          TRvals = fhHistData(b1)["TRvals"],
          QDvals = fhHistData(b1)["QDvals"],
          lower = 0, upper = 256, subdivisions = 1000)

b1b <- browseFlowHist(batch1)


vac2 <-
  FlowHist(
    file = "~/research/flow/gating examples/Vac.ON.DL.14.026",
    channel = "FL2.A")

browseFlowHist(vac2)

vac2ND <- dropComponents(vac2, c("SC", "AG"))
vac2NDAn <- browseFlowHist(vac2ND)

dat <- exprs(fhRaw(vac2))
thresh <- 0.05
test <- dat[ , "FL3.H"]/ dat[, "FL2.A"] < thresh

vac2g <- setGate(vac2, test)
## vac2g <- cleanPeaks(findPeaks(vac2g))
vac2g <- pickPeaks(vac2g)
vac2g <- addComponents(vac2g)
vac2g <- setLimits(vac2g)
vac2g <- makeModel(vac2g)
vac2g <- getInit(vac2g)
vac2g <- fhAnalyze(vac2g)


par(mfrow=c(2,2))
plot(vac2, sub = "ungated")
plot(dat[, "FL2.A"], dat[ , "FL3.H"]/ dat[, "FL2.A"], pch = 16,
     ylim = c(0, 0.08),
     col = "#11111130", cex = 0.5, main = fhFile(vac2))
abline(h = thresh, col = 2)
plot(vac2g, sub = "gated")
plotResid(vac2g)


## SSC plots
thresh <- 1
test <- dat[ , "SSC.H"]/ dat[, "FL2.A"] < thresh

vac2g <- setGate(vac2, test)
## vac2g <- cleanPeaks(findPeaks(vac2g))
vac2g <- pickPeaks(vac2g)
vac2g <- addComponents(vac2g)
vac2g <- setLimits(vac2g)
vac2g <- makeModel(vac2g)
vac2g <- getInit(vac2g)
vac2g <- fhAnalyze(vac2g)


par(mfrow=c(2,2))
plot(vac2, sub = "ungated")
plot(dat[, "FL2.A"], dat[ , "SSC.H"]/ dat[, "FL2.A"], pch = 16,
     ylim = c(0, 3),
     col = "#11111130", cex = 0.5, main = fhFile(vac2))
abline(h = thresh, col = 2)
plot(vac2g, sub = "gated")
plotResid(vac2g)



vac2g <- fhAnalyze(vac2g)

vac2gd <- dropComponents(vac2g, c("SC", "MC"))
vac2gd <- fhAnalyze(vac2gd)

vac2g <- findPeaks(vac2g)
vac2g <- cleanPeaks(vac2g, window = 20)
vac2g <- addComponents(vac2g)
vac2g <- makeModel(vac2g)
vac2g <- getInit(vac2g)
vac2g <- pickInit(vac2g)

plot(vac2g, init = TRUE)


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



