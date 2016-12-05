library(devtools)
library(flowPloidyData)
load_all()

batch1 <- batchFlowHist(flowPloidyFiles, channel="FL3.INT.LIN")
b1b <- browseFlowHist(batch1)


vac2 <-
  FlowHist(
    file = "~/research/flow/gating examples/Vac.ON.DL.14.026",
    channel = "FL2.A")

browseFlowHist(vac2)

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




