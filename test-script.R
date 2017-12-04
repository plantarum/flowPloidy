library(devtools)
library(flowPloidyData)
load_all()

badFiles <- list.files("/home/tws/research/flow/paul/badFiles", full.names
                       = TRUE)

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




