library(flowPloidy)
library(flowPloidyData)
context("models and components")

fh1 <- FlowHist(file = flowPloidyFiles["188-15.LMD"],
                channel = "FL3.INT.LIN")
fh2 <- fh1

## ignore peak at 100, such that both A and B components have no G2 peak
fhPeaks(fh2)[1, 1] <- fhPeaks(fh2)[1, 1] * 2
fh2 <- addComponents(fh2)
fh2 <- makeModel(fh2)
fh2 <- getInit(fh2)

## replace peak at 100, A component should now have G1 and G2 peak
fh3 <- fh2

fhPeaks(fh3)[1, 1] <- fhPeaks(fh3)[1, 1] / 2
fh3 <- addComponents(fh3)
fh3 <- makeModel(fh3)
fh3 <- getInit(fh3)

test_that("Peak components are updated when selecting peaks", {
  expect_true("fA1" %in% names(fhComps(fh1)))
  expect_true("fA2" %in% names(fhComps(fh1)))
  expect_true("fB1" %in% names(fhComps(fh1)))
  expect_false("fB2" %in% names(fhComps(fh1)))

  expect_true("fA1" %in% names(fhComps(fh2)))
  expect_false("fA2" %in% names(fhComps(fh2)))
  expect_true("fB1" %in% names(fhComps(fh2)))
  expect_false("fB2" %in% names(fhComps(fh2)))

  expect_true("fA1" %in% names(fhComps(fh3)))
  expect_true("fA2" %in% names(fhComps(fh3)))
  expect_true("fB1" %in% names(fhComps(fh3)))
  expect_false("fB2" %in% names(fhComps(fh3)))
})

