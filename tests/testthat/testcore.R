library(flowPloidy)
library(flowPloidyData)
context("core features")

fh1 <- FlowHist(file = flowPloidyFiles()["188-15.LMD"],
                channel = "FL3.INT.LIN", analyze = TRUE)

fh2 <- FlowHist(file = flowPloidyFiles()["222.LMD"],
                channel = "FL3.INT.LIN", analyze = FALSE)

fh3 <- FlowHist(file = flowPloidyFiles()["240+S.LMD"],
                channel = "FL3.INT.LIN", debris = "MC", analyze = TRUE)

test_that("FlowHist objects print without error", {
  expect_error(print(fh1), NA)
  expect_error(print(fh2), NA)
  expect_error(print(fh3), NA)
})

## Testing that analyzed, unanalyzed, and groups of possibly mixed analyzed
## and unanalyzed FlowHist objects are all handled politely by
## tabulateFlowHist:  
test_that("FlowHist objects tabulate without error", {
  expect_error(tabulateFlowHist(fh1), NA)
  expect_error(tabulateFlowHist(fh2), NA)
  expect_error(tabulateFlowHist(list(fh1, fh2), NA))
  expect_error(tabulateFlowHist(list(fh2, fh3)), NA)
  expect_error(tabulateFlowHist(list(fh1, fh2, fh3)), NA)
})

test_that("FlowHist plots work without error", {
  expect_error(plot(fh1, init = TRUE), NA)
  expect_error(plot(fh2, init = TRUE), NA)
  expect_error(plot(fh3, init = TRUE), NA)
})

test_that("Crappy FCM files can be loaded and plotted without error", {
  expect_error(fhBad <- FlowHist(fpBad(), channel = "FL2.A"), NA)
  expect_error(plot(fhBad, init = TRUE), NA)
  expect_error(print(fhBad), NA)
})
