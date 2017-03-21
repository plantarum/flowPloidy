## Functions for analyzing FlowHist datasets

#' @importFrom car deltaMethod
NULL

#' @importFrom minpack.lm nlsLM
NULL

#' Complete non-linear regression analysis of FlowHist histogram data
#'
#' Completes the NLS analysis, and calculates the modelled events and CVs
#' for the result.
#' 
#' @title fhAnalyze
#' @param fh a \code{\link{FlowHist}} object
#' @return a \code{\link{FlowHist}} object with the analysis (nls, counts, cv,
#'   RCS) slots filled.
#' @seealso \code{\link{FlowHist}}
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData)
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1 <- fhAnalyze(fh1)
#' @export
fhAnalyze <- function(fh){
  message("analyzing ", fhFile(fh))
  tryVal <- try(fh <- fhDoNLS(fh), silent = TRUE)
  if(inherits(tryVal, "try-error")){
    message("    *** Model Fit Needs Attention: ", fhFile(fh), " ***")
  } else {
    fh <- fhDoCounts(fh)
    fh <- fhDoCV(fh)
    fh <- fhDoRCS(fh)
  }
  return(fh)
}

fhDoNLS <- function(fh){
  model <- fhModel(fh)
  form1 <- paste("intensity ~ model(")
  args <- fhArgs(fh)
  pLims <- fhLimits(fh)
  pLims <- pLims[! names(pLims) %in% fhSpecialParams(fh)]
  lLims <- sapply(pLims, function(x) x[1])
  uLims <- sapply(pLims, function(x) x[2])
  args <- paste(args, collapse = ", ")
  form3 <- paste(", ", getSpecialParamArgs(fh), ")")
  form <- as.formula(paste(form1, args, form3))

  ## Ignore the lowest channels, before we start modeling the debris
  ## component.
  dat <- fhHistData(fh)
  start <- fhStart(dat$intensity)
  dat <- dat[-(1:(start - 1)), ]

  fhNLS(fh) <- nlsLM(formula = form, start = fhInit(fh),
                     data = dat, 
                     lower = lLims, upper = uLims,
                     control = list(ftol = .Machine$double.xmin,
                                    ptol = .Machine$double.xmin,
                                    maxiter = 1024))
  return(fh)
}
  
fhDoCounts <- function(fh){
  ## lower was originally an argument to fhCount, but I don't think it will
  ## ever be anything other than 0?
  lower = 0
  ## similarly, upper was an argument, but it should always be the number
  ## of bins 
  upper = nrow(fhHistData(fh))
  ## I think anything >= the number of bins should be fine for
  ## subdivisions: 
  subdivisions = upper * 2

  ## This integration is the single slowest step in the analysis, requiring
  ## 10-40 seconds or more to complete. For usability, I'm no longer
  ## supporting it. Perhaps there are ways to speed it up, but all the
  ## important count data is still available quickly.
  ## total <-
  ##   do.call(integrate,
  ##           c(substitute(fh$model),
  ##             as.list(coef(fh$nls)),
  ##             SCvals = substitute(fh$data$SCvals),
  ##             lower = lower, upper = upper,
  ##             subdivisions = subdivisions))
  firstPeak <-
    integrate(mcFunc(fhComps(fh)$fA1), a1 = coef(fhNLS(fh))["a1"],
              Ma = coef(fhNLS(fh))["Ma"],
              Sa = coef(fhNLS(fh))["Sa"],
              lower = lower, upper = upper,
              subdivisions = 1000)
  if("fB1" %in% names(fhComps(fh))){
    secondPeak <-
      integrate(mcFunc(fhComps(fh)$fB1), b1 = coef(fhNLS(fh))["b1"],
                Mb = coef(fhNLS(fh))["Mb"],
                Sb = coef(fhNLS(fh))["Sb"],
                lower = lower, upper = upper,
                subdivisions = 1000)
  } else {
    secondPeak <- NULL
  }

  if("fC1" %in% names(fhComps(fh))){
    thirdPeak <-
      integrate(mcFunc(fhComps(fh)$fC1), c1 = coef(fhNLS(fh))["c1"],
                Mc = coef(fhNLS(fh))["Mc"],
                Sc = coef(fhNLS(fh))["Sc"],
                lower = lower, upper = upper,
                subdivisions = 1000)
  } else {
    thirdPeak <- NULL
  }

  fhCounts(fh) <- list(firstPeak = firstPeak, secondPeak = secondPeak,
                       thirdPeak = thirdPeak)

  fh
}  

fhDoCV <- function(fh){
  CVa <- coef(fhNLS(fh))["Sa"]/coef(fhNLS(fh))["Ma"]
  if("fB1" %in% names(fhComps(fh))){
    CVb <- coef(fhNLS(fh))["Sb"]/coef(fhNLS(fh))["Mb"]
    AB <- deltaMethod(fhNLS(fh), "Ma/Mb")
  } else {
    CVb <- CI <- AB <- NULL
  }

  if("fC1" %in% names(fhComps(fh))){
    CVc <- coef(fhNLS(fh))["Sc"]/coef(fhNLS(fh))["Mc"]
    AC <- deltaMethod(fhNLS(fh), "Ma/Mc")
    BC <- deltaMethod(fhNLS(fh), "Mb/Mc")
  } else {
    CVc <- AC <- BC <- NULL
  }
  
  fhCV(fh) <- list(CVa = CVa, CVb = CVb, CVc = CVc,
                   AB = AB, AC = AC, BC = BC)
  fh
}

fhDoRCS <- function(fh){
  #########################################################################
  ## This may not be a useful measure of analysis quality. It is heavily ##
  ## influenced by the highest channels, where the expected value is     ##
  ## close to zero. Consequently, observing 2 or 3 stray events in one   ##
  ## of these channels, where the expected value may be < 0.02, produces ##
  ## a higher value of (obs - exp)^2 / exp than larger absolute          ##
  ## differences in the main region of the histogram.                    ##
  ##                                                                     ##
  ## Rabinovitch 1994:                                                   ##
  ##                                                                     ##
  ## The x2 is affected by a large number of variables, not all related  ##
  ## to goodness of the fit; these include the number of cells acquired  ##
  ## in the histogram and the end points of the analysis region used     ##
  ## within the histogram.                                               ##
  #########################################################################

  ## Ignoring the lowest channels, before the debris component starts. We
  ## calculate RCS based on the number of channels fit in the model, not
  ## the full data set, which includes a number of empty/unmodelled
  ## channels at the beginning.
  dat <- fhHistData(fh)
  start <- fhStart(fhHistData(fh)$intensity)
  dat <- dat[-(1:(start - 1)), ]
  obs <- dat$intensity

  exp <- predict(fhNLS(fh))
  chi <- sum(((obs - exp)^2) / exp)

  fhRCS(fh) <- chi/summary(fhNLS(fh))$df[2]
  fh
}
