# Functions for creating and viewing FlowHist objects.

#' @importFrom flowCore read.FCS exprs pData parameters
NULL

#' @importFrom knitr kable
NULL

#' @importFrom rmarkdown render
NULL

#' @importFrom graphics hist lines locator plot points polygon grconvertX
#'   grconvertY text abline par
NULL

#' @importFrom stats as.formula coef integrate predict pnorm
NULL

#' @importFrom utils write.table
NULL

#' @importFrom methods setClass setMethod new callNextMethod
NULL

setOldClass("nls")

#' FlowHist
#'
#' Creates a \code{\link{FlowHist}} object from an FCS file, setting up the
#' histogram data for analysis.
#'
#' Starting with a \code{\link{flowFrame}} object, read from a FCS file,
#' \code{\link{FlowHist}} will:
#'
#' \enumerate{
#' \item Extract the intensity data from \code{channel}.
#'
#' \item Remove the top bin, which contains off-scale readings we ignore
#' in the analysis.
#'
#' \item aggregates the raw data into the desired number of bins, as
#' specified with the \code{bins} argument. The default is 256, but you may
#' also try 128 or 512. Any integer is technically acceptable, but I
#' wouldn't stray from the default without a good reason.
#'
#' \item identify model components to include. All \code{\link{FlowHist}}
#' objects will have the single-cut debris model and the G1 peak for sample
#' A, and the broadened rectangle for the S-phase of sample A. Depending on
#' the data, additional components for the G2 peak and sample B (G1, G2,
#' s-phase) may also be added.
#' 
#' \item Build the NLS model. All the components are combined into a single
#' model. 
#'
#' \item Identify starting values for Gaussian (G1 and G2 peaks) model
#' components. For reasonably clean data, the built-in peak detection is
#' fine. You can evaluate this by plotting the \code{\link{FlowHist}}
#' object with the argument \code{init = TRUE}. If it doesn't look good,
#' you can play with the \code{window} and \code{smooth} arguments (which
#' is tedious!), or pick the peaks visually yourself with \code{pick =
#' TRUE}. }
#' 
#' @name FlowHist
#'
#' @param file character, the name of the file to load
#' @param files character, a vector of file names to load
#' @param channel character, the name of the data column to use
#' @param bins integer, the number of bins to use to aggregate events into
#'   a histogram
#' @param gate logical, a vector of events to exclude from analysis. In
#'   normal usage this will be set interactively, not as a function
#'   argument.
#' @param linearity character, either "variable", the default, or "fixed".
#'   If "fixed", linearity is fixed at 2; if "variable", linearity is fit
#'   as a model parameter.
#' @param debris character, either "SC", the default, or "MC", to set the
#'   debris model component to the Single-Cut or Multi-Cut models.
#' @param analyze boolean, if TRUE the  model will be analyzed
#'   immediately
#' @param opts list, currently not used, but maybe in future as a way to
#'   test additional model options
#' @param window the width of the moving window used to identify local
#'   maxima for peak detection via \code{\link{caTools::runmax}}
#' @param smooth the width of the moving window used to reduce noise in
#'   the histogram via \code{\link{caTools::runmean}}
#' @param pick boolean; if TRUE, the user will be prompted to select peaks
#'   to use for starting values. Otherwise (the default), starting values
#'   will be detected automatically.
#' @param samples integer; the number of samples in the data. Default is 2
#' (unknown and standard), but can be set to 3 if two standards are  used.
#' @param verbose boolean; if TRUE, \code{\link{histBatch}} will list files
#'   as it processes them. 
#' 
#' @slot raw a flowFrame object containing the raw data from the FCS file
#' @slot channel character, the name of the data column to use
#' @slot bins integer, the number of bins to use to aggregate events into a
#'   histogram
#' @slot gate logical, a vector indicating events to exclude from the
#'   analysis. In normal use, the gate will be modified via interactive
#'   functions, not set directly by users. 
#' @slot histdata data.frame, the columns are the histogram bin number
#'   (xx), florescence intensity (intensity), and the raw single-cut debris
#'   model values (SCVals, used in model fitting). Additional columns may
#'   be added if/when I add gating, so refer to columns by name, not
#'   position.
#' @slot peaks matrix, containing the coordinates used for peaks when
#'   calculcating initial parameter values.
#' @slot comps a list of \code{\link{modelComponent}} objects included for
#'   these data.
#' @slot model the function (built from \code{comps}) to fit to these data.
#' @slot init a list of initial parameter estimates to use in fitting the
#'   model.
#' @slot nls the nls object produced by the model fitting
#' @slot counts a list of cells counted in each peak of the fitted model
#' @slot CV a list of the coefficients of variation for each peak in the
#'   fitted model.
#' @slot RCS numeric, the residual chi-square for the fitted model.
#' @slot samples numeric, the number of samples included in the data. The
#'   default is 2 (i.e., unknown and standard), but if two standards are
#'   used it should be set to 3.
#' @return \code{\link{FlowHist}} returns a \code{\link{FlowHist}} object.
#' @author Tyler Smith
setClass(
  Class = "FlowHist",
  representation = representation(
    raw = "flowFrame", ## raw data, object defined in flowCore
    channel = "character", ## data channel to use for histogram
    samples = "integer", ## (maximum) number of sample peaks to fit
    bins = "integer", ## the number of bins to use
    linearity = "character", ## "fixed" or "variable", to determine whether
    ## or not linearity is fixed at 2, or allowed to vary as a model
    ## parameter 
    debris = "character", ## "SC" or "MC", to set the debris model. 
    gate = "logical", ## vector indicating which events to exclude from
    ## analysis, i.e., the gate
    histData = "data.frame", ## binned histogram data
    peaks = "matrix", ## peak coordinates for initial values
    opts = "list",    ## flags for selecting model components
    comps = "list", ## model components
    model = "function", ## model to fit
    init = "list", ## inital parameter estimates
    limits = "list", ## parameter limits
    nls = "nls", ## nls output
    counts = "list", ## cell counts in each peak
    CV = "list", ## CVs
    RCS = "numeric" ## residual chi-square
  ),
  prototype = prototype(
    ## TODO complete this?
  )
)

setMethod(
  f = "initialize",
  signature = "FlowHist",
  definition = function(.Object, file, channel, bins = 256,
                        window = 20, smooth = 20, pick = FALSE,
                        linearity = "variable", debris = "SC",
                        gate = logical(), samples = 2, opts = list(),
                        ...){ 
    .Object@raw <- read.FCS(file, dataset = 1, alter.names = TRUE)
    .Object@channel <- channel
    .Object@gate <- gate
    .Object@samples <- as.integer(samples)
    .Object <- setBins(.Object, bins)
    if(pick){
      .Object <- pickPeaks(.Object)
    } else {
      .Object <- findPeaks(.Object, window = window,
                           smooth = smooth)
      .Object <- cleanPeaks(.Object, window = window)
    }
    .Object@linearity <- linearity
    .Object@debris <- debris
    .Object@opts <- opts
    .Object <- addComponents(.Object)
    .Object <- setLimits(.Object)
    .Object <- makeModel(.Object)
    .Object <- getInit(.Object)
    callNextMethod(.Object, ...)
  })

###############
## Accessors ##
###############

fhGate <- function(fh){
  fh@gate
}

`fhGate<-` <- function(fh, value){
  fh@gate <- value
  fh
}

fhLimits <- function(fh){
  fh@limits
}

`fhLimits<-` <- function(fh, value){
  fh@limits <- value
  fh
}

fhSamples <- function(fh){
  fh@samples
}

`fhSamples<-` <- function(fh, value){
  fh@samples <- value
  fh
}

fhPeaks <- function(fh){
  fh@peaks
}

`fhPeaks<-` <- function(fh, value){
  fh@peaks <- value
  fh
}

fhInit <- function(fh){
  fh@init
}

`fhInit<-` <- function(fh, value){
  fh@init <- value
  fh
}

fhComps <- function(fh){
  fh@comps
}

`fhComps<-` <- function(fh, value){
  fh@comps <- value
  fh
}

fhModel <- function(fh){
  fh@model
}

`fhModel<-` <- function(fh, value){
  fh@model <- value
  fh
}

fhSpecialParams <- function(fh){
  names(getSpecialParams(fh))
}

fhArgs <- function(fh){
  res <- names(formals(fhModel(fh)))
  res <- res[!res %in% fhSpecialParams(fh)]
  res
}

fhNLS <- function(fh){
  fh@nls
}

`fhNLS<-` <- function(fh, value){
  fh@nls <- value
  fh
}

fhCounts <- function(fh){
  fh@counts
}

`fhCounts<-` <- function(fh, value){
  fh@counts <- value
  fh
}

fhCV <- function(fh){
  fh@CV
}

`fhCV<-` <- function(fh, value){
  fh@CV <- value
  fh
}

fhRCS <- function(fh){
  fh@RCS
}

`fhRCS<-` <- function(fh, value){
  fh@RCS <- value
  fh
}

fhFile <- function(fh){
  fh@raw@description$GUID
}

`fhFile<-` <- function(fh, value){
  warning("Changing the FlowHist Data File is a terrible idea")
  warning("-- it would obscure the link between your analysis")
  warning("and the raw data file, and you don't want to do that, do you?")

  fh
}

fhChannel <- function(fh){
  fh@channel
}

`fhChannel<-` <- function(fh, value){
  fh@channel <- value
  fh
}

fhBins <- function(fh){
  fh@bins
}

`fhBins<-` <- function(fh, value){
  fh@bins <- value
  fh
}

fhLinearity <- function(fh){
  fh@linearity
}

`fhLinearity<-` <- function(fh, value){
  fh@linearity <- value
  fh
}

fhDebris <- function(fh){
  fh@debris
}

`fhDebris<-` <- function(fh, value){
  fh@debris <- value
  fh
}

fhHistData <- function(fh){
  fh@histData
}

`fhHistData<-` <- function(fh, value){
  fh@histData <- value
  fh
}

fhRaw <- function(fh){
  fh@raw
}

`fhRaw<-` <- function(fh, value){
  warning("Changing the raw data is a terrible idea")
  warning("It is not supported")
  fh
}

resetFlowHist <- function(fh, from = "peaks"){
  ## Clear analysis slots
  ## Default is to clear everything from peaks onwards
  removeFrom <- c("gate", "peaks", "comps", "limits")

  ## coded to allow for further refinement, if/when additions to the
  ## FlowHist class makes it sensible to change the granularity of slot
  ## resetting. 
  
  removeNum <- which(removeFrom == from)

  rmF <- function(x)
    removeNum <= which(removeFrom == x)

  if(rmF("gate")){
    fhGate(fh) <- logical()
    fh <- setBins(fh)                   # recreate histData
  }
  if(rmF("peaks"))
    fhPeaks(fh) <- matrix()
  if(rmF("comps")){
    fhComps(fh) <- list()
  }
  if(rmF("limits")){    
    fhLimits(fh) <- list()
    fhModel(fh) <- function(){}
    fhInit(fh) <- list()
    fhNLS(fh) <- structure(list(), class = "nls")
    fhCounts(fh) <- list()
    fhCV(fh) <- list()
    fhRCS(fh) <- NA_real_
  }
  fh
}


#' @rdname FlowHist
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1
#' @export
FlowHist <- function(file, channel, bins = 256, window = 20, smooth = 20,
                     pick = FALSE, linearity = "variable", debris = "SC",
                     opts = list(), samples = 2, gate = logical(),
                     analyze = TRUE){
  fh <-  new("FlowHist", file = file, channel = channel,
             bins = as.integer(bins), window = window, smooth = smooth,
             pick = pick, linearity = linearity, debris = debris,
             opts = opts, samples = samples, gate = gate)
  if(analyze)
    fh <- fhAnalyze(fh)
  return(fh)
}

#' Displays the column names present in an FCS file
#'
#' A convenience function for viewing column names in a FCS data file, or a
#' FlowHist object. Used to select one for the \code{channel} argument
#' in \code{\link{FlowHist}}, or for viewing additional channels for use in
#' gating.
#' 
#' @title viewFlowChannels
#' @param file character, the name of a FCS data file; or the name of a
#'   FlowHist object.
#' @return A vector of column names from the FCS file/FlowHist object.
#' @seealso \code{\link{FlowHist}}
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' viewFlowChannels(flowPloidyFiles[1])
#' @export
viewFlowChannels <- function(file){
  if(class(file) == "FlowHist"){
    res <- colnames(exprs(fhRaw(file)))
  } else {
    tmp <- read.FCS(file, alter.names = TRUE, dataset = 1)
    res <- colnames(exprs(tmp))
  }
  names(res) <- NULL
  res
}


#' @rdname FlowHist
#' @examples
#' library(flowPloidyData) 
#' batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
#' batch1
#' @return
#' \code{\link{batchFlowHist}} returns a list of \code{\link{FlowHist}}
#'   objects. 
#' @export
batchFlowHist <- function(files, channel, bins = 256, verbose = TRUE,
                      window = 20, smooth = 20, linearity = "variable",
                      debris = "SC", samples = 2){ 
  res <- list()
  for(i in seq_along(files)){
    if(verbose) message("processing ", files[i])
    tmpRes <- FlowHist(file = files[i], channel = channel, bins = bins,
                       window = window, smooth = smooth, pick = FALSE,
                       linearity = linearity, debris = debris,
                       samples = samples)
    res[[fhFile(tmpRes)]] <- tmpRes
  }              
  return(res)
}

setMethod(
  f = "show",
  signature = "FlowHist",
  def = function(object){
    cat("FlowHist object '")
    cat(fhFile(object)); cat("'\n")
    cat("channel: "); cat(fhChannel(object)); cat("\n")
    cat(fhSamples(object)); cat(" samples"); cat("\n")
    cat("bins: "); cat(fhBins(object)); cat("\n")
    cat("linearity: "); cat(fhLinearity(object)); cat("\n")
    cat("debris: "); cat(fhDebris(object)); cat("\n")
    cat(length(fhComps(object))); cat(" model components: ")
    cat(paste(names(fhComps(object)), collapse = ", ")); cat("\n")
    pnames <- names(formals(fhModel(object)))
    specialP <- names(getSpecialParams(object))
    pnames <- pnames[which(! pnames %in% specialP)]
    cat(length(pnames)); cat(" parameters: ");
    cat(paste(pnames, collapse = ", ")); cat("\n")
    cat(length(specialP)); cat(" special parameters: ");
    cat(paste(specialP, collapse = ", ")); cat("\n")    
    if(length(fhNLS(object)) == 0)
      cat("Model fitting not complete\n")
    else
      cat("Model fit\n")

    if(length(fhCounts(object)) > 0){
      cat(paste("\nAnalysis\n========\n"))
      ## cat(paste("Modelled events: ",
      ##           round(object$counts$total$value, 1)))
      counts <- c(fhCounts(object)$firstPeak$value,
                  fhCounts(object)$secondPeak$value)
      size <- c(coef(fhNLS(object))["Ma"], coef(fhNLS(object))["Mb"])
      if(is.na(size[2])) size <- size[1]
    }
  
    if(length(object@CV) > 0){
    cvs <- c(fhCV(object)$CVa, fhCV(object)$CVb)
    if(!is.null(fhCV(object)$CVb)){
      cat(paste("\nRatio Peak A / Peak B: ", round(fhCV(object)$AB[1], 3),
                ", SE: ", round(fhCV(object)$AB[2], 5), sep = ""))
    }
  }

  if(length(fhCounts(object)) > 0 & length(fhCV(object)) > 0){
    if(length(counts) == 2)
      rnames <- c("Peak A", "Peak B")
    else if (length(counts) == 1)
      rnames <- "Peak A"
    print(kable(data.frame(counts = counts, size = size, cvs = cvs,
                           row.names = rnames), format = "markdown",
                digits = 3))
  }
  
  if(!is.null(fhRCS(object))){
    cat(paste("\nRCS:", round(fhRCS(object), 3), "\n"))
  }

  }
)

########################
## Plotting functions ##
########################
#' Plot the raw data for a FlowHist object
#'
#' Creates a simple plot of the raw histogram data. Used as a utility for
#' other plotting functions, and perhaps useful for users who wish to
#' create their own plotting routines.
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @param ... additional parameters passed to \code{\link{plot}}
#' @return Not applicable, used for plotting
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' plotFH(fh1)
#' @export
plotFH <- function(fh, main = fhFile(fh), ...){
  ## plots the raw data for a FlowHist object
  plot(fhHistData(fh)$intensity, type = 'n', main = main,
       ylab = "Intensity", xlab = fhChannel(fh), ...)
  polygon(x = c(fhHistData(fh)$xx, max(fhHistData(fh)$xx) + 1),
          y = c(fhHistData(fh)$intensity, 0),
          col = "lightgray", border = NA)
}

#' Plot histograms for FlowHist objects
#'
#' Plot histograms for FlowHist objects
#'
#' @param x a \code{\link{FlowHist}} object
#' @param init boolean; if TRUE, plot the regression model using the
#'   initial parameter estimates over the raw data. 
#' @param nls boolean; if TRUE, plot the fitted regression model over the
#'   raw data (i.e., using the final parameter values)
#' @param comps boolean; if TRUE, plot the individual model components
#'   over the raw data.
#' @param ... additional arguments passed on to plot()
#' @return Not applicable
#' @author Tyler Smith
#' @export
plot.FlowHist <- function(x, init = FALSE, nls = TRUE, comps = TRUE,
                          main = fhFile(x), ...){
  plotFH(x, main = main, ...)

  if(init){
    yy <- with(fhHistData(x),
               do.call(fhModel(x),
                       args = c(getSpecialParams(x), fhInit(x)))) 
    lines(x = fhHistData(x)$xx,
          y = yy, 
          col = "grey", lwd = 1, lty = 5)
    points(x = fhInit(x)$Ma,
           y  = fhHistData(x)$intensity[round(fhInit(x)$Ma, 0)],
           cex = 1.5, pch = 16, col = "blue")
    text("A", cex = 1,
         x = fhInit(x)$Ma, col = "blue", pos = 2,
         y = fhHistData(x)$intensity[round(fhInit(x)$Ma, 0)])
    points(x = 2 * fhInit(x)$Ma,
           y = fhHistData(x)$intensity[round(2 * fhInit(x)$Ma, 0)],
           col = "blue", cex = 1.5)
    if(! is.null(fhInit(x)$Mb)){
      points(x = fhInit(x)$Mb,
             y = fhHistData(x)$intensity[round(fhInit(x)$Mb, 0)],
             cex = 1.5, pch = 16, col = "orange")
      text("B", cex = 1,
           x = fhInit(x)$Mb, col = "orange", pos = 2,
           y = fhHistData(x)$intensity[round(fhInit(x)$Mb, 0)])
      points(x = 2 * fhInit(x)$Mb,
             y = fhHistData(x)$intensity[round(2 * fhInit(x)$Mb, 0)],
             col = "orange", cex = 1.5)
    }
    if(! is.null(fhInit(x)$Mc)){
      points(x = fhInit(x)$Mc,
             y = fhHistData(x)$intensity[round(fhInit(x)$Mc, 0)],
             cex = 1.5, pch = 16, col = "darkgreen")
      text("C", cex = 1,
           x = fhInit(x)$Mc, col = "darkgreen", pos = 2,
           y = fhHistData(x)$intensity[round(fhInit(x)$Mc, 0)])
      points(x = 2 * fhInit(x)$Mc,
             y = fhHistData(x)$intensity[round(2 * fhInit(x)$Mc, 0)],
             col = "green", cex = 1.5)
    }
  }

  if(nls & (length(fhNLS(x)) > 0)){
    dat <- tabulateFlowHist(x)
    lines(x = fhHistData(x)$xx[-(1:(fhStart(fhHistData(x)$intensity) -
                                    1))], 
          y = predict(fhNLS(x)), col = 2)
    yPos <- grconvertY(0.95, from = "npc", to = "user") # starting pos
    lHt <- par("cxy")[2]                # line height
    text(paste("RCS: ", round(dat$rcs, 3)), cex = 1, pos = 2,
         x = grconvertX(0.975, from = "npc", to = "user"),
         y = yPos)
    yPos <- yPos - lHt
    text(paste("A: ", round(dat$sizeA, 1), "/",
               round(dat$countsA, 1), "/",
               round(100 * dat$cvA, 1)),
         cex = 1, pos = 2, col = "blue",
         x = grconvertX(0.975, from = "npc", to = "user"),
         y = yPos)
    yPos <- yPos - lHt

    if(!is.na(dat$sizeB)){
      text(paste("B: ", round(dat$sizeB, 1), "/",
                 round(dat$countsB, 1), "/",
                 round(100 * dat$cvB, 1)),
           cex = 1, pos = 2, col = "orange",
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }

    if(!is.na(dat$sizeC)){
      text(paste("C: ", round(dat$sizeC, 1), "/",
                 round(dat$countsC, 1), "/",
                 round(100 * dat$cvC, 1)),
           cex = 1, pos = 2, col = "darkgreen",
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }

    if(is.na(dat$linearity))
      linval <- "fixed"
    else
      linval <- round(dat$linearity, 3)
    
    text(paste("Linearity: ", linval), cex = 1,
         pos = 2, 
         x = grconvertX(0.975, from = "npc", to = "user"),
         y = yPos) 

  }

  if(comps & (length(fhNLS(x)) > 0)){
    yy <- with(fhHistData(x),
               do.call(fhModel(x),
                       args = c(getSpecialParams(x), fhInit(x))))
  
    for(i in seq_along(fhComps(x))){
      params <-
        as.list(coef(fhNLS(x))[names(formals(mcFunc(fhComps(x)[[i]])))])
      params <- params[! is.na(names(params))]
      yy <-
        with(fhHistData(x),
             do.call(mcFunc(fhComps(x)[[i]]),
                      args = c(getSpecialParamsComp(fhComps(x)[[i]]),
                               params)))
      lines(x = fhHistData(x)$xx, y = yy,
            col = mcColor(fhComps(x)[[i]])) 
    }
  }
}

####################
## Exporting Data ##
####################

#' Extract analysis results from a FlowHist object
#'
#' A convenience function for extracting the results of the NLS
#'   curve-fitting analysis on a FlowHist object.
#'
#' If \code{fh} is a single FlowHist object, a data.frame with a single
#' row is returned. If \code{fh} is a list of \code{\link{FlowHist}} objects, a
#' row for each object will be added to the data.frame.
#'
#' If a file name is provided, the data will be saved to that file.
#' 
#' @title exportFlowHist
#' @param fh a FlowHist object, or a list of FlowHist objects.
#' @param file character, the name of the file to save data to
#' @return a data frame 
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1 <- fhAnalyze(fh1)
#' tabulateFlowHist(fh1)
#' @export
tabulateFlowHist <- function(fh, file = NULL){
  if(class(fh) == "FlowHist")
    res <- exFlowHist(fh)
  else if (class(fh) == "list" && all(sapply(fh, class) == "FlowHist")){
    res <- do.call(rbind, lapply(fh, exFlowHist))
  }
  row.names(res) <- res$file
  res$file <- NULL
  
  if(! is.null(file))
    write.table(x = res, file = file)

  res
}

exFlowHist <- function(fh){
  df <- data.frame(file = fhFile(fh), channel = fhChannel(fh),
                   components = paste(names(fhComps(fh)), collapse = ";"),
                   totalEvents = sum(fhHistData(fh)$intensity))
  
  if(length(fhNLS(fh)) > 0){
    df$countsA = fhCounts(fh)$firstPeak$value
    df$countsB = ifelse(is.null(fhCounts(fh)$secondPeak$value), NA,
                        fhCounts(fh)$secondPeak$value)
    df$countsC = ifelse(is.null(fhCounts(fh)$thirdPeak$value), NA,
                        fhCounts(fh)$thirdPeak$value)
    df$sizeA = coef(fhNLS(fh))["Ma"]
    df$sizeB = coef(fhNLS(fh))["Mb"]
    df$sizeC = coef(fhNLS(fh))["Mc"]
    df$cvA = fhCV(fh)$CVa
    df$cvB = ifelse(is.null(fhCV(fh)$CVb), NA, fhCV(fh)$CVb)
    df$cvC = ifelse(is.null(fhCV(fh)$CVc), NA, fhCV(fh)$CVc)
    df$AB = unlist(ifelse(is.null(fhCV(fh)$AB[1]), NA,
                          fhCV(fh)$AB[1]))
    df$ABse = unlist(ifelse(is.null(fhCV(fh)$AB[2]), NA,
                            fhCV(fh)$AB[2]))

    df$AC = unlist(ifelse(is.null(fhCV(fh)$AC[1]), NA,
                          fhCV(fh)$AC[1]))
    df$ACse = unlist(ifelse(is.null(fhCV(fh)$AC[2]), NA,
                            fhCV(fh)$AC[2]))

    df$BC = unlist(ifelse(is.null(fhCV(fh)$BC[1]), NA,
                                     fhCV(fh)$BC[1]))
    df$BCse = unlist(ifelse(is.null(fhCV(fh)$BC[2]), NA,
                            fhCV(fh)$BC[2]))

    df$rcs = fhRCS(fh)

    if(fhLinearity(fh) == "variable")
      df$linearity = coef(fhNLS(fh))["d"]
    else
      df$linearity = NA
    row.names(df) = NULL
  } else {
    df[, c("countsA", "countsB", "countsC", "sizeA", "sizeB", "sizeC",
           "cvA", "cvB", "cvC", "AB", "ABse", "AC", "ACse", "BC", "BCse",
           "rcs", "linearity")] <- NA 
  }
  df
}

#################################################
## Functions for initializing FlowHist objects ##
#################################################

#' (Re-) set the bins for a FlowHist object
#'
#' This function sets (or resets) the number of bins to use in aggregating
#' FCS data into a histogram, and generates the corresponding data matrix.
#'
#' The \code{histData} matrix also contains the \code{SCvals} column. This
#' is used to calculate the single-cut debris component in the NLS model.
#' 
#' @title setBins
#' @param fh a \code{\link{FlowHist}} object
#' @param bins integer, the number of bins to use in aggregating FCS data
#' @return a \code{\link{FlowHist}} object, with the \code{bins} slot set to
#'   \code{bins}, and the corresonding binned data stored in a matrix in
#'   the \code{histData} slot. Any previous analysis slots are removed:
#'   \code{peaks, comps, model, init, nls, counts, CV, RCS}.
#' @author Tyler Smith
#' @examples
#' ## defaults to 256 bins:
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' plot(fh1)
#' ## reset them to 512 bins:
#' fh1 <- setBins(fh1, 512)
#' plot(fh1)
#' @export
setBins <- function(fh, bins = 256){
  fhBins(fh) = as.integer(bins)

  ## Extract the data channel
  chanDat <- exprs(fhRaw(fh))[, fhChannel(fh)]
  gateResid <- NULL
  gate <- fhGate(fh)
  ## remove the top bin - this contains clipped values representing all
  ## out-of-range data, not true values
  chanTrim <- chanDat[chanDat < max(chanDat)]
  gate <- gate[chanDat < max(chanDat)]
  
  if(sum(fhGate(fh)) != 0){
    gateResid <- chanTrim[!gate]
    chanTrim <- chanTrim[gate]
  }
  
  metaData <- pData(parameters(fhRaw(fh)))
  maxBins <- metaData[which(metaData$name == fhChannel(fh)), "range"]
  
  ## aggregate bins: combine maxBins into bins via hist
  binAg <- floor(maxBins / bins)

  histBins <- hist(chanTrim, breaks = seq(from = 0, to = 1024, by = binAg),
                   plot = FALSE)

  intensity <- histBins$counts
  xx <- 1:length(intensity)
  startBin <- fhStart(intensity)
  SCvals <- getSingleCutVals(intensity, xx, startBin)
  MCvals <- getMultipleCutVals(intensity, startBin)
  DBvals <- getDoubletVals(intensity)
  TRvals <- getTripletVals(intensity, DBvals)
  QDvals <- getQuadrupletVals(intensity, DBvals, TRvals)

  if(!is.null(gateResid)){
    gateResid <- hist(gateResid, breaks = histBins$breaks,
                      plot = FALSE)$counts
  } else {
    gateResid <- numeric(length(intensity))
  }
  fhHistData(fh) <- data.frame(xx = xx, intensity = intensity,
                            SCvals = SCvals, MCvals = MCvals,
                            DBvals = DBvals, TRvals = TRvals,
                            QDvals = QDvals, gateResid = gateResid)
  fh <- resetFlowHist(fh)
  fh
}

fhStart <- function(intensity){
  ## Returns the first channel to include in the modelling process. We
  ## start on the first peak, ignoring any noise in lower channels. This
  ## is the same general principle applied in ModFit. I implement this idea
  ## by picking the highest point in the first 20 non-zero channels in the
  ## histogram.
  startMax <- max(intensity[which(intensity != 0)][1:10])
  startBin <- which(intensity == startMax)[1]
  startBin
}
  

#' @importFrom caTools runmean runmax
NULL

##############################
## Peak Detection/Selection ##
##############################

#' findPeaks
#'
#' Locate potential peaks in histogram data
#'
#' Peaks are defined as local maxima in the vector of values, using a
#' moving window. Note that these are used in the context of finding
#' starting values - accuracy isn't important, we just need something
#' `close-enough' that the nls algorithm will be able to find the correct
#' value.
#'
#' Utility functions for use internally by flowPloidy; not exported and
#' won't be visible to users. Usually invoked from within
#' \code{\link{FlowHist}}. 
#'
#' Note that there is a trade-off between accuracy in detected peaks, and
#' avoiding noise. Increasing the value of \code{smooth} will reduce the
#' amount of 'noise' that is included in the peak list. However, increasing
#' smoothing shifts the location of the actual peaks. Most of the time the
#' default values provide an acceptable compromise, given we only need to
#' get 'close enough' for the NLS optimization to find the true parameter
#' values. If you'd like to explore this, the internal (unexported)
#' function \code{fhPeakPlot} may be useful.
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @param window an integer, the width of the moving window to use in
#'   identifying local maxima via \code{\link{caTools::runmax}}
#' @param smooth an integer, the width of the moving window to use in
#'   removing noise via \code{\link{caTools::runmean}}
#' 
#' @return Returns a matrix with two columns:
#' \describe{
#' \item{mean}{the index position of each potential peak}
#' \item{height}{the height (intensity) of the peak at that index position}
#' }
#' 
#' @author Tyler Smith
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' test.dat <-(cumsum(runif(1000, min = -1)))
#' plot(test.dat, type = 'l')
#' test.peaks <- flowPloidy::findPeaks(test.dat, window = 20)
#' points(test.peaks, col = 'red', cex = 2)
#' }
#'
#' @name findPeaks
findPeaks <- function(fh, window = 20, smooth = 20){
  ## extract all peaks from data
  ## smoothing removes most of the noisy peaks
  dat <- fhHistData(fh)[, "intensity"]

  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat

  ## This odd section loops over each "maximum value", and checks to see if
  ## adjacent values in the raw data are actually greater. They may be,
  ## because of the smoothing used above to filter out noise distorts the
  ## position of peak maxima. So we counter this, partially, by sliding the
  ## identified maxima in the smoothed data to the local maxima in the raw
  ## data. I'm only checking the values against their immediate neighbours.
  ## It might work a little better by checking a broader neighbourhood, but
  ## then we may find ourselves reintroducing noise that we only just
  ## finished screening out.

  ## Finding peaks is tricky. I've tried more sophisticated smoothing
  ## algorithms, including FFT, but in removing noise the peaks always end
  ## up shifting, so there's going to be a trade-off involved no matter
  ## what approach I take.
  for(i in which(isMax)){
    while(i == 1 ||                     # skip second test on first bin
          (i < length(dat) && dat[i] < max(dat[c(i - 1, i + 1)], TRUE))){
            if(i == 1)
              break
            else
              if(dat[i + 1] > dat[i -1]){
                isMax[i] <- FALSE
                isMax[i + 1] <- TRUE
                i <- i + 1
              } else {
                isMax[i] <- FALSE
                isMax[i - 1] <- TRUE
                i <- i - 1
              }
          }}

  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)
  fhPeaks(fh) <- res
  fh
}

#' @rdname findPeaks
#'
#' @details \code{\link{cleanPeaks}} filters the output of
#'   \code{\link{findPeaks}} to:   
#' \itemize{
#'
#' \item remove duplicates, ie., peaks with the same intensity that occur
#' within \code{window} positions of each other. Otherwise,
#' \code{\link{findPeaks}} will consider noisy peaks without a single highest
#' point to be multiple distinct peaks.
#'
#' \item drop G2 peaks. In some cases the G2 peak for one sample will have
#' greater intensity than the G1 peak for another sample. We correct for
#' this by removing detected peaks with means close to twice that of other
#' peaks.
#'
#' \item ignore noise, by removing peaks with \code{intensity} < 40. A
#' somewhat arbitrary value. It's tricky to deal with this issue when the
#' debris field is large.
#' }
#' 
cleanPeaks <- function(fh, window = 20){
  ## Remove ties and multiple peaks for histogram analysis

  ## Screen out any ties - if two peaks have the same height, and are
  ## within the same 'window', we need to drop one.
  
  ## If a peak has a 'match' at half the size, use the smaller peak (ie.,
  ## take the G1 peak in cases where the G2 peak is higher) 

  ## After the first peak is selected, only consider peaks that are not a
  ## multiple of the size of this peak when selecting the next one.

  peaks <- fhPeaks(fh)
  peaks <- peaks[order(peaks[,2], decreasing = TRUE), ]

  ## eliminate the debris field?
  peaks <- peaks[which(peaks[, "mean"] > 40), , drop = FALSE]

  drop <- numeric()
  if(nrow(peaks) > 1)
    for(i in 2: nrow(peaks)){
      if((peaks[i-1, "height"] == peaks[i, "height"]) &
         (abs(peaks[i-1, "mean"] - peaks[i, "mean"]) <= window)){ 
        ## It's a tie!
        drop <- c(drop, i)
      }
    }

  if(length(drop) > 0){                  # there was at least one tie 
    peaks <- peaks[-drop, ]
  }
  
  out <- matrix(NA, nrow = 0, ncol = 2)

  while(nrow(peaks) > 0){
    ## which peaks are half or double the size of the first peak:
    paircheck <-
      which(((peaks[, "mean"] < 0.53 * peaks[1, "mean"]) &
             (peaks[, "mean"] > 0.47 * peaks[1, "mean"])) |
            ((peaks[, "mean"] < 2.13 * peaks[1, "mean"]) &
             (peaks[, "mean"] > 1.89 * peaks[1, "mean"])))
    ## Add the first peak to that list:
    paircheck <- c(1, paircheck)
    if(length(paircheck) == 1){            # no pairs
      out <- rbind(out, peaks[1, ])
      peaks <- peaks[-1, , drop = FALSE]              # remove peak
    } else if (length(paircheck == 2)) {              # pick the smallest
                                        # of the pair 
      out <- rbind(out,
                   peaks[paircheck[which.min(peaks[paircheck, "mean"])], ])
      peaks <- peaks[-paircheck, , drop = FALSE]      # remove pair
    } else {
      warning("paircheck found more than 2 peaks")
    }

  }

  if(is.vector(peaks))
    out <- rbind(out, peaks)

  rownames(out) <- NULL

  ## out <- out[1:min(2, nrow(out)), , drop = FALSE]
  if(nrow(out) > 1){
    out <- out[order(out[, "mean"]), ]
  }
  fhPeaks(fh) <- out
  fh
}

#' @title Interactively select model starting values
#'
#' @description Prompts the user to select the peaks to use as initial
#'   values for non-linear regression on a plot of the histogram data. 
#'
#' @details The raw histogram data are plotted, and the user is prompted to
#'   select the peak positions to use as starting values in the NLS
#'   procedure. This is useful when the automated peak-finding algorithm
#'   fails to discriminate between overlapping peaks, or is confused by
#'   noise.
#'
#' The normal use, \code{\link{pickPeaks}} is called from
#'   \code{\link{pickInit}}, rather than directly by the user.
#'
#' @param fh A \code{\link{FlowHist}} object
#' 
#' @return \code{\link{pickInit}} returns the \code{\link{FlowHist}} object
#'   with its initial value slot updated.
#'
#' \code{\link{pickPeaks}} returns a matrix with each peak as a row, with
#'   the mean (position) in the first column, and the height (intensity) in
#'   the second column. 
#'
#' @author Tyler Smith
#'
#' @examples
#' library(flowPloidyData) 
#' fh2 <- FlowHist(file = flowPloidyFiles[12], channel = "FL3.INT.LIN")
#' plot(fh2, init = TRUE) ## automatic peak estimates
#' \dontrun{
#' fh2 <- pickInit(fh2)   ## hand-pick peak estimates
#' }
#' plot(fh2, init = TRUE) ## revised starting values
#' @export
pickInit <- function(fh){
  fhPeaks(fh) <- matrix()
  fhComps(fh) <- list()
  fhModel(fh) <- function(){}
  fhInit(fh) <- list()
  fhNLS(fh) <- structure(list(), class = "nls")
  fhCounts(fh) <- list()
  fhCV(fh) <- list()
  fhRCS(fh) <- NA_real_

  fh <- pickPeaks(fh)
  fh <- addComponents(fh)
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)
  fh
}

pickPeaks <- function(fh){
  ## Does the work of actually plotting and selecting peaks for
  ##   \code{\link{pickInit}}
  if(class(fh) != "FlowHist")
    stop("fh must be a FlowHist object")
  message("plotting data...")
  plotFH(fh)
  message("select sample A peak:")
  peakA <- unlist(locator(1))
  points(peakA[1], peakA[2], col = 2, cex = 3)
  message("select sample B peak:")
  peakB <- unlist(locator(1))
  points(peakB[1], peakB[2], col = 3, cex = 3)
  res <- rbind(peakA, peakB)
  colnames(res) <- c("mean", "height")
  rownames(res) <- NULL
  fhPeaks(fh) <- res
  fh
}

##########################
## Change Model Options ##
##########################
#' Update, and optionally re-analyze, a FlowHist object
#'
#' Allows users to switch the debris model from Single-Cut to Multi-Cut (or
#'   vice-versa), or to toggle linearity between fixed and variable.
#' @title updateFlowHist
#' @param fh a \code{\link{FlowHist}} object
#' @param linearity character, either "variable", the default, or "fixed".
#'   If "fixed", linearity is fixed at 2; if "variable", linearity is fit
#'   as a model parameter.
#' @param debris character, either "SC", the default, or "MC", to set the
#'   debris model component to the Single-Cut or Multi-Cut models.
#' @param analyze boolean, if TRUE the updated model will be analyzed
#'   immediately
#' @param samples integer, the number of samples in the data
#' @return a \code{\link{FlowHist}} object with the modified values of linearity
#'   and/or debris, and, if \code{analyze} was TRUE, a new NLS fitting
#' @author Tyler Smith
#' @examples
#' ## defaults to 256 bins:
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' ## default is Single-Cut, change that to Multi-Cut:
#' fh1mc <- updateFlowHist(fh1, debris = "MC")
#' plot(fh1)
#' @export
updateFlowHist <- function(fh, linearity = NULL, debris = NULL,
                           samples = NULL, analyze = TRUE){
  ## keep the existing peaks, as they may have already been tweaked by the
  ## user
  message("updating FlowHist")

  if(!is.null(linearity))
    if(linearity %in% c("fixed", "variable"))
      fhLinearity(fh) <- linearity
    else
      stop("Invalid linearity value")
  if(!is.null(debris))
    if(debris %in% c("SC", "MC"))
      fhDebris(fh) <- debris
    else
      stop("Invalid debris value")
  if(!is.null(samples))
    if(samples > 0 && samples < 4)
      fhSamples(fh) <- as.integer(samples)
    else
      stop("Invalid sample number: must be between 1 and 3")
  
  fh <- resetFlowHist(fh, from = "comps")
  
  fh <- addComponents(fh)
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)
  if(analyze)
    fh <- fhAnalyze(fh)
  fh
}
