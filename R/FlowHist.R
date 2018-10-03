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

#' @importFrom utils write.table file_test
NULL

#' @importFrom methods setClass setMethod new callNextMethod is
NULL

setOldClass("nls")

#' An S4 class to represent internal standard details for
#' \code{\link{FlowHist}} objects
#'
#' The \code{sizes} slot is set in \code{\link{FlowHist}} or
#' \code{\link{batchFlowHist}}. The other values are updates through
#' interaction with the \code{\link{browseFlowHist}} GUI.
#'
#' @name FlowStandards
#'
#' @return \code{\link{stdSizes}}, \code{\link{stdSelected}} and
#'   \code{\link{stdPeak}} return the corresponding slot values
#' 
#' @slot sizes numeric, the size (in pg) of the internal size standard. Can
#'   be a vector of multiple values, if the sample is part of a set that
#'   included different standards for different samples.
#' @slot selected numeric, the size (in pg) of the internal size standard
#'   actually used for this sample. Must be one of the values in the
#'   \code{sizes} slot.
#' @slot peak character, "A" or "B", indicating which of the histogram
#'   peaks is the size standard.
#'
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN",
#'                 standards = c(1.96, 5.43))
#' fhStandards(fh1)  ## display standards included in this object
#' stdSizes(fhStandards(fh1))  ## list standard sizes
setClass(
  Class = "FlowStandards",
  representation = representation(
    sizes = "numeric",
    selected = "numeric",
    peak = "character"
  )
)

#' @rdname FlowStandards
#' @param std a \code{\link{FlowStandards}} object
#' @export
stdSizes <- function(std){
  std@sizes
}

#' @rdname FlowStandards
#' @export
stdSelected <- function(std){
  std@selected
}

`stdSelected<-` <- function(std, value){
  if(value %in% stdSizes(std))
    std@selected <- value
  else
    stop("selected standard size not in list of standard sizes!")
  std
}

#' @rdname FlowStandards
#' @export
stdPeak <- function(std){
  std@peak
}

`stdPeak<-` <- function(std, value){
  std@peak <- value
  std
}

setMethod(
  f = "show",
  signature = "FlowStandards",
  def = function(object){
    sizes <- stdSizes(object)
    sizes <- sizes[!is.na(sizes)]
    selected <- stdSelected(object)
    peak <- stdPeak(object)
    if(length(sizes) == 1){
      cat("standard: ", sizes, "pg")
    } else {
      cat(length(sizes), "standards: ")
      cat(sizes, "\n");
      if(is.numeric(selected))
        cat("set to: ", selected, "\n")
      else
        cat("not set\n")
    }
    if(is.character(peak))
      cat("standard peak: ", peak, "\n")
    else
      cat("standard peak not identified\n")
  }
)

FlowStandards <- function(sizes, selected = 0, peak = "X"){
  if((selected != 0) && ! selected %in% sizes)
    stop("Selected standard size must be in the sizes vector!")

  if(! 0 %in% sizes)
    sizes <- c(0, sizes)
  
  if(length(sizes[sizes != 0]) == 1)
    selected <- sizes[sizes != 0]
  
  new("FlowStandards", sizes = sizes, selected = as.numeric(selected),
      peak = peak)
}


#' FlowHist
#'
#' Creates a \code{\link{FlowHist}} object from an FCS file, setting up the
#' histogram data for analysis.
#'
#' For most uses, simpling calling \code{\link{FlowHist}} with a
#' \code{file}, \code{channel}, and \code{standards} argument will do what
#' you need. The other arguments are provided for optional tuning of this
#' process. In practice, it's easier to correct the model fit using
#' \code{\link{browseFlowHist}} than to determine 'perfect' values to pass
#' in as arguments to \code{\link{FlowHist}}.
#'
#' Similarly, \code{\link{batchFlowHist}} is usually used with only the
#' \code{files}, \code{channel}, and \code{standards} arguments.
#'
#' In operation, \code{\link{FlowHist}} starts by reading an FCS file
#' (using the function \code{\link{read.FCS}} internally). This produces a
#' \code{\link{flowFrame}} object, which we extend to a
#' \code{\link{FlowHist}} object as follows:
#'
#' \enumerate{
#' \item Extract the fluorescence data from \code{channel}.
#'
#' \item Remove the top bin, which contains off-scale readings we ignore in
#' the analysis.
#'
#' \item Remove negative fluorescence values, which are artifacts of
#' instrument compensation
#'
#' \item Removes the first 5 bins, which often contain noisy values,
#' probably further artifacts of compensation.
#'
#' \item aggregates the raw data into the desired number of bins, as
#' specified with the \code{bins} argument. The default is 256, but you may
#' also try 128 or 512. Any integer is technically acceptable, but I
#' wouldn't stray from the default without a good reason. (I've never had a
#' good reason!)
#'
#' \item identify model components to include. All \code{\link{FlowHist}}
#' objects will have the single-cut debris model and the G1 peak for sample
#' A, and the broadened rectangle for the S-phase of sample A. Depending on
#' the data, additional components for the G2 peak and sample B (G1, G2,
#' s-phase) may also be added. The \code{debris} argument can be used to
#' select the Multi-Cut debris model instead, or this can be toggled in
#' \code{\link{browseFlowHist}}
#'
#' \item Build the NLS model. All the components are combined into a single
#' model.
#'
#' \item Identify starting values for Gaussian (G1 and G2 peaks) model
#' components. For reasonably clean data, the built-in peak detection is
#' ok. You can evaluate this by plotting the \code{\link{FlowHist}} object
#' with the argument \code{init = TRUE}. The easiest way to fix bad peak
#' detection is via the \code{\link{browseFlowHist}} interface. You can
#' also play with the \code{window} and \code{smooth} arguments (which is
#' tedious!), or pick the peaks visually yourself with \code{pick = TRUE}.
#'
#' \item Finally, we fit the model and calculate the fitted parameters.
#' Model fitting is suppressed if the \code{analyze} argument is set as
#' \code{FALSE}
#' }
#'
#' @name FlowHist
#'
#' @param file character, the name of the single file to load
#' @param files character, a vector of file names to load, or a single
#'   character value giving the path to a directory; if the latter, all
#'   files in the directory will be loaded 
#' @param channel character, the name of the data column to use
#' @param bins integer, the number of bins to use to aggregate events into
#'   a histogram
#' @param linearity character, either "variable", the default, or "fixed".
#'   If "fixed", linearity is fixed at 2; if "variable", linearity is fit
#'   as a model parameter.
#' @param debris character, either "SC", the default, "MC", or "none", to
#'   set the debris model component to the Single-Cut or Multi-Cut models,
#'   or to not include a debris component (such as for gated data).
#' @param analyze logical, if TRUE the model will be analyzed immediately
#' @param pick logical; if TRUE, the user will be prompted to select peaks
#'   to use for starting values. Otherwise (the default), starting values
#'   will be detected automatically.
#' @param samples integer; the number of samples in the data. Default is 2
#'   (unknown and standard), but can be set to 3 if two standards are used,
#'   or up to 6 for endopolyploidy analysis.
#' @param standards numeric; the size of the internal standard in pg. When
#'   loading a data set where different samples have different standards, a
#'   vector of all the standard sizes. If set to 0, calculation of pg for
#'   the unknown sample will not be done.
#' @param verbose logical; if TRUE, \code{\link{batchFlowHist}} will list
#'   files as it processes them.
#' @param debrisLimit an integer value, default is 40. Passed to
#'   \code{\link{cleanPeaks}}. Peaks with fluorescence values less than
#'   \code{debrisLimit} will be ignored by the automatic peak-finding
#'   algorithm. Used to ignore the debris often found at the left side of
#'   the histogram.
#' @param g2 a logical value, default is TRUE. Should G2 peaks be included
#'   in the model?
#' @param ... additional arguments passed from \code{\link{batchFlowHist}}
#'   to \code{FlowHist}, or to assorted helper functions. See
#'   \code{\link{findPeaks}} (arguments \code{window} and \code{smooth})
#'
#' @slot raw a flowFrame object containing the raw data from the FCS file
#' @slot channel character, the name of the data column to use
#' @slot bins integer, the number of bins to use to aggregate events into a
#'   histogram
#' @slot linearity character, either "fixed" or "variable" to indicate if
#'   linearity is fixed at 2 or fit as a model parameter
#' @slot debris character, either "SC" or "MC" to indicate if the model
#'   should include the single-cut or multi-cut model
#' @slot gate logical, a vector indicating events to exclude from the
#'   analysis. In normal use, the gate will be modified via interactive
#'   functions, not set directly by users.
#' @slot histdata data.frame, the columns are the histogram bin number
#'   (xx), florescence intensity (intensity), and the raw single-cut and
#'   multi-cut debris model values (SCvals and MCvals), and the raw
#'   doublet, triplet and quadruplet aggregate values (DBvals, TRvals, and
#'   QDvals). The debris and aggregate values are used in the NLS fitting
#'   procedures.
#' @slot peaks matrix, containing the coordinates used for peaks when
#'   calculcating initial parameter values.
#' @slot opts list, currently unused. A convenient place to store flags
#'   when trying out new options.
#' @slot comps a list of \code{ModelComponent} objects included for these
#'   data.
#' @slot model the function (built from \code{comps}) to fit to these
#'   data.
#' @slot limits list, a list of lower and upper bounds for model
#'   parameters
#' @slot init a list of initial parameter estimates to use in fitting the
#'   model.
#' @slot nls the nls object produced by the model fitting
#' @slot counts a list of cells counted in each peak of the fitted model
#' @slot CV a list of the coefficients of variation for each peak in the
#'   fitted model.
#' @slot RCS numeric, the residual chi-square for the fitted model.
#' @slot samples numeric, the number of samples included in the data. The
#'   default is 2 (i.e., unknown and standard), but if two standards are
#'   used it should be set to 3. It can be up to 6 for endopolyploidy
#'   analysis, and can be interactively increased (or decreased) via
#'   \code{\link{browseFlowHist}}
#' @slot standards a \code{\link{FlowStandards}} object.
#' @slot g2 logical, if TRUE the model will include G2 peaks for each
#'   sample (as long as the G1 peak is less than half-way across the
#'   histogram). Set to FALSE to drop the G2 peaks for endopolyploidy
#'   analyses.
#' @slot annotation character, user-added annotation for the sample.
#' @slot fail logical, set by the user via the \code{\link{browseFlowHist}}
#'   interface to indicate the sample failed and no model fitting should be
#'   done.
#'
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
    RCS = "numeric", ## residual chi-square
    standards = "FlowStandards", ## a FlowStandards object
    g2 = "logical", ## should G2 peaks be included in the model?
    annotation = "character",
    fail = "logical"
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
                        gate = logical(), samples = 2, standards = 0,
                        opts = list(), debrisLimit = 40, g2 = TRUE,
                        fail = FALSE, ...){
    .Object@raw <- read.FCS(file, dataset = 1, alter.names = TRUE)
    .Object@channel <- channel
    .Object@gate <- gate
    .Object@samples <- as.integer(samples)
    .Object@standards <- FlowStandards(sizes = standards)
    .Object <- setBins(.Object, bins)
    if(pick){
      .Object <- pickPeaks(.Object)
    } else {
      .Object <- findPeaks(.Object, window = window,
                           smooth = smooth)
      .Object <- cleanPeaks(.Object, window = window,
                            debrisLimit = debrisLimit) 
    }
    .Object@linearity <- linearity
    .Object@debris <- debris
    .Object@opts <- opts
    .Object@g2 <- g2
    .Object@fail <- FALSE
    if(! any(is.na(fhPeaks(.Object)))){
      ## We have good peaks:
      .Object <- addComponents(.Object)
      .Object <- setLimits(.Object)
      .Object <- makeModel(.Object)
      .Object <- getInit(.Object)
    } else {
      message("WARNING: couldn't find peaks for ", fhFile(.Object))
    }
    callNextMethod(.Object, ...)
  })

###############
## Accessors ##
###############

#' Functions to access slot values in \code{\link{FlowHist}} objects
#'
#' For normal users, these functions aren't necessary. Overly curious
#' users, or those wishing to hack on the code, may find these useful for
#' inspecting the various bits and pieces inside a \code{\link{FlowHist}}
#' object.
#'
#' The versions of these functions that allow modification of the
#' \code{\link{FlowHist}} object are not exported. Functions are provided
#' for users to update \code{\link{FlowHist}} objects in a safe way.
#'
#' @name fhAccessors
#' 
#' @title FlowHist Accessors
#' @param fh a \code{\link{FlowHist}}
#' @return Used to access a slot, returns the value of the slot. Used to
#'   update the value of a slot, returns the updated \code{\link{FlowHist}}
#'   object.
#' @author Tyler Smith
#' @rdname fhAccessors
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fhModel(fh1) ## prints the model to screen
#' @export
fhGate <- function(fh){
  fh@gate
}

`fhGate<-` <- function(fh, value){
  fh@gate <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhLimits <- function(fh){
  fh@limits
}

`fhLimits<-` <- function(fh, value){
  fh@limits <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhSamples <- function(fh){
  fh@samples
}

`fhSamples<-` <- function(fh, value){
  fh@samples <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhPeaks <- function(fh){
  fh@peaks
}

`fhPeaks<-` <- function(fh, value){
  fh@peaks <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhInit <- function(fh){
  fh@init
}

`fhInit<-` <- function(fh, value){
  fh@init <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhComps <- function(fh){
  fh@comps
}

`fhComps<-` <- function(fh, value){
  fh@comps <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhModel <- function(fh){
  fh@model
}

`fhModel<-` <- function(fh, value){
  fh@model <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhSpecialParams <- function(fh){
  names(getSpecialParams(fh))
}

#' @rdname fhAccessors
#' @export
fhArgs <- function(fh){
  res <- names(formals(fhModel(fh)))
  res <- res[!res %in% fhSpecialParams(fh)]
  res
}

#' @rdname fhAccessors
#' @export
fhNLS <- function(fh){
  fh@nls
}

`fhNLS<-` <- function(fh, value){
  fh@nls <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhCounts <- function(fh){
  fh@counts
}

`fhCounts<-` <- function(fh, value){
  fh@counts <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhCV <- function(fh){
  fh@CV
}

`fhCV<-` <- function(fh, value){
  fh@CV <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhRCS <- function(fh){
  fh@RCS
}

`fhRCS<-` <- function(fh, value){
  fh@RCS <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhFile <- function(fh){
  fh@raw@description$GUID
}

`fhFile<-` <- function(fh, value){
  warning("Changing the FlowHist Data File is a terrible idea")
  warning("-- it would obscure the link between your analysis")
  warning("and the raw data file, and you don't want to do that, do you?")

  fh
}

#' @rdname fhAccessors
#' @export
fhChannel <- function(fh){
  fh@channel
}

`fhChannel<-` <- function(fh, value){
  fh@channel <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhBins <- function(fh){
  fh@bins
}

`fhBins<-` <- function(fh, value){
  fh@bins <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhLinearity <- function(fh){
  fh@linearity
}

`fhLinearity<-` <- function(fh, value){
  fh@linearity <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhDebris <- function(fh){
  fh@debris
}

`fhDebris<-` <- function(fh, value){
  fh@debris <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhHistData <- function(fh){
  fh@histData
}

`fhHistData<-` <- function(fh, value){
  fh@histData <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhRaw <- function(fh){
  fh@raw
}

`fhRaw<-` <- function(fh, value){
  warning("Changing the raw data is a terrible idea")
  warning("It is not supported")
  fh
}

#' @rdname fhAccessors
#' @export
fhStandards <- function(fh){
  fh@standards
}

`fhStandards<-` <- function(fh, value){
  fh@standards <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhStdPeak <- function(fh){
  stdPeak(fhStandards(fh))
}

`fhStdPeak<-` <- function(fh, value){
  stdPeak(fhStandards(fh)) <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhStdSelected <- function(fh){
  stdSelected(fhStandards(fh))
}

`fhStdSelected<-` <- function(fh, value){
  stdSelected(fhStandards(fh)) <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhStdSizes <- function(fh){
  stdSizes(fhStandards(fh))
}

#' @rdname fhAccessors
#' @export
fhOpts <- function(fh){
  fh@opts
}

`fhOpts<-` <- function(fh, value){
  fh@opts <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhG2 <- function(fh){
  fh@g2
}

`fhG2<-` <- function(fh, value){
  fh@g2 <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhAnnotation <- function(fh){
  fh@annotation
}

`fhAnnotation<-` <- function(fh, value){
  fh@annotation <- value
  fh
}

#' @rdname fhAccessors
#' @export
fhFail <- function(fh){
  fh@fail
}

`fhFail<-` <- function(fh, value){
  fh@fail <- value
  fh
}

#' Reset the values in a \code{\link{FlowHist}} object
#'
#' NB: This function isn't required for normal use, and isn't exported for
#' general use. It's provided as a convenience for anyone interested in
#' tweaking model construction and associated parameters. Regular users
#' don't need to do this!
#' 
#' This function provides a safe way to reset the values in a
#' \code{\link{FlowHist}} object. This is important because changing
#' something early in the process will require updating all the dependent
#' values in the appropriate order.
#'
#' The dependency relationships are:
#'
#' \code{gate} <- \code{peaks} <- \code{comps} <- \code{limits}
#'
#' Consequently, changing the \code{gate} requires updating \code{peaks},
#' \code{comps} and \code{limits}. Changing \code{components} only requires
#' updating the \code{limits}. Updating \code{limits} implicitly updates
#' the model and subsequent analysis (i.e., NLS, CV, counts and RCS).
#'
#' In practice, this means that if you change the components, you should
#' call \code{resetFlowHist} to update the dependencies. i.e.,
#' \code{resetFlowHist(fh, from = "limits")}.
#'
#' @param fh a \code{\link{FlowHist}} object.
#' @param from character, the point in the \code{\link{FlowHist}} process
#'   to reset from (see details).
#' @return the updated \code{\link{FlowHist}} object.
#' @author Tyler Smith
#' @keywords internal
resetFlowHist <- function(fh, from = "peaks"){
  ## Clear analysis slots
  ## Default is to clear everything from peaks onwards
  removeFrom <- c("gate", "peaks", "comps", "limits", "analysis")

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
  }
  if(rmF("analysis")){
    fhNLS(fh) <- structure(list(), class = "nls")
    fhCounts(fh) <- list()
    fhCV(fh) <- list()
    fhRCS(fh) <- NA_real_
  }
  fh
}

failFlowHist <- function(fh){
  fh <- resetFlowHist(fh, from = "analysis")
  fhFail(fh) <- TRUE
  fh
}

passFlowHist <- function(fh){
  fhFail(fh) <- FALSE
  fh <- updateFlowHist(fh)
  fh
}

#' @rdname FlowHist
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1
#' @export
FlowHist <- function(file, channel, bins = 256, analyze = TRUE,
                     linearity = "variable", debris = "SC", samples = 2,
                     pick = FALSE, standards = 0, g2 = TRUE,
                     debrisLimit = 40, ...){ 
  fh <-  new("FlowHist", file = file, channel = channel,
             bins = as.integer(bins), linearity = linearity,
             debris = debris, samples = samples, pick = pick,
             standards = standards, g2 = g2, debrisLimit = debrisLimit,
             ...) 
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
#' @param file character, the name of an FCS data file; or the name of a
#'   FlowHist object.
#' @return A vector of column names from the FCS file/FlowHist object.
#' @seealso \code{\link{FlowHist}}
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' viewFlowChannels(flowPloidyFiles[1])
#' @export
viewFlowChannels <- function(file){
  if(is(file, "FlowHist")){
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
#' batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
#' batch1
#' @return
#' \code{\link{batchFlowHist}} returns a list of \code{\link{FlowHist}}
#'   objects. 
#' @export
batchFlowHist <- function(files, channel, verbose = TRUE, ...){ 
  if(length(files) == 1 && file_test("-d", files)){
    files <- list.files(files, full.names = TRUE)
  }
  
  res <- list()
  failures <- character()
  for(i in seq_along(files)){
    if(verbose) message(i, ": processing ", files[i])
    tryVal <- try(tmpRes <- FlowHist(file = files[i], channel = channel,
                                     ...), silent = TRUE)
    if(inherits(tryVal, "try-error")){
      message("    ** PROBLEM: I couldn't import ", files[i], " **")
      failures <- c(failures, files[i])
    } else {
      res[[fhFile(tmpRes)]] <- tmpRes
      if(verbose) message(" ")
    }
  }
  if(length(failures) > 0){
    message("The following files failed to load:")
    cat(paste(failures, collapse = "\n"), "\n")
  }
  return(res)
}

setMethod(
  f = "show",
  signature = "FlowHist",
  def = function(object){
    cat("FlowHist object '")
    cat(fhFile(object)); cat("'\n")
    if(length(fhAnnotation(object)) > 0 && fhAnnotation(object) != "")
      cat("#", fhAnnotation(object), "\n\n")

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
      cat(paste("\nCounts\n------\n"))
      print(vapply(fhCounts(object), function(x) round(x[[1]], 2),
                   numeric(1))) 
      
      cat(paste("\nPeak Sizes\n----------\n"))
      sizeInd <- grep("_mean", names(coef(fhNLS(object))))
      print(round(coef(fhNLS(object))[sizeInd], 2))
      
      cat(paste("\nCVs\n---\n"))
      print(round(100 * vapply(fhCV(object), function(x) x[[1]],
                               numeric(1)), 1))  

    }

    if(!is.null(fhRCS(object))){
      cat(paste("\nRCS:", round(fhRCS(object), 3), "\n"))
    }

  }
)

####################
## Exporting Data ##
####################

#' Extract analysis results from a FlowHist object
#'
#' A convenience function for extracting the results of the NLS
#' curve-fitting analysis on a FlowHist object.
#'
#' If \code{fh} is a single FlowHist object, a data.frame with a single row
#' is returned. If \code{fh} is a list of \code{\link{FlowHist}} objects, a
#' row for each object will be added to the data.frame.
#'
#' If a file name is provided, the data will be saved to that file.
#'
#' The columns of the returned data.frame are:
#'
#' \describe{
#'    \item{countsA, countsB, countsC: }{the cell counts for the G1 peak of
#'    each sample}
#'    \item{sizeA, sizeB, sizeC: }{the peak position for the G1 peak of each
#'    sample}
#'    \item{countsA2, countsB2, countsC2: }{the cell counts for the G2 peak
#'    of each sample}
#'    \item{sizeA2, sizeB2, sizeC2: }{the peak position for the G2 peak of
#'    each sample} 
#'    \item{countsSA, countsSB, countsSC: }{the cell counts for the S-phase
#'    for each sample}
#'    \item{cvA, cvB, cvC: }{the coefficient of variation for each sample}
#'    \item{AB, AC, BC: }{the ratios of sizeA/sizeB, sizeA/sizeC, and
#'    sizeB/sizeC} 
#'    \item{ABse, ACse, BCse: }{the standard error for each ratio}
#'    \item{rcs: }{the residual Chi-Square for the model fit}
#'    \item{linearity: }{the linearity value, if not fixed at 2}
#'    \item{pg: }{genome size estimate, if the sample peak was identified}
#' }
#'
#' @title exportFlowHist
#' @param fh a \code{\link{FlowHist}} object, or a list of
#'   \code{\link{FlowHist}} objects.
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
  if(is(fh, "FlowHist")){
    fhName <- fhFile(fh)
    fh <- list(fh)
    names(fh) <- fhName
  }

  res <- exFlowHist(fh)
  
  if(! is.null(file))
    write.table(x = res, file = file)

  res
}

exFlowHist <- function(fhList, file = NULL){
  samples <- names(fhList)
  if(is.null(samples)){
    ## If the user passes a list of FlowHist objects they created manually,
    ## they need to have names that match fhFile!

    ## fhFile and the names of the list elements are redundant at this
    ## point, it may make sense to drop the names and rely soley on fhFile
    ## for consistency?
    if (is(fhList, "list") && all(vapply(fhList, class, character(1)) ==
                                       "FlowHist")){
      ## list of flowHist objects, but without proper names:
      samples <- names(fhList) <- vapply(fhList, fhFile, character(1))
    } else {
      stop("tabulateFlowHist needs a FlowHist object, or a list of FlowHist objects!")
    }
  }
  outFields <- character()
  coefs <- character()
  counts <- character()
  CVs <- character()
  for(i in fhList){
    if(fhStdPeak(i) != "X")
      outFields <- union(outFields, c("StdPeak", "ratio"))
    if(stdSelected(fhStandards(i)) != 0)
      outFields <- union(outFields, c("StdSize"))
    if(fhStdPeak(i) != "X" && stdSelected(fhStandards(i)) != 0)
      outFields <- union(outFields, c("pg"))
    if(length(fhNLS(i)) > 0){
      coefs <- union(names(coef(fhNLS(i))), coefs)
      counts <- union(names(fhCounts(i)), counts)
      CVs <- union(names(fhCV(i)), CVs)
      outFields <- union(outFields, c("RCS"))
    }
    if(! fhG2(i) && length(fhNLS(i)) > 0){
      outFields <- union(outFields, "EI")
    }
  }
  if(length(fhAnnotation(i)) > 0 && fhAnnotation(i) != "")
    outFields <- union(outFields, c("annotation"))
  
  coefs <- coefs[grep(".*P$", coefs, invert = TRUE)]

  out <- data.frame(row.names = samples)
  for(i in c(outFields, coefs, counts, CVs))
    out[[i]] <- NA

  for(i in fhList){
    if(fhStdPeak(i) != "X" && length(fhNLS(i)) > 0){
      out[fhFile(i), "StdPeak"] <- fhStdPeak(i)
      if(fhStdPeak(i) == "A"){
        out[fhFile(i), "ratio"] <-
          coef(fhNLS(i))["b_mean"]/coef(fhNLS(i))["a_mean"]
      } else if(fhStdPeak(i) == "B"){
        out[fhFile(i), "ratio"] <-
          coef(fhNLS(i))["a_mean"]/coef(fhNLS(i))["b_mean"]
      } else {
        message("More than three peaks, can't calculate ratio")
      }
    }
    if(stdSelected(fhStandards(i)) != 0)
      out[fhFile(i), "StdSize"] <- stdSelected(fhStandards(i))
    if(fhStdPeak(i) != "X" && stdSelected(fhStandards(i)) != 0)
      out[fhFile(i), "pg"] <- stdSelected(fhStandards(i)) *
        out[fhFile(i), "ratio"]
    if(length(fhNLS(i)) > 0){
      coefValues <- coef(fhNLS(i))
      coefNames <- names(coefValues)
      P <- grep(".*P$", coefNames, invert = TRUE)
      coefValues <- coefValues[P]
      coefNames <- coefNames[P]
      out[fhFile(i), coefNames] <- coefValues
      countValues <- vapply(fhCounts(i), function(x) x[[1]], numeric(1))
      countNames <- names(countValues)
      out[fhFile(i), countNames] <- countValues      
      CVs <- fhCV(i)
      CVNames <- names(CVs)
      out[fhFile(i), CVNames] <- CVs
      out[fhFile(i), "RCS"] <- fhRCS(i)
    }
    if(! fhG2(i) && length(fhNLS(i)) > 0){
      means <- sort(names(coef(fhNLS(i)))[grep("_mean",
                                               names(coef(fhNLS(i))))])
      EI <- sum(coef(fhNLS(i))[means] * (seq_len(length(means)) - 1)) /
        sum(coef(fhNLS(i))[means] * (seq_len(length(means))))

      out[fhFile(i), "EI"] <- EI
    }

    if(length(fhAnnotation(i)) > 0 && fhAnnotation(i) != "")
      out[fhFile(i), "annotation"] <- fhAnnotation(i)
  }
  out
}

#################################################
## Functions for initializing FlowHist objects ##
#################################################

#' (Re-)set the bins for a FlowHist object
#'
#' This function sets (or resets) the number of bins to use in aggregating
#' FCS data into a histogram, and generates the corresponding data matrix.
#' Not exported for general use.
#'
#' The \code{histData} matrix also contains the columns corresponding to
#' the raw data used in calculating the single-cut and multiple-cut debris
#' components, as well as the doublet, triplet, and quadruplet aggregate
#' values. (i.e., \code{SCvals}, \code{MCvals}, \code{DBvals},
#' \code{TRvals}, and \code{QDvals}).
#'
#' \code{\link{setBins}} includes a call to \code{\link{resetFlowHist}}, so
#' all the model components that depend on the bins are updated in the
#' process (as you want!).
#' 
#' @title setBins
#' @param fh a \code{\link{FlowHist}} object
#' @param bins integer, the number of bins to use in aggregating FCS data
#' @return a \code{\link{FlowHist}} object, with the \code{bins} slot set
#'   to \code{bins}, and the corresonding binned data stored in a matrix in
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
  if(sum(is.na(chanDat) > 0))
    stop("Error: Flow Data contains missing values!")
  gateResid <- NULL
  gate <- fhGate(fh)
  ## remove the top bin - this contains clipped values representing all
  ## out-of-range data, not true values
  chanTrim <- chanDat[chanDat < max(chanDat)]
  gate <- gate[chanDat < max(chanDat)]
  ## remove values < 0
  ## negative values are artifacts produced by compensation in the
  ## instrument
  posVals <- chanTrim > 0
  chanTrim <- chanTrim[posVals]
  gate <- gate[posVals]

  if(sum(fhGate(fh)) != 0){
    gateResid <- chanTrim[!gate]
    chanTrim <- chanTrim[gate]
  }
  
  metaData <- pData(parameters(fhRaw(fh)))
  maxBins <- metaData[which(metaData$name == fhChannel(fh)), "range"]
  
  ## aggregate bins: combine maxBins into bins via hist
  binAg <- floor(maxBins / bins)

  histBins <- hist(chanTrim, breaks = seq(from = 0, to = maxBins,
                                          by = binAg),
                   plot = FALSE)

  intensity <- histBins$counts
  ## remove the first 5 bins, eliminating noisy artifacts produced by
  ## instrument compensation. This is below the level of actual data, so
  ## shouldn't cause any problems with analysis.
  intensityExclude <- 1:5 ## BiocCheck is overly pedantic on this point
  intensity[intensityExclude] <- 0
  xx <- seq_along(intensity)
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

#' Calculate the where to start analysis for a \code{\link{FlowHist}}
#' histogram 
#'
#' We exclude the first five bins at the outset (as part of the function
#' \code{\link{setBins}}. For some flow cytometers, these values contain
#' very high spikes that are an artifact of compensation, and are not
#' useful data.
#'
#' After that, we call \code{\link{fhStart}} to skip to the highest value
#' in the first 10 non-zero bins, and ignore everything below that. The
#' motivation here is the same - to get out beyond the noisy bins and into
#' the actual data we're trying to fit.
#' 
#' @param intensity numeric, the fluorescence intensity channel bins
#' @return an integer, the index of the first intensity element to include
#'   in the actual model fitting. That is, everything from \code{startBin}
#'   to the end of \code{intensity} gets fit in the model, everything below
#'   \code{startBin} is ignored.
#' @author Tyler Smith
#' @keywords internal
fhStart <- function(intensity){
  ## Returns the first channel to include in the modelling process. We
  ## start on the first peak, ignoring any noise in lower channels. This is
  ## the same general principle applied in ModFit (although I don't know
  ## how they actually do this!). I implement this idea by picking the
  ## highest point in the first 10 non-zero channels in the histogram.
  startRange <- 1:10                    # BiocCheck is dumb
  startMax <- max(intensity[which(intensity != 0)][startRange])
  startBin <- which(intensity == startMax)[1]
  startBin
}

fhStop <- function(intensity){
  ## Returns the last channel to include in the modelling process. We want
  ## to constrain our model to the range of the data. If there is no data
  ## at the highest fluorescence values, including that in the model will
  ## distort our RCS value. We'll consider singleton bins as 'empty', to
  ## clear up stray noise at the top end as well.
  Position(f = function(x) {x > 1}, x = intensity, right = TRUE)
}

fhRange <- function(intensity){
  startM <- fhStart(intensity)
  stopM <- fhStop(intensity)
  startM:stopM
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
#'   identifying local maxima via \code{\link{runmax}}
#' @param smooth an integer, the width of the moving window to use in
#'   removing noise via \code{\link{runmean}}
#' 
#' @return Returns a matrix with two columns:
#' \describe{
#' \item{mean}{the index position of each potential peak}
#' \item{height}{the height (intensity) of the peak at that index position}
#' }
#' 
#' @author Tyler Smith
#' @keywords internal
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
  res <- cbind(mean = (seq_along(dat))[isMax], height = maxVals)
  fhPeaks(fh) <- res
  fh
}

#' @rdname findPeaks
#'
#' @param debrisLimit an integer value. Peaks with fluorescence values less than
#'   \code{debrisLimit} will be ignored by the automatic peak-finding algorithm.
#' 
#' @details \code{\link{cleanPeaks}} filters the output of
#'   \code{\link{findPeaks}} to:   
#' \itemize{
#'
#' \item remove duplicates, ie., peaks with the same intensity that occur
#' within \code{window} positions of each other. Otherwise,
#' \code{\link{findPeaks}} will consider noisy peaks without a single
#' highest point to be multiple distinct peaks.
#'
#' \item drop G2 peaks. In some cases the G2 peak for one sample will have
#' greater intensity than the G1 peak for another sample. We correct for
#' this by removing detected peaks with means close to twice that of other
#' peaks.
#'
#' \item ignore noise, by removing peaks with \code{fluorescence} <
#' \code{debrisLimit}. The default is 40, which works well for
#' moderate-to-large debris fields. You may need to reduce this value if
#' you have clean histograms with peaks below 40. Note that this value does
#' not affect peaks selected manually. }
#' 
cleanPeaks <- function(fh, window = 20, debrisLimit = 40){
  ## Remove ties and multiple peaks for histogram analysis

  ## debrisLimit sets the lower bounds, on the x-axis,k for peak detection
  ## - anything smaller than this is ignored. Needs to be lowered for
  ## histograms with peaks closer to the y axis.

  ## Screen out any ties - if two peaks have the same height, and are
  ## within the same 'window', we need to drop one.
  
  ## If a peak has a 'match' at half the size, use the smaller peak (ie.,
  ## take the G1 peak in cases where the G2 peak is higher) 

  ## After the first peak is selected, only consider peaks that are not a
  ## multiple of the size of this peak when selecting the next one.

  peaks <- fhPeaks(fh)
  peaks <- peaks[order(peaks[,2], decreasing = TRUE), ]

  ## eliminate the debris field?
  if(is.matrix(peaks)){
    peaks <- peaks[which(peaks[, "mean"] > debrisLimit), , drop = FALSE]
  } else {
    ## We only have one or zero potential peaks!
    if (is.vector(peaks) && peaks["mean"] > debrisLimit){
      ## one peak, big enough to keep, convert it to a matrix:
      peaks <- t(as.matrix(peaks))
    } else {
    ## else no peaks, so peaks is now empty...
      peaks <- matrix(NA, nrow = 1, ncol = 2)
      colnames(peaks) <- c("mean", "height")
    }
  }
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

  if (is.vector(peaks)){
    ## only one peak, need to convert the vector back into a matrix:
    peaks <- t(as.matrix(peaks))
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
#' Note that the A peak must be lower (smaller mean, further left) than the
#' B peak. If the user selects the A peak with a higher mean than the B
#' peak, the peaks will be swapped to ensure A is lower.
#'
#' @param fh A \code{\link{FlowHist}} object
#' 
#' @return \code{\link{pickInit}} returns the \code{\link{FlowHist}} object
#'   with its initial value slot updated.
#'
#' @author Tyler Smith
#'
#' @examples
#' library(flowPloidyData) 
#' fh2 <- FlowHist(file = flowPloidyFiles[2], channel = "FL3.INT.LIN")
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
  if(!is(fh, "FlowHist"))
    stop("fh must be a FlowHist object")
  message("plotting data...")
  plotFH(fh)
  message("select sample A peak:")
  peakA <- unlist(locator(1))
  points(peakA[1], peakA[2], col = 2, cex = 3)
  message("select sample B peak:")
  peakB <- unlist(locator(1))
  points(peakB[1], peakB[2], col = 3, cex = 3)

  selectPeaks(fh, peakA[1], peakB[1], NULL)
}

selectPeaks <- function(fh, peakA, peakB = NULL, peakC = NULL,
                        peakD = NULL, peakE = NULL, peakF = NULL){ 
  pA <- fhHistData(fh)[round(peakA, 0), c("xx", "intensity")]
  pB <- pC <- pD <- pE <- pF <- NULL
  if(is.numeric(peakB))                 
    pB <- fhHistData(fh)[round(peakB, 0), c("xx", "intensity")]
  if(is.numeric(peakC))                 
    pC <- fhHistData(fh)[round(peakC, 0), c("xx", "intensity")]
  if(is.numeric(peakD))                 
    pD <- fhHistData(fh)[round(peakD, 0), c("xx", "intensity")]
  if(is.numeric(peakE))                 
    pE <- fhHistData(fh)[round(peakE, 0), c("xx", "intensity")]
  if(is.numeric(peakF))                 
    pF <- fhHistData(fh)[round(peakF, 0), c("xx", "intensity")]
  
  fh <- resetFlowHist(fh)
  newPeaks <- as.matrix(rbind(pA, pB, pC, pD, pE, pF))
  
  colnames(newPeaks) <- c("mean", "height")
  newPeaks <- newPeaks[order(newPeaks[, "mean"]), ]

  ## if we have a single row, the previous selection will return a numeric
  ## vector, which needs to be converted back into a matrix with 1 row:
  if(is.numeric(newPeaks) && ! is.matrix(newPeaks))
    newPeaks <- t(as.matrix(newPeaks))
  
  fhPeaks(fh) <- newPeaks
  
  fh <- addComponents(fh)
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)

  return(fh)
}

##########################
## Change Model Options ##
##########################
#' Update, and optionally re-analyze, a \code{\link{FlowHist}} object
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
#' @param analyze logical, if TRUE the updated model will be analyzed
#'   immediately
#' @param samples integer, the number of samples in the data
#' @return a \code{\link{FlowHist}} object with the modified values of
#'   linearity and/or debris, and, if \code{analyze} was TRUE, a new NLS
#'   fitting
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
    if(debris %in% c("SC", "MC", "none"))
      fhDebris(fh) <- debris
    else
      stop("Invalid debris value")
  if(!is.null(samples))
    if(samples > 0 && samples < 7)
      fhSamples(fh) <- as.integer(samples)
    else
      stop("Invalid sample number: must be between 1 and 6")

  fh <- addComponents(fh)
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)
  if(analyze)
    fh <- fhAnalyze(fh)
  fh
}
