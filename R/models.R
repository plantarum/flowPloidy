## Contains the various model components used in the flowPloidy analysis.

##########################
## ModelComponent Class ##
##########################
setClass(
  Class = "ModelComponent",
  representation = representation(
    name = "character",
    desc = "character", 
    color = "character",
    includeTest = "function",
    ## function with one argument, the FlowHist object.
    ## Return TRUE if the component should be included, FALSE otherwise.
    func = "function",
    ## a single-line function that returns the value of the component.
    ## Can take multiple arguments, usually one of which will be 'xx'
    initParams = "function",
        ## a function that returns a named list of initial parameter
    ## estimates, based on the single argument of the FlowHist object 
    ## list(param1 = param1, ...)
    specialParams = "list",
    ## A named list, the names are parameters to exclude from the default
    ## argument list, as they aren't variables to fit in the NLS procedure.
    ## the body of the list element is the object to insert into the model
    ## formula to account for that variable. e.g., in the singleCut
    ## component, the SCvals parameter is not a variable, and is instead
    ## assigned to the SCvals column in the histData slot. Therefore, it
    ## has a specialParams slot value of `list(SCvals =
    ## substitute(SCvals))`
    specialParamSetter = "function",
    ## function with one argument, the FlowHist object, used to set the
    ## value of specialParams. This allows parameters to be declared
    ## "special" based on values in the fh. i.e., if it fh@linearity is
    ## "fixed", we can declare the parameter d special with a set value of
    ## 2; with linearity == "variable", d is a regular parameter to fit in
    ## the model.
    paramLimits = "list"
    ## A list with the lower and upper limits of each parameter in the
    ## function
  )
)

setMethod(
  f = "show",
  signature = "ModelComponent",
  def = function(object){
    cat("** flowHist model component: ")
    cat(mcName(object)); cat(" ** \n")
    cat(mcDesc(object)); cat(" \n")
    cat("Parameters: ")
    pnames <- names(formals(mcFunc(object)))
    pnames <- pnames[which(pnames != "xx")]
    cat(paste(pnames, collapse = ", "))
    cat("\n")
    if(length(mcSpecialParams(object)) > 0){
      cat("Special Parameters: ")
      cat(paste(names(mcSpecialParams(object)), collapse = ", "))
      cat("\n")
    }
  }
)

###############
## Accessors ##
###############

mcFunc <- function(mc){
  mc@func
}

mcColor <- function(mc){
  mc@color
}

mcName <- function(mc){
  mc@name
}

mcDesc <- function(mc){
  mc@desc
}

mcSpecialParams <- function(mc){
  mc@specialParams
}

`mcSpecialParams<-` <- function(mc, value){
  mc@specialParams <- value
  mc
}

mcSpecialParamSetter <- function(mc){
  mc@specialParamSetter
}

mcIncludeTest <- function(mc){
  mc@includeTest
}

mcInitParams <- function(mc){
  mc@initParams
}

mcParamLimits <- function(mc){
  mc@paramLimits
}

ModelComponent <- function(name, color, desc, includeTest, func,
                           initParams,
                           specialParamSetter = function(fh)
                             list(xx = substitute(xx)),
                           paramLimits = list()){ 
  new("ModelComponent", name = name, color = color, desc = desc,
      includeTest = includeTest, func = func, initParams = initParams,
      specialParamSetter = specialParamSetter, paramLimits = paramLimits)
}

######################
## Model Components ##
######################

## Store all the components in a single, unexported list. This serves as
## our 'menu', which we will search through for each dataset, selecting the
## components that pass the includeTest to add to the components for that
## dataset. 
fhComponents <- list()

## Define new components with the following template:
##
## fhComponents$<name> <-
##   new("ModelComponent", name = "<name>", color = "<colour>",
##       desc = "<one-line description>",
##       includeTest = function(fh){
##
##       },
##       func = function(){
##
##       },
##       initParams = function(fh){
##
##       }
##       )
## 
## specialParamSetter is optional - it will default to a function that
## returns "xx = xx", indicating that all other parameters will be fit. If
## the component doesn't include xx, or includes other fixed parameters,
## then specialParamSetter will need to be provided.

#' Gaussian model components
#'
#' Components for modeling Gaussian features in flow histograms
#'
#' Typically the complete models will contain fA1 and fB2, which model the
#' G1 peaks of the sample and the standard. In most cases, they will also
#' contain fA2 and fB2, which model the G2 peaks. The G2 peaks are linked
#' to the G1 peaks, in that they require some of the parameters from the
#' G1 peaks as well (mean and standard deviation).
#'
#' @param a1,a2,b1,b2 area parameters
#' @param Ma,Mb curve mean parameter
#' @param Sa,Sb curve standard deviation parameter
#' @param xx vector of histogram intensities
#' @return NA
#' @author Tyler Smith
#' @name gauss
fhComponents$fA1 <-
  ModelComponent(
    name = "fA1", color = "blue",
    desc = "Gaussian curve for G1 peak of sample A",
    includeTest = function(fh) {TRUE},
    func = function(a1, Ma, Sa, xx){
      (a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)))
    },
    initParams = function(fh){
      Ma <- as.numeric(fhPeaks(fh)[1, "mean"])
      Sa <- as.numeric(Ma / 20)
      a1 <- as.numeric(fhPeaks(fh)[1, "height"] * Sa / 0.45)
      list(Ma = Ma, Sa = Sa, a1 = a1)
    },
    paramLimits = list(Ma = c(0, Inf), Sa = c(0, Inf), a1 = c(0, Inf))
  )

setLinearity <- function(fh){
  ## Helper function for components that include the linearity parameter d
  if(fh@linearity == "fixed")
    return(list(xx = substitute(xx), d = 2))
  else if(fh@linearity == "variable")
    return(list(xx = substitute(xx)))
  else
    stop("Incorrect setting for linearity")
}

fhComponents$fA2 <-
  ModelComponent(
    name = "fA2", color = "blue",
    desc = "Gaussian curve for G2 peak of sample A",
    includeTest = function(fh){
      (fhPeaks(fh)[1, "mean"] * 2) <= nrow(fhHistData(fh))
    },
    func = function(a2, Ma, Sa, d, xx){
      (a2 / (sqrt(2 * pi) * Sa * 2) *
       exp(-((xx - Ma * d)^2)/(2 * (Sa * 2)^2))) 
    },
    initParams = function(fh){
      Ma <- as.numeric(fhPeaks(fh)[1, "mean"])
      Sa <- as.numeric(Ma / 20)
      a2 <- as.numeric(fhHistData(fh)[Ma * 2, "intensity"] *
                       Sa * 2 / 0.45)
      res <- list(a2 = a2)
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
      res
    },
    paramLimits = list(Ma = c(0, Inf), Sa = c(0, Inf), a2 = c(0, Inf),
                       d = c(1.9, 2.1)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

fhComponents$fB1 <-
  ModelComponent(
    name = "fB1", color = "orange",
    desc = "Gaussian curve for G1 peak of sample B",
    includeTest = function(fh){
      nrow(fhPeaks(fh)) > 1
    },
    func = function(b1, Mb, Sb, xx){
      (b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)))
    },
    initParams = function(fh){
      Mb <- as.numeric(fhPeaks(fh)[2, "mean"])
      Sb <- as.numeric(Mb / 20)
      b1 <- as.numeric(fhPeaks(fh)[2, "height"] * Sb / 0.45)
      list(Mb = Mb, Sb = Sb, b1 = b1)
    },
    paramLimits = list(Mb = c(0, Inf), Sb = c(0, Inf), b1 = c(0, Inf)),
  )

fhComponents$fB2 <-
  ModelComponent(
    name = "fB2", color = "orange",
    desc = "Gaussian curve for G2 peak of sample B",
    includeTest = function(fh){
      if(nrow(fhPeaks(fh)) > 1)
        (fhPeaks(fh)[2, "mean"] * 2) <= nrow(fhHistData(fh))
      else
        FALSE
    },
    func = function(b2, Mb, Sb, d, xx){
      (b2 / (sqrt(2 * pi) * Sb * 2) *
       exp(-((xx - Mb * d)^2)/(2 * (Sb * 2)^2))) 
    },
    initParams = function(fh){
      Mb <- fhPeaks(fh)[2, "mean"]
      Sb <- Mb / 20
      b2 <- as.numeric(fhHistData(fh)[fhPeaks(fh)[2, "mean"] * 2,
                                      "intensity"]
                       * Sb * 2 / 0.45)
      res <- list(b2 = b2)
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
      res
    },
    paramLimits = list(Mb = c(0, Inf), Sb = c(0, Inf), b2 = c(0, Inf),
                       d = c(1.9, 2.1)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

## Single-cut debris model
##
## S(x) = a \sum{j = x + 1}^{n} \sqrt[3]{j} Y_j P_s(j, x)
## P_s (j, x) = \frac{2}{(\pi j \sqrt{(x/j) (1 - x/j)}
##
## a = amplitude parameter
## Y_j = intensity in channel j
## P_s(j, x) = probability of a nuclei from channel j falling into channel x
## when cut.
##
## By this formula, the debris intensity in a channel/bin is a function of
## the intensity in all the subsequent bins. This recursive relationship is
## tricky to code; in order to take full advantage of all of the R tools
## that support nls, the mean function needs to return one fitted value for
## one predictor value. The following implementation of singleCut therefore
## takes the entire vector of the response vector (intensity), necessary to
## calculate the debris curve, and returns only the value for a single
## predictor value.

## Single-cut debris model
##
## Models debris using the single-cut model described by Bagwell et al.
## (1991).
##
## The model is:
## \deqn{S(x) = a \sum{j = x + 1}^{n} \sqrt[3]{j} Y_j P_s(j, x)}
##
## x is the histogram channel that we're estimating the debris value for
## SCa is the amplitude parameter
## Y_j is the histogram intensity for channel j.
##
## where P_s(j, x) is the probability of a nuclei from channel j falling
## into channel x when cut. That is, for j > x, the probability that
## fragmenting a nuclei from channel j with a single cut will produce a
## fragment of size x. This probability is calculated as:
##
## \deqn{P_s (j, x) = \frac{2}{(\pi j \sqrt{(x/j) (1 - x/j)}}}
##
## This model involves a recursive calculation, since the fitted value
## for channel x depends not just on the intensity for channel x, but
## also the intensities at all channels > x. Consequently, this is coded
## with an internal loop, and then vectorized to produce a well-behaved
## function that we can use with the standard nls toolchain.
##
## @name singleCut
##
## @param SCa a numeric value, the single-cut amplitude parameter
## @param intensity a numeric vector, the histogram intensity in each
##   channel 
## @param xx an integer vector, the ordered channels corresponding to the
##   values in `intensity'.
## @param SCvals a numeric vector, stored in the \code{\link{FlowHist}} object
##   slot `SCvals`. Users shouldn't need this.
## @return NA
##
## @references Bagwell, C. B., Mayo, S. W., Whetstone, S. D., Hitchcox,
##   S. 
##   A., Baker, D. R., Herbert, D. J., Weaver, D. L., Jones, M. A. and
##   Lovett, E. J. (1991), DNA histogram debris theory and compensation.
##   Cytometry, 12: 107-118. doi: 10.1002/cyto.990120203
##
## @author Tyler Smith
## @rdname singleCut
getSingleCutValsBase <- function(intensity, xx, first.channel){
  ## compute the single cut debris model values
  
  ## Do not extend the model below/beyond the data
  ## Modfit appears to cut off the debris slightly above the lowest data,
  ## which gives a better fit. Perhaps set first.channel to 2-4? Need to
  ## test this and determine best fit. Possibly use an extra parameter to
  ## tune this for each data set individually.
  ##first.channel <- which(intensity > 0)[2]

  res <- 0
  if(xx >= first.channel & xx < length(intensity)){
    channels = (xx + 1):length(intensity)
    denominator = pi * channels * sqrt(xx/channels * (1 - xx/channels))
    res <- res + sum(channels^(1/3) * intensity[channels] * 2 /
                        denominator)
  }
  res
}

getSingleCutVals <- Vectorize(getSingleCutValsBase, "xx")

fhComponents$SC <-
  ModelComponent(
    name = "SC", color = "green",
    desc = "The single-cut debris model.",
    includeTest = function(fh){
      fhDebris(fh) == "SC"
    },
    func = function(SCa, SCvals){
      SCa * SCvals
    },
    initParams = function(fh){
      list(SCa = 0.1)
    },
    paramLimits = list(SCa = c(0, Inf)),
    specialParamSetter = function(fh){
      list(SCvals = substitute(SCvals))
    }
  )

getMultipleCutVals <- function(intensity, startBin){
  tmpI <- intensity
  res <- sum(tmpI) - cumsum(tmpI)
  res[1:(startBin - 1)] <- 0
  res
}

## getMultipleCutValsBase <- function(a, xx, k){
##   a * exp(-k * xx)
## }

fhComponents$MC <-
  ModelComponent(
    name = "MC", color = "green",
    desc = "The single-cut debris model.",
    includeTest = function(fh){
      fhDebris(fh) == "MC"
    },
    func = function(xx, MCa, k, MCvals){
      MCa * exp(-k * xx) * MCvals[xx]
    },
    initParams = function(fh){
      list(MCa = 0.01, k = 0.001)
    },
    paramLimits = list(MCa = c(0, Inf), k = c(0, Inf)),
    specialParamSetter = function(fh){
      list(xx= substitute(xx), MCvals = substitute(MCvals))
    }
  )

getDoubletVals <- function(intensity){
  doublets <- numeric(length(intensity))
  for(i in seq_along(intensity)[-1]){
    j <- 1:floor(i/2)
    doublets[i] <-
      sum(intensity[j] * intensity[i-j] * (j * (i - j))^(2/3))
  }
  doublets
}

getTripletVals <- function(intensity, doublets){
  triplets <- numeric(length(intensity))
  for(i in seq_along(intensity)[-1]){
    j <- 1:floor(i/2)
    triplets[i] <- 
      sum(intensity[j] * doublets[i-j] * (j * (i - j))^(2/3))
  }
  triplets
}

getQuadrupletVals <- function(intensity, doublets, triplets){
  quadruplets <- numeric(length(intensity))
  for(i in seq_along(intensity)[-1]){
    j <- 1:floor(i/2)
    quadruplets[i] <-
      sum(intensity[j] * triplets[i - j] * (j * (i - j))^(2/3) +
          doublets[j] + doublets[i - j] * (j * (i - j))^(2/3))
  }
  quadruplets
}

fhComponents$AG <-
  ModelComponent(
    name = "AG", color = "purple",
    desc = "Continuous Aggregate",
    includeTest = function(fh){
      TRUE
    },
    func = function(Ap, DBvals, TRvals, QDvals){
      Ap * DBvals + Ap * Ap * TRvals + Ap * Ap * Ap * QDvals
    },
    initParams = function(fh){
      list(Ap = 1e-9)
    },
    paramLimits = list(Ap = c(0, Inf)),
    specialParamSetter = function(fh){
      list(DBvals = substitute(DBvals), TRvals = substitute(TRvals),
           QDvals = substitute(QDvals))
    }
  )

## Broadened rectangles:
## simplified with a fixed sd of 1. Very little change in results with more
## flexible models, so keeping it simple.
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

fhComponents$brA <-
  ModelComponent(
    name = "brA", color = "magenta",
    desc = "Broadened rectangle for S-phase of sample A",
    includeTest = function(fh){
      TRUE
    },
    func = function(BRA, Ma, xx){
      ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
      ## adding sd as a separate parameter
      BRA * ((flowPloidy:::erf(((2 * Ma) - xx)/sqrt(2 * 1)) -
              flowPloidy:::erf((Ma - xx)/sqrt(2 * 1))) / 2)
    },
    initParams = function(fh){
      list(BRA = 10)
    },
    paramLimits = list(BRA = c(0, Inf))
  )

fhComponents$brB <-
  ModelComponent(
    name = "brB", color = "turquoise",
    desc = "Broadened rectangle for S-phase of sample B",
    includeTest = function(fh){
      nrow(fhPeaks(fh)) > 1        
    },
    func = function(BRB, Mb, xx){
      ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
      ## adding sd as a separate parameter
      BRB * ((flowPloidy:::erf(((2 * Mb) - xx)/sqrt(2 * 1)) -
              flowPloidy:::erf((Mb - xx)/sqrt(2 * 1))) / 2)
    },
    initParams = function(fh){
      list(BRB = 10)
    },
    paramLimits = list(BRB = c(0, Inf))
  )

## for testing the influence of sd:
## This isn't used in any other code, retained here for further study if
## needed. 
## brA1 <- function(BRA, Ma, xx, sd){
##   BRA * ((flowPloidy::erf(((2 * Ma) - xx)/sqrt(2 * sd)) -
##           erf((Ma - xx)/sqrt(2 * sd))) / 2)
## }

## The basic broadened trapezoid functions
## Retained here for study, but the complexity doesn't provide much/any
## useful improvement in the model fit.
## broadenedTrapezoid <- function(BTt1, BTt2, BTx1, BTx2, BTs1, BTs2, xx){
##   ((BTt2 - BTt1) / (BTx2 - BTx1) * (xx - BTx2) + BTt2) *
##     ((erf((BTx2 - xx)/sqrt(2 * BTs2)) -
##       erf((BTx1 - xx)/sqrt(2 * BTs1))) / 2)
## }

## ## Translated into model components:
## btA <- function(BTt1A, BTt2A, Ma, xx){
##   ## Simplified to use a fixed sd (5), and bounded to the G1 mean value
##   ## (and by extension the G2 mean value).
##   ((BTt2A - BTt1A) / Ma * (xx - (2 * Ma)) + BTt2A) *
##     ((erf(((2 * Ma) - xx)/sqrt(2 * 5)) -
##       erf((Ma - xx)/sqrt(2 * 5))) / 2)
## }

## btB <- function(BTt1B, BTt2B, Mb, xx){
##   ## Simplified to use a fixed sd (5), and bounded to the G1 mean value
##   ## (and by extension the G2 mean value).
##   ((BTt2B - BTt1B) / Mb * (xx - (2 * Mb)) + BTt2B) *
##     ((erf(((2 * Mb) - xx)/sqrt(2 * 5)) -
##       erf((Mb - xx)/sqrt(2 * 5))) / 2)
## }

##############################
## Model Building Functions ##
##############################
addComponents <- function(fh){
  for(i in fhComponents)
    if(mcIncludeTest(i)(fh)){
      newComp <- i
      mcSpecialParams(newComp) <- mcSpecialParamSetter(newComp)(fh)
      fhComps(fh)[[mcName(i)]] <- newComp
      lims <- mcParamLimits(i)
      newLims <- fhLimits(fh)
      for(j in names(lims))
        newLims[[j]] <- lims[[j]]
      fhLimits(fh) <- newLims
    }
  fh
}

makeModel <- function(fh, env = parent.frame()){
  components <- fhComps(fh)
  names(components) <- NULL
  args <- unlist(lapply(components, FUN = function(x) formals(mcFunc(x))))
  args <- args[unique(names(args))]

  bodList <- lapply(components, FUN = function(x) body(mcFunc(x)))
  bod <- bodList[[1]]
  bodList <- bodList[-1]

  while(length(bodList) > 0){
    bod <- call("+", bod, bodList[[1]])
    bodList <- bodList[-1]
  }

  fhModel(fh) <- eval(call("function", as.pairlist(args), bod), env)
  fhNLS(fh) <- structure(list(), class = "nls")
  fh
}

getInit <- function(fh){
  fhInit(fh) <- list()
  for(i in fhComps(fh)){
    fhInit(fh) <- c(fhInit(fh), mcInitParams(i)(fh))
  }
  fhInit(fh) <- fhInit(fh)[unique(names(fhInit(fh)))]
  fh
}

getSpecialParams <- function(fh){
  res <- list()
  for(i in fhComps(fh))
    res <- c(res, mcSpecialParams(i))
  res[-1 * which(duplicated(res))]
}

getSpecialParamArgs <- function(fh){
  params <- getSpecialParams(fh)
  res <- character()
  for(i in seq_along(params))
    res <- c(res, paste(names(params)[i], " = ", params[[i]]))
  paste(res, collapse = ", ")
}

getSpecialParamsComp <- function(comp){
  mcSpecialParams(comp)
}
