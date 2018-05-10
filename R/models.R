## Contains the various model components used in the flowPloidy analysis.

##########################
## ModelComponent Class ##
##########################
#' An S4 class to represent model components
#'
#' \code{\link{ModelComponent}} objects bundle the actual mathematical
#' function for a particular component with various associated data
#' necesarry to incorporate them into a complete NLS model.
#'
#' To be included in the automatic processing of potential model
#' components, a \code{\link{ModelComponent}} needs to be added to the
#' variable \code{fhComponents}.
#' 
#' @name ModelComponent
#'
#' @slot name character, a convenient name with which to refer to the
#'   component
#'
#' @slot desc character, a short description of the component, for human
#'   readers
#'
#' @slot color character, the color to use when plotting the component
#'
#' @slot includeTest function, a function which takes a single argument, a
#'   \code{\link{FlowHist}} object, and returns \code{TRUE} if the
#'   component should be included in the model for that object.
#'
#' @slot function function, a single-line function that returns the value
#'   of the component. The function can take multiple arguments, which
#'   usually will include \code{xx}, the bin number (i.e., x value) of the
#'   histogram. The other arguments are model parameters, and should be
#'   included in the \code{initParams} function.
#'
#' @slot initParams function, a function with a single argument, a
#'   \code{\link{FlowHist}} object, which returns named list of model
#'   parameters and their initial estimates.
#'
#' @slot specialParams list, a named list. The names are variables to
#'   exclude from the default argument list, as they aren't parameters to
#'   fit in the NLS procedure, but are actually fixed values. The body of
#'   the list element is the object to insert into the model formula to
#'   account for that variable. Note that this slot is not set directly,
#'   but should be provided by the value returned by
#'   \code{specialParamSetter} (which by default is \code{list(xx =
#'   substitute(xx))}).
#'
#' @slot specialParamSetter function, a function with one argument, the
#'   \code{\link{FlowHist}} object, used to set the value of
#'   \code{specialParams}. This allows parameters to be declared 'special'
#'   based on values in the \code{\link{FlowHist}} object. The default
#'   value for this slot is a function which returns \code{list(xx =
#'   substitute(xx))} 
#'
#' @slot paramLimits list, a named list with the upper and lower limits of
#'   each parameter in the function.
#'
#' @section Coding Concepts:
#'
#' See the source code file \code{models.R} for the actual code used in
#' defining model components. Here are a few examples to illustrate
#' different concepts.
#'
#' We'll start with the G1 peaks. They are modelled by the components
#' \code{fA1} and \code{fB1} (for the A and B samples). The
#' \code{includeTest} for \code{fA1} is simply \code{function(fh) TRUE},
#' since there will always be at least one peak to fit. \code{fB1} is
#' included if there is more than 1 detected peak, and the setting
#' \code{samples} is more than 1, so the \code{includeTest} is
#' \preformatted{function(fh) nrow(fhPeaks(fh)) > 1 && fhSamples(fh) > 1}
#'
#' The G1 component is defined by the function
#' \preformatted{(a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 *
#'   Sa^2)))} 
#'
#' with the arguments \code{a1, Ma, Sa, xx}. \code{xx} is treated
#' specially, by default, and we don't need to deal with it here. The
#' initial estimates for the other parameters are calculated in
#' \code{initParams}:
#' \preformatted{function(fh){
#'   Ma <- as.numeric(fhPeaks(fh)[1, "mean"])
#'   Sa <- as.numeric(Ma / 20)
#'   a1 <- as.numeric(fhPeaks(fh)[1, "height"] * Sa / 0.45)
#'   list(Ma = Ma, Sa = Sa, a1 = a1)
#' }
#' }
#'
#' \code{Ma} is the mean of the distribution, which should be very close to
#' the peak. \code{Sa} is the standard distribution of the distribution.
#' If we assume the CV is 5\%, that means the \code{Sa} should be 5\% of
#' the distribution mean, which gives us a good first estimate.
#' \code{a1} is a scaling parameter, and I came up with the initial
#' estimate by trial-and-error. Given the other two values are going to be
#' reasonably close, the starting value of \code{a1} doesn't seem to be
#' that crucial.
#' 
#' The limits for these values are provided in \code{paramLimits}.
#' \preformatted{paramLimits = list(Ma = c(0, Inf), Sa = c(0, Inf), a1 =
#'   c(0, Inf))}
#'
#' They're all bound between 0 and Infinity. The upper bound for \code{Ma}
#' and \code{Sa} could be lowered to the number of bins, but I haven't had
#' time or need to explore this yet.
#' 
#' The G2 peaks include the \code{d} argument, which is the ratio of the G2
#' peak to the G1 peak. That is, the linearity parameter:
#' \preformatted{func = function(a2, Ma, Sa, d, xx){
#'   (a2 / (sqrt(2 * pi) * Sa * 2) *
#'     exp(-((xx - Ma * d)^2)/(2 * (Sa * 2)^2))) 
#' }
#' }
#' 
#' \code{d} is the ratio between the G2 and G1 peaks. If \code{linearity =
#' "fixed"}, it is set to 2. Otherwise, it is fit as a model parameter.
#' This requires special handling. First, we check the \code{linearity}
#' value in \code{initParams}, and provide a value for \code{d} if needed:
#' \preformatted{res <- list(a2 = a2)
#' if(fhLinearity(fh) == "variable")
#'     res <- c(res, d = 2)
#' }
#' 
#' Here, \code{a2} is always treated as a parameter, and \code{d} is
#' appended to the initial paramter list only if needed.
#'
#' We also need to use the \code{specialParamSetter} function, in this case
#' calling the helper function \code{setLinearity(fh)}. This function
#' checks the value of \code{linearity}, and returns the appropriate object
#' depending on the result.
#'
#' Note that we use the arguments \code{Ma} and \code{Sa} appear in the
#' \code{function} slot for \code{fA2}, but we don't need to provide their
#' initial values or limits. These values are already supplied in the
#' definition of \code{fA1}, which is always present when \code{fA2} is.
#'
#' NB.: This isn't checked in the code! I know \code{fA1} is always
#' present, but there is no automated checking of this fact. If you create
#' a \code{ModelComponent} that has parameters that are not defined in that
#' component, and are not defined in other components (like \code{Ma} is in
#' this case), you will cause problems. There is also nothing to stop you
#' from defining a parameter multiple times. That is, you could define
#' initial estimates and limits for \code{Ma} in \code{fA1} and \code{fA2}.
#' This may also cause problems. It would be nice to do some
#' sanity-checking to protect against using parameters without defining
#' initial estimates or limits, or providing multiple/conflicting
#' definitions.
#'
#' The Single-Cut Debris component is unusual in two ways. It doesn't
#' include the argument \code{xx}, but it uses the pre-computed values
#' \code{SCvals}. Consequently, we must provide a function for
#' \code{specialParamSetter} to deal with this:
#' \preformatted{specialParamSetter = function(fh){ list(SCvals =
#' substitute(SCvals)) } }
#' 
#' The Multi-Cut Debris component \code{MC} is similar, but it needs to
#' include \code{xx} as a special parameter. The aggregate component
#' \code{AG} also includes several special parameters.
#'
#' For more discussion of the debris components, see
#' \code{\link{DebrisModels}}. 
#' 
#' The code responsible for this is in the file \code{models.R}. Accessor
#' functions are provided (but not exported) for getting and setting
#' \code{\link{ModelComponent}} slots. These functions are named
#'   \code{mcSLOT}, and include \code{mcFunc}, \code{mcColor},
#'   \code{mcName}, \code{mcDesc}, \code{mcSpecialParams},
#'   \code{mcSpecialParamSetter}, \code{mcIncludeTest},
#'   \code{mcInitParams}.  
#'
#' @examples
#' ## The 'master list' of components is stored in fhComponents:
#' flowPloidy:::fhComponents ## outputs a list of component summaries
#'
#' ## adding a new component to the list:
#' \dontrun{
#' fhComponents$pois <-
#'   new("ModelComponent", name = "pois", color = "bisque",
#'       desc = "A poisson component, as a silly example",
#'       includeTest = function(fh){
#'           ## in this case, we check for a flag in the opt slot
#'           ## We could also base the test on some feature of the
#'           ## data, perhaps something in the peaks or histData slots
#'           "pois" %in% fh@opt
#'       },
#'       func = function(xx, plam){
#'           ## The function needs to be complete on a single line, as it
#'           ## will be 'stitched' together with other functions to make
#'           ## the complete model.
#'           exp(-plam)*plam^xx/factorial(xx)
#'       },
#'       initParams = function(fh){
#'           ## If we were to use this function for one of our peaks, we
#'           ## could use the peak position as our initial estimate of
#'           ## the Poisson rate parameter:
#'           plam <- as.numeric(fhPeaks(fh)[1, "mean"])
#'       },
#'       ## bound the search for plam between 0 and infinity. Tighter
#'       ## bounds might be useful, if possible, in speeding up model
#'       ## fitting and avoiding local minima in extremes.
#'       paramLimits = list(plam = c(0, Inf)) 
#'   )
#'
#'   ## specialParamSetter is not needed here - it will default to a
#'   ## function that returns "xx = xx", indicating that all other
#'   ## parameters will be fit. That is what we need for this example. If
#'   ## the component doesn't include xx, or includes other fixed
#'   ## parameters, then specialParamSetter will need to be provided.  
#' 
#'   ## Note that if our intention is to replace an existing component with
#'   ## a new one, we either need to explicitly change the includeTest for
#'   ## the existing component to account for situations when the new one
#'   ## is used instead. As a temporary hack, you could add both and then
#'   ## manually remove one with \code{dropComponents}. 
#'   }
setClass(Class = "ModelComponent",
         representation =
           representation(name = "character",
                          desc = "character",
                          color = "character",
                          includeTest = "function", 
                          func = "function",
                          initParams = "function",
                          specialParams = "list", 
                          specialParamSetter = "function",
                          paramLimits = "list"
                          ))

setMethod(f = "show", signature = "ModelComponent",
          def = function(object){
            cat("** flowHist model component: ")
            cat(mcName(object)); cat(" ** \n")
            cat(mcDesc(object)); cat(" \n")
            cat("Parameters: ")
            pnames <- mcParams(object)
            cat(paste(pnames, collapse = ", "))
            cat("\n")
            if(length(mcSpecialParams(object)) > 0){
              cat("Special Parameters: ")
              cat(paste(names(mcSpecialParams(object)), collapse = ", "))
              cat("\n")
            }
          })

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

mcParams <- function(mc){
  pnames <- names(formals(mcFunc(mc)))
  pnames <- pnames[which(pnames != "xx")]
  pnames
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

###########################################
## Helper Functions for Model Components ##
###########################################

## Set the lower and upper bounds of linearity. We can't leave it
## unconstrained, as in particularly messy histograms it may drift so far
## that the G2 peak is lower than the G1 peak, which is impossible. Not
## sure what the best values to use here actually are. The original range
## of 1.9-2.1 is too constraining for some users.

linL <- 1.5
linH <- 2.5

setLinearity <- function(fh){
  ## Helper function for components that include the linearity parameter d
  if(fh@linearity == "fixed")
    return(list(xx = substitute(xx), d = 2))
  else if(fh@linearity == "variable")
    return(list(xx = substitute(xx)))
  else
    stop("Incorrect setting for linearity")
}

erf <- function(x) {
  ## Helper function used in broadened rectangle features
  2 * pnorm(x * sqrt(2)) - 1
}

##################################################
## Notes on broadened rectangles and trapezoids ##
##################################################

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

makeG1 <- function(l, clr, desc, num){
  ## Generate a G1 peak model component
  v1 <- paste(l, 1, sep = "")
  vM <- paste("M", l, sep = "")
  vS <- paste("S", l, sep = "")

  pL <- list(c(0, Inf), c(0, Inf), c(0, Inf))
  names(pL) <- c(vM, vS, v1)

  makeFun <- function(l){
    tmp <- function(){}
    formals(tmp) <-
      eval(parse(text = sprintf("alist(%s = , %s = , %s = , xx =  )",
                                v1, vM, vS)))
    body(tmp) <-
      parse(text =
        sprintf("(%s / (sqrt(2 * pi) * %s) * exp(-((xx - %s)^2)/(2 * %s^2)))",
                v1, vS, vM, vS))
    return(tmp)
  }
  
  newComp <- ModelComponent(
    name = paste("f", l, 1, sep = ""), color = clr,
    desc = desc,
    includeTest = function(fh){
      nrow(fhPeaks(fh)) >= num && fhSamples(fh) >= num
    },
    func = makeFun(l),
    initParams = function(fh){
      mI <- as.numeric(fhPeaks(fh)[num, "mean"])
      sI <- as.numeric(mI / 20)
      iI <- as.numeric(fhPeaks(fh)[num, "height"] * sI / 0.45)
      res <- list(mI, sI, iI)
      names(res) <- c(vM, vS, v1)
      res
    },
    paramLimits = pL
  )

  return(newComp)
}

makeG2 <- function(l, clr, desc, num){
  ## Generate a G2 peak model component
  v2 <- paste(l, 2, sep = "")
  vM <- paste("M", l, sep = "")
  vS <- paste("S", l, sep = "")

  pL <- list(c(0, Inf), c(0, Inf), c(0, Inf), c(linL, linH))
  names(pL) <- c(vM, vS, v2, "d")

  makeFun <- function(l){
    tmp <- function(){}
    formals(tmp) <-
      eval(parse(text = sprintf("alist(%s = , %s = , %s = , d = , xx = )",
                                v2, vM, vS)))
    body(tmp) <-
      parse(text =
    sprintf("(%s / (sqrt(2 * pi) * %s * 2) * exp(-((xx - %s * d)^2)/(2 * (%s * 2)^2)))",
            v2, vS, vM, vS))
    return(tmp)
  }
  
  newComp <- ModelComponent(
    name = sprintf("f%s2", l), color = clr,
    desc = desc,
    includeTest = function(fh){
      fhG2(fh) && (fhPeaks(fh)[num, "mean"] * 2) <= nrow(fhHistData(fh))
    },
    
    func = makeFun(l),
    initParams = function(fh){
      mI <- as.numeric(fhPeaks(fh)[num, "mean"])
      sI <- as.numeric(mI / 20)
      iI <- as.numeric(fhHistData(fh)[mI * 2, "intensity"] *
                       sI * 2 / 0.45)
      res <- list(iI)
      names(res) <- v2
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
      res
    },
    paramLimits = pL,
    specialParamSetter = function(fh){
      setLinearity(fh)
    }

  )

  return(newComp)
}

#' Gaussian model components
#'
#' Components for modeling Gaussian features in flow histograms
#'
#' Typically the complete models will contain fA1 and fB1, which model the
#' G1 peaks of the sample and the standard. In many cases, they will also
#' contain fA2 and fB2, which model the G2 peaks. The G2 peaks are linked
#' to the G1 peaks, in that they require some of the parameters from the
#' G1 peaks as well (mean and standard deviation).
#'
#' If the linearity parameter is set to "fixed", the G2 peaks will be fit
#' as exactly 2 times the mean of the G1 peaks. If linearity is set to
#' "variable", the ratio of the G2 peaks to the G1 peaks will be fit as a
#' model parameter with an initial value of 2, and constrained to the range
#' 1.5 -- 2.5. (The range is coded as linL and linH. If in doubt, check the
#' values of those, i.e., flowPloidy:::linL, flowPloidy:::linH, to be sure
#' Tyler hasn't changed the range without updating this documentation!!)
#'
#' Additionally, for each set of peaks (sample and standard(s)), a
#' broadened rectangle component is included to model the S-phase. At
#' present, this is component has a single parameter, the height of the
#' rectangle. The standard deviation is fixed at 1. Allowing the SD to vary
#' in the model fitting doesn't make an appreciable difference in my tests
#' so far, so I've left it simple.
#' 
#' @param a1,a2,b1,b2,c1,c2 area parameters
#' @param Ma,Mb,Mc curve mean parameter
#' @param Sa,Sb,Sc curve standard deviation parameter
#' @param xx vector of histogram intensities
#' @param d numeric, the ratio of G2/G1 peak means. When linearity is
#'   fixed, this is set to 2. Otherwise, it is fit as a model parameter
#'   bounded between flowPloidy:::linL and flowPloidy:::linH.
#' @return NA
#' @author Tyler Smith
#' @name gauss
#' @aliases GaussianComponents
fhComponents$fA1 <-
  makeG1("a", "blue", "Gaussian curve for G1 peak of sample A", 1)

fhComponents$fA2 <-
  makeG2("a", "blue", "Gaussian curve for G2 peak of sample A", 1)

## fhComponents$fA2 <-
##   ModelComponent(
##     name = "fA2", color = "blue",
##     desc = "Gaussian curve for G2 peak of sample A",
##     includeTest = function(fh){
##       fhG2(fh) && (fhPeaks(fh)[1, "mean"] * 2) <= nrow(fhHistData(fh))
##     },
##     ## This hard-codes the StdDev of the G2 peak to be exactly twice the
##     ## StdDev of the G1 peak, as suggested by bagwell_1993. Maybe this
##     ## should actually be set to the parameter d?
##     func = function(a2, Ma, Sa, d, xx){
##       (a2 / (sqrt(2 * pi) * Sa * 2) *
##        exp(-((xx - Ma * d)^2)/(2 * (Sa * 2)^2))) 
##     },
##     initParams = function(fh){
##       Ma <- as.numeric(fhPeaks(fh)[1, "mean"])
##       Sa <- as.numeric(Ma / 20)
##       a2 <- as.numeric(fhHistData(fh)[Ma * 2, "intensity"] *
##                        Sa * 2 / 0.45)
##       res <- list(a2 = a2)
##       if(fhLinearity(fh) == "variable")
##         res <- c(res, d = 2)
##       res
##     },
##     paramLimits = list(Ma = c(0, Inf), Sa = c(0, Inf), a2 = c(0, Inf),
##                        d = c(linL, linH)),
##     specialParamSetter = function(fh){
##       setLinearity(fh)
##     }
##   )

fhComponents$brA <-
  ModelComponent(
    name = "brA", color = "magenta",
    desc = "Broadened rectangle for S-phase of sample A",
    includeTest = function(fh){
      TRUE
    },
    func = function(BRA, Ma, d, xx){
      ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
      ## adding sd as a separate parameter
      BRA * ((flowPloidy:::erf(((d * Ma) - xx)/sqrt(2 * 1)) -
              flowPloidy:::erf((Ma - xx)/sqrt(2 * 1))) / 2)
    },
    initParams = function(fh){
      res <- list(BRA = 10)
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
      res
    },
    ## paramLimits = list(BRA = c(0, Inf))
    paramLimits = list(BRA = c(0, Inf), d = c(linL, linH)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

fhComponents$fB1 <-
  makeG1("b", "orange", "Gaussian curve for G1 peak of sample B", 2)

fhComponents$fB2 <-
  ModelComponent(
    name = "fB2", color = "orange",
    desc = "Gaussian curve for G2 peak of sample B",
    includeTest = function(fh){
      if(fhG2(fh) && nrow(fhPeaks(fh)) > 1 && fhSamples(fh) > 1)
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
                       d = c(linL, linH)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

fhComponents$brB <-
  ModelComponent(
    name = "brB", color = "turquoise",
    desc = "Broadened rectangle for S-phase of sample B",
    includeTest = function(fh){
      nrow(fhPeaks(fh)) > 1 && fhSamples(fh) > 2
    },
    func = function(BRB, Mb, d, xx){
      ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
      ## adding sd as a separate parameter
      BRB * ((flowPloidy:::erf(((d * Mb) - xx)/sqrt(2 * 1)) -
              flowPloidy:::erf((Mb - xx)/sqrt(2 * 1))) / 2)
    },
    initParams = function(fh){
      res <- list(BRB = 10)
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
      res
    },
    paramLimits = list(BRB = c(0, Inf), d = c(linL, linH)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

fhComponents$fC1 <-
  makeG1("c", "darkgreen", "Gaussian curve for G1 peak of sample C", 3)

fhComponents$fC2 <-
  ModelComponent(
    name = "fC2", color = "darkgreen",
    desc = "Gaussian curve for G2 peak of sample C",
    includeTest = function(fh){
      fhG2(fh) && nrow(fhPeaks(fh)) > 2 && fhSamples(fh) > 2 &&
        (fhPeaks(fh)[3, "mean"] * 2) <= nrow(fhHistData(fh))
    },
    func = function(c2, Mc, Sc, d, xx){
      (c2 / (sqrt(2 * pi) * Sc * 2) *
       exp(-((xx - Mc * d)^2)/(2 * (Sc * 2)^2))) 
    },
    initParams = function(fh){
      Mc <- as.numeric(fhPeaks(fh)[3, "mean"])
      Sc <- as.numeric(Mc / 20)
      c2 <- as.numeric(fhHistData(fh)[Mc * 2, "intensity"] *
                       Sc * 2 / 0.45)
      res <- list(c2 = c2)
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
      res
    },
    paramLimits = list(Mc = c(0, Inf), Sc = c(0, Inf), c2 = c(0, Inf),
                       d = c(linL, linH)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

fhComponents$brC <-
  ModelComponent(
    name = "brC", color = "magenta",
    desc = "Broadened rectangle for S-phase of sample C",
    includeTest = function(fh){
      FALSE
    },
    func = function(BRC, Mc, d, xx){

      ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
      ## adding sd as a separate parameter
      
      ## WARNING: do I need to replace '2' with 'd' to account
      ## for variable linearity? Bagwell uses "2" in his formula in
      ## bagwell_1993, but I'm not sure if that's because d is usually 2,
      ## or if this is a coincidence. 
      BRC * ((flowPloidy:::erf(((d * Mc) - xx)/sqrt(2 * 1)) -
              flowPloidy:::erf((Mc - xx)/sqrt(2 * 1))) / 2)
    },
    initParams = function(fh){
      res <- list(BRC = 10)
      if(fhLinearity(fh) == "variable")
        res <- c(res, d = 2)
    },
    paramLimits = list(BRC = c(0, Inf), d = c(linL, linH)),
    specialParamSetter = function(fh){
      setLinearity(fh)
    }
  )

#' Histogram Debris Models
#'
#' Implementation of debris models described by Bagwell et al. (1991).
#'
#' @section Single Cut Model:
#'
#' This is the theoretical probability distribution of the size of pieces
#' formed by a single random cut through an ellipsoid. In other words, we
#' assume that the debris is composed of nuclei pieces generated by cutting
#' a subset of the nuclei in a sample into two pieces.
#' 
#' The model is:
#' \deqn{S(x) = a \sum_{j = x + 1}^{n} \sqrt[3]{j} Y_j P_s(j, x)}
#'
#' \enumerate{
#' \item \code{x} the histogram channel that we're estimating the debris
#' value for.
#' \item \code{SCa} the amplitude parameter.
#' \item \code{Y_j} the histogram intensity for channel j.
#' }
#' 
#' where P_s(j, x) is the probability of a nuclei from channel j falling
#' into channel x when cut. That is, for j > x, the probability that
#' fragmenting a nuclei from channel j with a single cut will produce a
#' fragment of size x. This probability is calculated as:
#'
#' \deqn{P_s(j, x) = \frac{2}{(\pi j \sqrt{(x/j) (1 - x/j)}}}
#'
#' This model involves a recursive calculation, since the fitted value
#' for channel x depends not just on the intensity for channel x, but also
#' the intensities at all channels > x. I deal with this by pre-calculating
#' the raw values, which don't actually depend on the only parameter,
#' \code{SCa}. These raw values are stored in the \code{histData} matrix
#' (which is a slot in the \code{\link{FlowHist}} object). This must be
#' accomodated by treating \code{SCvals} as a 'special parameter' in the
#' \code{\link{ModelComponent}} definition. See that help page for details.
#'
#' @section Multiple-Cut Model:
#'
#' The Multiple-Cut model extends the Single-Cut model by assuming that a
#' single nuclei may be cut multiple times, thus creating more than two
#' fragments.
#'
#' The model is:
#' \deqn{S(x) = MCa e^{-kx}\sum_{j = x + 1}^{n} Y_j}
#'
#' \enumerate{
#' \item \code{x} the histogram channel that we're estimating the debris
#' value for.
#' \item \code{k} an exponential fitting parameter
#' \item \code{MCa} the amplitiude parameter
#' \item \code{Y_j} the histogram intensity for channel j.
#' }
#'
#' This model involves another recursive or "histogram-dependent"
#' component. Again, the sum is independent of the fitted parameters, so we
#' can pre-compute that and add it to the \code{histData} slot, in the
#' column \code{MCvals}. This is treated as a 'special parameter' when the
#' Multiple-Cut model is applied, so we only need to fit the parameters k
#' and MCa.
#'
#' @section Debris Models and Gating:
#'
#' The debris models assume that all debris is composed of nuclei (G1 and
#' G2), that have been cut into 2 or more fragments. In actual practice, at
#' least when working with plant cells, the debris likely also includes
#' other cellular debris, including secondary compounds. This non-nuclear
#' debris may take up, and interact with, the stain in unpredictable ways.
#' In extreme cases, such as the Vaccinium example in the ``flowPloidy
#' Getting Started'' vignette, this cellular debris can completely obscure
#' the G1 and G2 peaks, requiring gating.
#'
#' The ideal gate would be one that excludes all of the non-nuclear debris,
#' and none of the nuclear debris (i.e., the nuclei fragments). If we could
#' accomplish this, then gating would improve our model-fitting. Leaving
#' non-nuclear debris in the data will result in it getting fit by some
#' combination of the model components, with a negative impact on their
#' accuracy. On the other hand, excluding nuclear debris will reduce the
#' information used to fit the SC or MC components, which will also reduce
#' model accuracy.
#'
#' Of course, we can't define an ideal gate, anymore than we can optimize
#' our sample preparation such that our histograms are completely free of
#' debris. As a practical approach, we recommend avoiding gating whenever
#' possible, and taking a conservative approach when it is unavoidable.
#'
#' @name DebrisModels
#'
#' @param intensity a numeric vector, the histogram intensity in each
#'   channel
#' @param xx an integer vector, the ordered channels corresponding to the
#'   values in `intensity'.
#' @param first.channel integer, the lowest bin to include in the modelling
#'   process. Determined by the internal function \code{fhStart}.
#' @return \code{getSingleCutVals}, the vectorized function built from
#'   getSingleCutValsBase, returns the fixed \code{SCvals} for the
#'   histogram.
#'
#' \code{getMultipleCutVals}, a vectorized function, returns the
#'   fixed \code{MCvals} for the histogram.
#' 
#' @references Bagwell, C. B., Mayo, S. W., Whetstone, S. D., Hitchcox, S.
#' A., Baker, D. R., Herbert, D. J., Weaver, D. L., Jones, M. A. and
#' Lovett, E. J. (1991), DNA histogram debris theory and compensation.
#' Cytometry, 12: 107-118. doi: 10.1002/cyto.990120203
#'
#' @author Tyler Smith
#'
#' @keywords internal
#' 
#' @examples
#' ## This is an internal function, called from setBins()
#' \dontrun{
#'   ## ...
#'   SCvals <- getSingleCutVals(intensity, xx, startBin)
#'   MCvals <- getMultipleCutVals(intensity, startBin)
#'   ## ...
#'   fhHistData(fh) <- data.frame(xx = xx, intensity = intensity,
#'                            SCvals = SCvals, MCvals = MCvals,
#'                            DBvals = DBvals, TRvals = TRvals,
#'                            QDvals = QDvals, gateResid = gateResid)
#'   ## ...
#' }
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

#' @rdname DebrisModels
#' @name SingleCut
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

#' @rdname DebrisModels
#' @name getMultipleCutVals
getMultipleCutVals <- function(intensity, first.channel){
  tmpI <- intensity
  res <- sum(tmpI) - cumsum(tmpI)
  res[1:(first.channel - 1)] <- 0
  res
}

#' @rdname DebrisModels
#' @name MultipleCut
fhComponents$MC <-
  ModelComponent(
    name = "MC", color = "green",
    desc = "The single-cut debris model.",
    includeTest = function(fh){
      fhDebris(fh) == "MC"
    },
    func = function(xx, MCa, k, MCvals){
      MCa * exp(-k * xx) * MCvals ##[xx]
    },
    initParams = function(fh){
      list(MCa = 0.01, k = 0.001)
    },
    paramLimits = list(MCa = c(1e-10, Inf), k = c(1e-10, Inf)),
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


##############################
## Model Building Functions ##
##############################

#' Functions for assembling non-linear regression models for
#' \code{\link{FlowHist}} objects.
#'
#' \code{\link{addComponents}} examines the model components in
#' \code{fhComponents} and includes the ones that pass their
#' \code{includeTest}.
#'
#' \code{\link{dropComponents}} removes a component from the
#' \code{\link{FlowHist}} model
#'
#' \code{\link{setLimits}} collates the parameter limits for the model
#'   components included in a \code{\link{FlowHist}} object. (could be
#'   called automatically from \code{\link{addComponents}}, as it already
#'   is from \code{\link{dropComponents}}?)
#'
#' \code{\link{makeModel}} creates a model out of all the included
#' components. 
#' 
#' @title Building Flow Histogram Models
#'
#' @name fhModels
#'
#' @aliases flowModels
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @return The updated \code{\link{FlowHist}} object.
#' @author Tyler Smith
addComponents <- function(fh){
  ## make sure old components are flushed!
  fh <- resetFlowHist(fh, from = "comps")
  for(i in fhComponents)
    if(mcIncludeTest(i)(fh)){
      newComp <- i
      mcSpecialParams(newComp) <- mcSpecialParamSetter(newComp)(fh)
      fhComps(fh)[[mcName(i)]] <- newComp
      ## lims <- mcParamLimits(i)
      ## newLims <- fhLimits(fh)
      ## for(j in names(lims))
      ##   newLims[[j]] <- lims[[j]]
      ## fhLimits(fh) <- newLims
    }
  if(fhLinearity(fh) == "variable")
    if(sum(c("fA2", "fB2", "fC2") %in% names(fhComps(fh))) == 0){
      message("No G2 peaks, using fixed linearity")
      fh <- updateFlowHist(fh, linearity = "fixed")
    }
  fh
}

#' @rdname fhModels
#' @param components character, a vector of \code{\link{ModelComponent}}
#'   names.
dropComponents <- function(fh, components){
  fh <- resetFlowHist(fh, "limits")  
  fhComps(fh) <- fhComps(fh)[! names(fhComps(fh)) %in% components]
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)
  fh
}

#' @rdname fhModels
setLimits <- function(fh){
  fhLimits(fh) <- list()
  for(i in fhComps(fh)){
    lims <- mcParamLimits(i)
    for(j in names(lims)){
      fhLimits(fh)[[j]] <- lims[[j]]
    }
  }
  fh
}

#' @rdname fhModels
#' @param env an R environment. Don't change this, it's R magic to keep the
#'   appropriate environment in scope when building our model.
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
