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
#' @param main character; the plot title. Defaults to the filename of the
#'   \code{\link{FlowHist}} object.
#' @param ... additional parameters passed to \code{\link{plot}}
#' @return Not applicable, used for plotting
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles()[1], channel = "FL3.INT.LIN")
#' plotFH(fh1)
#' @export
plotFH <- function(fh, main = fhFile(fh), ...){
  ## plots the raw data for a FlowHist object
  plot(fhHistData(fh)$intensity, type = 'n', main = main,
       ylab = "Nuclei", xlab = "Fluorescence", ...)
  polygon(x = c(fhHistData(fh)$xx, max(fhHistData(fh)$xx), 0),
          y = c(fhHistData(fh)$intensity, 0, 0),
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
#' @param main character; the plot title. Defaults to the filename of the
#'   \code{\link{FlowHist}} object.
#' @param ... additional arguments passed on to plot()
#' @return Not applicable
#' @author Tyler Smith
#' @export
plot.FlowHist <- function(x, init = FALSE, nls = TRUE, comps = TRUE,
                          main = fhFile(x), ...){
  plotFH(x, main = main, ...)

  if(init && is.list(fhInit(x)) && length(fhInit(x))){
    yy <- with(fhHistData(x),
               do.call(fhModel(x),
                       args = c(getSpecialParams(x), fhInit(x)))) 
    lines(x = fhHistData(x)$xx,
          y = yy, 
          col = "grey", lwd = 1, lty = 5)

    for(comp in letters[seq_len(6)]){
      NAME <- paste(comp, 1, sep = "")
      if(NAME %in% names(fhComps(x))){
        MEAN <- paste(comp, "mean", sep = "_")
        COL <- mcColor(fhComps(x)[[NAME]])
        
        points(x = fhInit(x)[[MEAN]],
               y  = fhHistData(x)$intensity[round(fhInit(x)[[MEAN]], 0)],
               cex = 1.5, pch = 16, col = COL)
        text(toupper(comp), cex = 1,
             x = fhInit(x)[[MEAN]], col = COL, pos = 2,
             y = fhHistData(x)$intensity[round(fhInit(x)[[MEAN]], 0)])

        if(fhG2(x))
          points(x = 2 * fhInit(x)[[MEAN]],
                 y = fhHistData(x)$intensity[round(2 * fhInit(x)[[MEAN]], 0)],
                 col = COL, cex = 1.5)
      }
    }
  } else {
    if(init)
      message("no init values available to plot!!")
  }

  if(fhFail(x)){
    text("** FAIL! **", cex = 1, pos = 2, col = 2,
         x = grconvertX(0.975, from = "npc", to = "user"),
         y = grconvertY(0.95, from = "npc", to = "user"))
  } else if(nls && (length(fhNLS(x)) > 0)){
    dat <- tabulateFlowHist(x)
    lines(x = fhHistData(x)$xx[fhRange(fhHistData(x)$intensity)], 
          y = predict(fhNLS(x)), col = 2, lwd = 2)
    yPos <- grconvertY(0.95, from = "npc", to = "user") # starting pos
    lHt <- par("cxy")[2]                # line height
    if(isGated(x)){
      text("-- GATED --", cex = 1, pos = 2, col = 2,
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }
    text(paste("RCS: ", round(dat$RCS, 2)), cex = 1, pos = 2,
         x = grconvertX(0.975, from = "npc", to = "user"),
         y = yPos)
    yPos <- yPos - lHt

    MEANS <- names(coef(fhNLS(x)))[grep(pattern = "[a-z]_mean",
                                        names(coef(fhNLS(x))))]

    for(i in MEANS){
      L = substring(i, 1, 1)
      COUNT = paste(L, "1_count", sep = "")
      CV = paste(L, "_CV", sep = "")
      NAME = paste(L, 1, sep = "")
      COL <- mcColor(fhComps(x)[[NAME]])
      text(paste(toupper(L), ": ", round(dat[i], 1), "/",
               round(dat[COUNT], 1), "/",
               round(100 * dat[CV], 1)),
           cex = 1, pos = 2, col = COL,
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }

    if(is.null(dat$linearity))
      linval <- "fixed"
    else
      linval <- round(dat$linearity, 2)

    if(fhG2(x)){
      text(paste("Linearity: ", linval), cex = 1,
           pos = 2,
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }

    if(!is.null(dat$EI) && !fhG2(x)){
      text(paste("EI: ", round(dat$EI, 2)), cex = 1,
           pos = 2,
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }
    if(!is.null(dat$pg) && !is.na(dat$pg)){
      text(paste("Sample ", ifelse(dat$StdPeak == "A", "B", "A"),
                 ": ", round(dat$pg, 3), "pg"),
           cex = 1, pos = 2,
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }

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
