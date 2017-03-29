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
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
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
    if(isGated(x)){
      text("-- GATED --", cex = 1, pos = 2, col = 2,
           x = grconvertX(0.975, from = "npc", to = "user"),
           y = yPos)
      yPos <- yPos - lHt
    }
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
    yPos <- yPos - lHt
    
    if(!is.na(dat$pg)){
      text(paste("Sample ", ifelse(dat$stdpeak == "A", "B", "A"),
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
