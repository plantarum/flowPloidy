## Initially this was a much larger file containing a number of approaches
## to gating. Now that we've settled on one, it doesn't actually require
## much code. We could probably move this back into FlowHist.R

## This function isn't actually used anymore in FlowHist. It might be
## useful for providing users with a way to visualize gates from the
## console? Or at least, could be extended/repurposed for that use?  
plotGate <- function(fh, x, y, ...){
  chans <- viewFlowChannels(fh)
  if(! x %in% chans)
    stop("Selected flow channel: ", x, " does not exist!")
  if(! y %in% chans)
      stop("Selected flow channel: ", y, " does not exist!")

  dat <- exprs(fhRaw(fh))[, c(x, y)]
  dat[,y] <- dat[,y] / dat[,x]
  plot(dat, pch = 16, col = "#05050510")
}

## Called from the shiny interface
plotResid <- function(fh, main = fhFile(fh), sub = "Gate Residuals", ...){
  plot(fhHistData(fh)$gateResid, type = 'n', main = main,
       sub = sub, xlab = fhChannel(fh), ylab = "Intensity",
       ...) 
  polygon(x = c(fhHistData(fh)$xx, max(fhHistData(fh)$xx), 0), 
          y = c(fhHistData(fh)$gateResid, 0, 0),
          col = "lightgray", border = NA)
}

setGate <- function(fh, gate, refresh = TRUE){
  fhGate(fh) <- gate
  ## We save and restore the peaks here. In most cases, you need to
  ## relocate peaks after setBins, since you're changing the underlying
  ## data the peaks are located in. When setting a gate, however, the peak
  ## position is often/usually going to remain in the same place, the gate
  ## only removes marginal noisy areas (debris). setBins could be modified
  ## to leave the peaks unaltered, but that might lead to abuse in
  ## situations where we want to enforce re-finding peaks.
  peaks <- fhPeaks(fh)
  fh <- setBins(fh, fhBins(fh))
  fhPeaks(fh) <- peaks
  if(refresh){
    fh <- addComponents(fh)
    fh <- setLimits(fh)
    fh <- makeModel(fh)
    fh <- getInit(fh)
    fh <- fhAnalyze(fh)
  }
  fh
}

isGated <- function(fh){
  ## returns TRUE if the FlowHist data is gated
  sum(fhGate(fh)) != 0
}

