#' @importFrom shiny selectInput

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

plotResid <- function(fh, main = fhFile(fh), sub = "Gate Residuals", ...){
  plot(fhHistData(fh)$gateResid, type = 'n', main = main,
       sub = sub, xlab = fhChannel(fh), ylab = "Intensity",
       ...) 
  polygon(x = c(fhHistData(fh)$xx, max(fhHistData(fh)$xx) + 1),
          y = c(fhHistData(fh)$gateResid, 0),
          col = "lightgray", border = NA)
}

setGate <- function(fh, gate, refresh = TRUE){
  fhGate(fh) <- gate
  fh <- setBins(fh, fhBins(fh))
  if(refresh){
    fh <- findPeaks(fh)
    fh <- cleanPeaks(fh)
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

gateFlowHist <- function(fh){
  raw <- exprs(fhRaw(fh))
  chan1 <- fhChannel(fh)
  chan2 <- viewFlowChannels(fh)[2]      # why 2? sometimes, by chance, the
                                        # second column is the SS value.
  
  initGateData <- data.frame(x = raw[, chan1],
                             y = raw[, chan2] / raw[, chan1]) 

  ui <- fluidPage(
    fluidRow(
      column(2,
             selectInput('xcol', 'X Variable', viewFlowChannels(fh),
                         selected = chan1)),
      column(2,
             selectInput('ycol', 'Y Variable', viewFlowChannels(fh),
                         selected = chan2)),
      column(2,
             sliderInput("yrange", "Zoom", min = 0, ticks = FALSE,
                         step = max(4,
    ceiling(log(max(initGateData$y))))/20, 
                         max = max(4, ceiling(log(max(initGateData$y)))),
                         value = 0,
                         dragRange = FALSE))),
       #      ),
    fluidRow(
      column(4,
             plotOutput("gatePlot", #height = 300,
                        click = "gatePlot_click",
                        brush = brushOpts(
                          id = "gatePlot_brush",
                          resetOnNew = TRUE
                        )
                        )
             ),
      column(4,
             plotOutput("gatedData") #, height = 300)
             ),
    ## fluidRow(
      column(4,
             plotOutput("gateResiduals"))) #, height = 300)))

  )

  server <- function(input, output, session) {
    observe({
        dat <- gateData()
        # Control the value, min, max, and step.
        # Step size is 2 when input value is even; 1 when value is odd.
        updateSliderInput(session, "yrange", 
                          step = max(4, ceiling(log(max(dat[,2]))))/20,
                          max = max(4, ceiling(log(max(dat[,2])))),
                          value = 0, min = 0)
    })

    ##browser()                           
    gateData <- reactive({
      chan1 <- input$xcol
      chan2 <- input$ycol
      
      df <- data.frame(x = raw[, chan1],
                       y = raw[, chan2] / raw[, chan1])
      names(df) <- c(chan1, paste(chan1, chan2, sep = "/"))
      df
    })

    output$gatePlot <- renderPlot({
      dat <- gateData()
      plot(dat, ylim = c(0, exp(log(max(dat[, 2])) - input$yrange)),
           pch = 16, col = "#05050510") 
    })

    output$gateResiduals <- renderPlot({
      plotResid(fhHistPlot(), main = "Gate Residuals")
    })

    output$gatedData <- renderPlot({
      plot(fhHistPlot(), nls = FALSE, init = FALSE, comps = FALSE)
    })
    
    fhHistPlot <- reactive({
      dat <- gateData()
      bp <- brushedPoints(dat, xvar = names(dat)[1], yvar = names(dat)[2],
                          input$gatePlot_brush, allRows = TRUE)$selected_
      fh <<- setGate(fh, bp)
      fh
    })
    
  }
  runApp(shinyApp(ui = ui, server = server))
}
