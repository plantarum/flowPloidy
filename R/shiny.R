#' @importFrom shiny fluidPage nearPoints reactive radioButtons
#'   actionButton plotOutput reactive eventReactive shinyApp titlePanel
#'   sidebarLayout sidebarPanel htmlOutput fluidRow tags mainPanel
#'   renderPrint renderTable renderPlot renderText column observe runApp
#'   stopApp wellPanel updateRadioButtons HTML numericInput sliderInput
#'   brushOpts eventReactive
NULL

#' @importFrom utils str
NULL

#' Visually assess histogram fits
#'
#' Visually assess histogram fits
#' @title browseFlowHist
#' @param flowList either a \code{\link{FlowHist}} object, or a list of
#'   \code{\link{FlowHist}} objects 
#' @param debug boolean, turns on debugging messages
#' @return Returns the list of \code{\link{FlowHist}} objects, updated by any
#'   changes made in the GUI.
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData)
#' batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
#' \dontrun{
#' batch1 <- browseFlowHist(batch1)
#' }
#' @export
browseFlowHist <- function(flowList, debug = FALSE){
  if(class(flowList) == "FlowHist"){
    flowList <- list(flowList)
    names(flowList) <- fhFile(flowList[[1]])
  }
  
  .fhI <- 1
  .fhList <- flowList

  initialLinearity <- fhLinearity(.fhList[[1]])

  if(debug) message("init Linearity: ", initialLinearity)
  
  initialDebris <- fhDebris(.fhList[[1]])

  initialSamples <- fhSamples(.fhList[[1]])
  
  if(debug) message("init Debris: ", initialDebris)

  raw <- exprs(fhRaw(.fhList[[.fhI]]))
  chan1 <- fhChannel(.fhList[[.fhI]])
  chan2 <- viewFlowChannels(.fhList[[.fhI]])[2]      # why 2? sometimes, by
                                        # chance, the 
                                        # second column is the SS value.
  
  initGateData <- data.frame(x = raw[, chan1],
                             y = raw[, chan2] / raw[, chan1]) 
  
  ui <- fluidPage(
    tags$head(
           tags$style(HTML("
      .sidepanel {
          max-width: 300px;
      }
      #gatePlot, #gatedData, #gateResiduals {
          max-height: 300px;
      }
    "))),
    fluidRow(
      column(width = 3,
             tags$div(class = "sidepanel", 
             wellPanel(
               fluidRow(
                 column(4,
                        actionButton("exit", label = "Exit")),
                 column(8, 
                        fluidRow(
                          column(12, htmlOutput("flowNumber",
                                                align = "center"))),
                        fluidRow(
                          column(5, actionButton("prev", label = "Prev")),
                          column(7, actionButton("nxt",
                                                 label = "Next"))))),
               tags$hr(),
               fluidRow(
                 column(4, 
                        numericInput('sampSelect', 'Samples',
                                     initialSamples, min = 2, max = 3)),
                 column(8,
                        radioButtons(inputId = "peakPicker",
                                     label = "Move peak:", 
                                     choices = list("A" = "A",
                                                    "B" = "B",
                                                    "C" = "C"), 
                                     selected = "A", inline = TRUE))),
               tags$hr(),
               radioButtons(inputId = "linearity",
                            label = "Linearity",
                            choices = list("Fixed" = "fixed",
                                           "Variable" = "variable"), 
                            inline = TRUE, selected = initialLinearity),
               tags$hr(),
               radioButtons(inputId = "debris",
                            label = "Debris Model",
                            choices = list("MC" = "MC", "SC" = "SC"),  
                            inline = TRUE, selected = initialDebris)
             ))),
      column(width = 9,
             plotOutput("fhHistogram", click = "pointPicker"))
    ),
    fluidRow(
      column(width = 3,
             tags$div(class = "sidepanel",
             wellPanel(
               fluidRow(
                 column(6,
                        selectInput('xcol', 'X Variable',
                                    viewFlowChannels(.fhList[[.fhI]]), 
                                    selected = chan1)),
                 column(6,
                        selectInput('ycol', 'Y Variable',
                                    viewFlowChannels(.fhList[[.fhI]]), 
                                    selected = chan2))),
               sliderInput("yrange", "Zoom", min = 0, ticks = FALSE,
                           step =
                             max(4, ceiling(log(max(initGateData$y))))/20, 
                           max = max(4, ceiling(log(max(initGateData$y)))),
                           value = 0, dragRange = FALSE),
               actionButton("setGate", label = "Set Gate")))),
      column(width = 3,
             plotOutput("gatePlot",
                        click = "gatePlot_click",
                        brush = brushOpts(id = "gatePlot_brush",
                                          resetOnNew = FALSE
                                          )
                        )
             ),
      column(width = 3,
             plotOutput("gatedData")
             ),
      column(width = 3,
             plotOutput("gateResiduals")))

  )
    
  server <- function(input, output, session){
    rv <- reactiveValues(fhI = 1, FH = .fhList[[1]])

    ## eventReactive would require an invalid reactive value in the
    ## expression to trigger the calculation; observeEvent will simply do
    ## the calculation:
    fhNext <- observeEvent(input$nxt, {
      if(rv$fhI < length(.fhList))
        rv$fhI <- rv$fhI + 1
    })

    fhPrev <- observeEvent(input$prev, {
      if(rv$fhI > 1)
        rv$fhI <- rv$fhI - 1
    })

    fhCurrent <- eventReactive(rv$fhI, {
      ## When navigating to a new FlowHist object via Prev/Next, update the
      ## radio buttons before updating rv$FH, so the replotting isn't
      ## triggered until the Radio Buttons are set to the values of the
      ## current FH object. The reaction chain is:

      ## fhNext/fhPrev --> rv$fhI --> fhCurrent --> rv$FH
      
      updateRadioButtons(session, "linearity",
                         selected = fhLinearity(.fhList[[rv$fhI]]))

      updateRadioButtons(session, "debris",
                         selected = fhDebris(.fhList[[rv$fhI]]))

      updateNumericInput(session, "sampSelect",
                         value = fhSamples(.fhList[[rv$fhI]]))
      rv$FH <- .fhList[[rv$fhI]]
      rv$fhI
    })      

    fhSetGate <- observeEvent(input$setGate, { 
      if(debug) message(prefix, "Setting gate")
      dat <- gateData()
      bp <- brushedPoints(dat, xvar = names(dat)[1],
                          yvar = names(dat)[2], input$gatePlot_brush,
                          allRows = TRUE)$selected_
      if(sum(bp) == 0) bp <- logical()
      
      .fhList[[fhCurrent()]] <<- setGate(.fhList[[fhCurrent()]], bp)
      rv$FH <- .fhList[[fhCurrent()]]
    })

    fhPickPeaks <- observeEvent(input$pointPicker, {
      xPt <- nearPoints(fhHistData(.fhList[[fhCurrent()]]),
                        input$pointPicker, "xx", "intensity",
                        threshold = 25, maxpoints = 1)
      if(nrow(xPt) > 0){
        if(input$peakPicker == "A"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        xPt[1,1],
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        fhInit(.fhList[[fhCurrent()]])$Mc) 
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else if(input$peakPicker == "B"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        xPt[1,1],
                        fhInit(.fhList[[fhCurrent()]])$Mc)
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else if(input$peakPicker == "C"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        }
        rv$FH <- .fhList[[fhCurrent()]]
      }
    })

    ## Not sure why the following toggle events don't respond as
    ## eventReactives?
    fhToggleLinearity <- observeEvent(input$linearity, {
      .fhList[[fhCurrent()]] <<-
        updateFlowHist(.fhList[[fhCurrent()]],
                       linearity = input$linearity, analyze = TRUE)
      rv$FH <- .fhList[[fhCurrent()]]
    })
    
    fhToggleDebris <- observeEvent(input$debris, {
      .fhList[[fhCurrent()]] <<-
        updateFlowHist(.fhList[[fhCurrent()]],
                       debris = input$debris, analyze = TRUE)
      rv$FH <- .fhList[[fhCurrent()]]
    })

    fhToggleSamples <- observeEvent(input$sampSelect, {
      .fhList[[fhCurrent()]] <<-
        updateFlowHist(.fhList[[fhCurrent()]],
                         samples = input$sampSelect, analyze = TRUE)
      rv$FH <- .fhList[[fhCurrent()]]
    })
    
    observe({
      if(input$exit > 0){
        stopApp()
      }
      
    })

    output$fhHistogram <- renderPlot({
      plot(rv$FH, init = TRUE, nls = TRUE, comps = TRUE)
    })

    output$flowNumber <- renderText({
      paste("File", tags$b(fhCurrent()), "of", length(.fhList))
    })

    observe({
      dat <- gateData()
      updateSliderInput(session, "yrange", 
                        step = max(6, ceiling(log(max(dat[,2]))))/20,
                        max = max(6, ceiling(log(max(dat[,2])))),
                        value = 0, min = 0)
    })

    ##browser()                           
    gateData <- reactive({
      chan1 <- input$xcol
      chan2 <- input$ycol

      raw <<- exprs(fhRaw(.fhList[[fhCurrent()]]))
      
      df <- data.frame(x = raw[, chan1],
                       y = raw[, chan2] / raw[, chan1])
      names(df) <- c(chan1, paste(chan1, chan2, sep = "/"))
      df
    })

    ################
    ## Gate Plots ##
    ################

    gateMar <- c(5, 4, 1, 0)
    output$gatePlot <- renderPlot({
      forceChange <- input$setGate
      dat <- gateData()
      op = par(mar = gateMar)
      if(isGated(.fhList[[fhCurrent()]])){
      plot(dat, ylim = c(0, exp(log(max(dat[, 2])) - input$yrange)),
           type = 'n')
      points(dat[!fhGate(.fhList[[fhCurrent()]]), ], pch = 16,
             col = "#05050510")
      points(dat[fhGate(.fhList[[fhCurrent()]]), ], pch = 16,
             col = "#80000010") 
      } else
        plot(dat, ylim = c(0, exp(log(max(dat[, 2])) - input$yrange)),
             pch = 16, col = "#05050510")

      par(op)
    })

    output$gateResiduals <- renderPlot({
      op = par(mar = gateMar)
      plotResid(fhHistPlot(), main = "Gate Residuals", sub = "")
      par(op)
    })

    output$gatedData <- renderPlot({
      op = par(mar = gateMar)
      plot(fhHistPlot(), nls = FALSE, init = FALSE, comps = FALSE,
           main = "")
      par(op)
    })
    
    fhHistPlot <- reactive({
      ## returns the gated data, but doesn't set it permanently!
      dat <- gateData()
      bp <- brushedPoints(dat, xvar = names(dat)[1], yvar = names(dat)[2],
                          input$gatePlot_brush, allRows = TRUE)$selected_
      setGate(.fhList[[fhCurrent()]], bp, refresh = FALSE)
    })

  }
  runApp(shinyApp(ui = ui, server = server))
  return(.fhList)
}

selectPeaks <- function(fh, peakA, peakB, peakC){
  pA <- fhHistData(fh)[round(peakA, 0), c("xx", "intensity")]
  if(is.numeric(peakB))                 
    pB <- fhHistData(fh)[round(peakB, 0), c("xx", "intensity")]
  if(is.numeric(peakC))                 
    pC <- fhHistData(fh)[round(peakC, 0), c("xx", "intensity")]
  
  fh <- resetFlowHist(fh)

  if(is.numeric(peakC))
    newPeaks <- as.matrix(rbind(pA, pB, pC))
  else if(is.numeric(peakB))
    newPeaks <- as.matrix(rbind(pA, pB))
  else
    newPeaks <- as.matrix(rbind(pA))
  
  colnames(newPeaks) <- c("mean", "height")

  fhPeaks(fh) <- newPeaks
  
  fh <- addComponents(fh)
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)

  return(fh)
}
