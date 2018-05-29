#' @importFrom shiny actionButton brushOpts brushedPoints column
#'   eventReactive eventReactive fluidPage fluidRow HTML htmlOutput
#'   mainPanel nearPoints numericInput observe observeEvent plotOutput
#'   radioButtons reactive reactiveValues renderPlot renderPrint
#'   renderTable renderText runApp selectInput shinyApp sidebarLayout
#'   sidebarPanel sliderInput stopApp tags titlePanel updateNumericInput
#'   updateRadioButtons updateSelectInput updateSliderInput wellPanel
#'   tabsetPanel tabPanel textInput
NULL

#' @importFrom utils str
NULL

#' Visually assess and correct histogram fits
#'
#' Visually assess histogram fits, correcting initial values, and selecting
#' model components.
#'
#' This function will open a browser tab displaying the first
#' \code{\link{FlowHist}} object from the argument \code{flowList}. Using
#' the interface, the user can modify the starting values for the histogram
#' peaks, select different debris model components, toggle the linearity
#' option, select which peak to treat as the standard, and, if multiple
#' standard sizes are available, select which one to apply.
#'
#' See the "Getting Started" vignette for a tutorial introduction.
#' 
#' @title browseFlowHist
#' @param flowList either a \code{\link{FlowHist}} object, or a list of
#'   \code{\link{FlowHist}} objects
#' @param debug boolean, turns on debugging messages
#' @return Returns the list of \code{\link{FlowHist}} objects, updated by
#'   any changes made in the GUI.
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
    ## if flowList is a single FlowHist object, wrap it in a list so we can
    ## use the same code for processing single and multiple FlowHist
    ## objects. Convert back to a single object on return.
    flowList <- list(flowList)
    names(flowList) <- fhFile(flowList[[1]])
  }
  
  .fhI <- 1
  .fhList <- flowList

  initialLinearity <- fhLinearity(.fhList[[1]])

  if(debug) message("init Linearity: ", initialLinearity)

  initialStdSelected <- fhStdSelected(.fhList[[1]])
  if(debug) message("initialStdSelected: ", initialStdSelected)
  standardList <- fhStdSizes(.fhList[[1]])
  if(debug) message("standardList: ", standardList)
  initialStdPeak <- fhStdPeak(.fhList[[1]])
  
  initialDebris <- fhDebris(.fhList[[1]])

  initialSamples <- fhSamples(.fhList[[1]])

  initialAnnotation <- fhAnnotation(.fhList[[1]])
  
  if(debug) message("init Debris: ", initialDebris)

  raw <- exprs(fhRaw(.fhList[[.fhI]]))
  chan1 <- fhChannel(.fhList[[.fhI]])
  if(length(viewFlowChannels(.fhList[[.fhI]])) > 1){
    ## why 2? sometimes, by chance, the second column is the SS value.
    chan2 <- viewFlowChannels(.fhList[[.fhI]])[2]
  } else {
    ## there's only one channel, so gating is pointless
    chan2 <- chan1
  }
  
  initGateData <- data.frame(x = raw[, chan1],
                             y = raw[, chan2] / raw[, chan1]) 

  maxY <- max(initGateData$y[!is.infinite(initGateData$y)],
              na.rm = TRUE)

  if(is.nan(maxY))
    maxY <- 6
  
  ui <- fluidPage(
    tags$head(
           tags$style(HTML("
      .sidepanel {
          max-width: 300px;
      }
      #gatePlot, #gatedData, #gateResiduals {
          max-height: 300px;
      }
      #exit {
          margin-top: 1.4em;
      }
      #setGate {
          margin-top: 1.8em;
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
                 column(6,
                        ## If we use numericInput, the user is able to
                        ## enter invalid values in the textbox. With
                        ## selectInput we constrain the possible input
                        ## values to those that won't break anything.
                        selectInput('sampSelect', 'Samples',
                                    initialSamples,
                                    choices = list("1" = 1,
                                                   "2" = 2,
                                                   "3" = 3,
                                                   "4" = 4,
                                                   "5" = 5,
                                                   "6" = 6))), 
                 column(6,
                        selectInput(inputId = "peakPicker",
                                    label = "Peak",
                                    selected = "A",
                                    choices = list("A", "B", "C", "D", "E",
                                                   "F")))),
               fluidRow(
                 column(6, 
                        selectInput(inputId = 'standardSelect',
                                    label = 'Standard Value',
                                    choices = standardList,
                                    selected = initialStdSelected)),
                 column(6,
                        selectInput(inputId = "standardPeak",
                                    label = "Standard Peak",
                                    selected = initialStdPeak,
                                    choices = list("X", "A", "B")))), 
               fluidRow(
                 column(6,
                        selectInput(inputId = "linearity",
                                    label = "Linearity",
                                    selected = initialLinearity,
                                    choices = list("Fixed" = "fixed",
                                                   "Variable" =
                                                     "variable"))),  
                 column(6, 
                        selectInput(inputId = "debris",
                                    label = "Debris Model",
                                    choices = list("MC" = "MC",
                                                   "SC" = "SC", 
                                                   "none" = "none"),  
                                    selected = initialDebris)))
             ))),
      column(width = 9,
             plotOutput("fhHistogram", click = "pointPicker"))
    ),
    tabsetPanel(type = "tabs",
      tabPanel("Annotation",
               textInput("annotation", "",
                         initialAnnotation, width = "100%"), 
               actionButton("updateAnnotation", label = "Annotate")),
      tabPanel("Gating", 
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
                             max(6, ceiling(log(maxY)))/20, 
                           max = max(6, ceiling(log(maxY)), na.rm = TRUE),
                           value = 0, dragRange = FALSE),
               fluidRow(
                 column(6, 
                        selectInput('yType', 'Y axis', c("Y/X", "Y"),
                                    selected = "Y/X")),
                 column(6,
                        actionButton("setGate", label = "Set Gate"))
               )))),
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
    )
    
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
      
      updateSelectInput(session, "linearity",
                         selected = fhLinearity(.fhList[[rv$fhI]]))
      updateSelectInput(session, "debris",
                         selected = fhDebris(.fhList[[rv$fhI]]))
      updateNumericInput(session, "sampSelect",
                         value = fhSamples(.fhList[[rv$fhI]]))

      updateSelectInput(session, "standardSelect",
                         selected = fhStdSelected(.fhList[[rv$fhI]]))
      updateSelectInput(session, "standardPeak",
                        selected = fhStdPeak(.fhList[[rv$fhI]]))
      updateTextInput(session, "annotation",
                      value = fhAnnotation(.fhList[[rv$fhI]]))
      rv$FH <- .fhList[[rv$fhI]]
      rv$fhI
    })      

    fhSetGate <- observeEvent(input$setGate, { 
      dat <- gateData()
      bp <- brushedPoints(dat, xvar = names(dat)[1],
                          yvar = names(dat)[2], input$gatePlot_brush,
                          allRows = TRUE)$selected_
      if(sum(bp) == 0) bp <- logical()
      
      .fhList[[fhCurrent()]] <<- setGate(.fhList[[fhCurrent()]], bp)
      rv$FH <- .fhList[[fhCurrent()]]
    })

    fhUpdateAnnotation <- observeEvent(input$updateAnnotation, {
      fhAnnotation(.fhList[[fhCurrent()]]) <<- input$annotation
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
        } else if(input$peakPicker == "D"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        fhInit(.fhList[[fhCurrent()]])$Mc,
                        xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else if(input$peakPicker == "E"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        fhInit(.fhList[[fhCurrent()]])$Mc,
                        fhInit(.fhList[[fhCurrent()]])$Md,
                        xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else if(input$peakPicker == "F"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        fhInit(.fhList[[fhCurrent()]])$Mc,
                        fhInit(.fhList[[fhCurrent()]])$Md,
                        fhInit(.fhList[[fhCurrent()]])$Me,
                        xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        }

        rv$FH <- .fhList[[fhCurrent()]]
      }
    })

    ## Not sure why the following toggle events don't respond as
    ## eventReactives?
    fhToggleLinearity <- observeEvent(input$linearity, {
      if(fhLinearity(.fhList[[fhCurrent()]]) != input$linearity){
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         linearity = input$linearity, analyze = TRUE)
        rv$FH <- .fhList[[fhCurrent()]]
      }
    })
    
    fhToggleDebris <- observeEvent(input$debris, {
      if(fhDebris(.fhList[[fhCurrent()]]) != input$debris){
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         debris = input$debris, analyze = TRUE)
        rv$FH <- .fhList[[fhCurrent()]]
      }
    })

    fhToggleSamples <- observeEvent(input$sampSelect, {
      if(fhSamples(.fhList[[fhCurrent()]]) != input$sampSelect){
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         samples = input$sampSelect, analyze = TRUE)
        rv$FH <- .fhList[[fhCurrent()]]
      }
    })

    fhUpdateStdPeak <- observeEvent(input$standardPeak, {
      if(fhStdPeak(.fhList[[fhCurrent()]]) != input$standardPeak){
        fhStdPeak(.fhList[[fhCurrent()]]) <<- input$standardPeak
        rv$FH <- .fhList[[fhCurrent()]]
      }
    })

    fhUpdateStdSelected <- observeEvent(input$standardSelect, {
      ## don't update the fh object if the input is the same as the actual
      ## value (may come up when switching to a new object)
      if(fhStdSelected(.fhList[[fhCurrent()]]) != input$standardSelect){
        fhStdSelected(.fhList[[fhCurrent()]]) <<-
          as.numeric(input$standardSelect)
        rv$FH <- .fhList[[fhCurrent()]]
      }
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
      maxY <- max(dat[,2][!is.nan(dat[,2]) & is.finite(dat[,2])],
                  na.rm = TRUE)
      updateSliderInput(session, "yrange", 
                        step = max(6, ceiling(log(maxY)),
                                   na.rm = TRUE )/20,
                        max = max(6, ceiling(log(maxY)), na.rm = TRUE),
                        value = 0, min = 0)
    })

    ##browser()                           
    gateData <- reactive({
      chan1 <- input$xcol
      chan2 <- input$ycol

      raw <<- exprs(fhRaw(.fhList[[fhCurrent()]]))
      if(input$yType == "Y/X"){
        yvals <- raw[, chan2] / raw[, chan1]
        yName <- paste(chan2, chan1, sep = "/")
      } else if(input$yType == "Y"){
        yvals <- raw[, chan2]
        yName <- chan2
      }
      df <- data.frame(x = raw[, chan1],
                       y = yvals)
      names(df) <- c(chan1, yName)
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
      ## need to account for infinite values when setting plot ranges.
      ## Infinite values generated when the dat[, 1] contains 0 values. 
      maxY <- max(dat[!is.nan(dat[,2]) & is.finite(dat[,2]), 2],
                  na.rm = TRUE)
      plot(dat, ylim = c(0, exp(log(maxY) - input$yrange)),
           type = 'n')
      if(isGated(.fhList[[fhCurrent()]])){
        points(dat[!fhGate(.fhList[[fhCurrent()]]), ], pch = 16,
               col = "#05050510")
        points(dat[fhGate(.fhList[[fhCurrent()]]), ], pch = 16,
               col = "#80000010") 
      } else
        points(dat, pch = 16, col = "#05050510")

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
  if(length(.fhList) == 1){
    ## if we started with one flowHist object, return one flowHist object
    return(.fhList[[1]])
  } else {
    ## we started with a list of flowHist objects, return a list of
    ## modified objects:
    return(.fhList)
  }
}

