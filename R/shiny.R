#' @importFrom shiny fluidPage nearPoints reactive radioButtons
#'   actionButton plotOutput reactive eventReactive shinyApp titlePanel
#'   sidebarLayout sidebarPanel htmlOutput fluidRow tags mainPanel
#'   renderPrint renderTable renderPlot renderText column observe runApp
#'   stopApp wellPanel updateRadioButtons
NULL

#' @importFrom utils str
NULL

#' Visually assess histogram fits
#'
#' Visually assess histogram fits
#' @title browseFlowHist
#' @param flowList a list of \code{\link{FlowHist}} objects
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
  .fhI <- 1
  .fhList <- flowList

  initialLinearity <- fhLinearity(.fhList[[1]])
  if(initialLinearity == "fixed")
    initialLinearity <- "ON"
  else
    initialLinearity <- "OFF"

  if(debug) message("init Linearity: ", initialLinearity)
  
  initialDebris <- fhDebris(.fhList[[1]])

  if(debug) message("init Debris: ", initialDebris)
  
  ui <- fluidPage(
    ## sidebarLayout(
    fluidRow(
      column(2,
             fluidRow(
               wellPanel(
                 htmlOutput("flowNumber", align = "center"),
                 tags$br(),
                 fluidRow(
                   column(6, actionButton("prev", label = "Prev")),
                   column(6, actionButton("nxt", label = "Next")))))
             ),
      ## tags$hr(),
      column(10,
             fluidRow(
               column(3, 
                      wellPanel(radioButtons(inputId = "peakPicker",
                                             label = "Move peak:", 
                                             choices = list("A" = "A",
                                                            "B" = "B"), 
                                             selected = "A", inline =
                                                               TRUE))), 
               column(4,
                      wellPanel(radioButtons(inputId = "linearity",
                                             label = "Linearity",
                                             choices =
                                               list("Fixed" = "ON",
                                                    "Variable" = "OFF"), 
                                             inline = TRUE,
                                             selected =
                                               initialLinearity))),
               column(3,
                      wellPanel(radioButtons(inputId = "debris",
                                             label = "Debris Model",
                                             choices =
                                               list("MC" = "MC",
                                                    "SC" = "SC"),  
                                             inline = TRUE,
                                             selected =
                                               initialDebris))),
               column(2,
                      tags$br(),
                      actionButton("exit", label = "Return to R")
                      )))),
        ## tags$hr(),
    fluidRow(
      plotOutput("init", click = "pointPicker"))
  )

  server <- function(input, output, session){
    nxtVal <- 0
    prevVal <- 0
    prefix <- ""
    
    fhInitPlot <- reactive({
      if(debug){
        message(prefix, "fhInitPlot ",
                        environmentName(environment())) 
        prefix <<- paste(prefix, " ", sep = "")}
      tmp <- .fhList[[fhCurrent()]]
      if(debug) message(prefix, "fh@linearity = ", fhLinearity(tmp))
      if(debug) message(prefix, "button value = ", input$linearity)
      xPt <- nearPoints(fhHistData(.fhList[[fhCurrent()]]),
                        input$pointPicker, "xx", "intensity",
                        threshold = 25, maxpoints = 1)
      if(debug) message(prefix, "xPt set")
      if(nrow(xPt) > 0){
        if(debug) message(prefix, "peak Picker")
        if(input$peakPicker == "A"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]], xPt[1,1],
                        fhInit(.fhList[[fhCurrent()]])$Mb) 
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else {
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma, xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        }
      }

      if(debug) message(prefix, "checking linearity for element ", .fhI)
      if(input$linearity == "ON" &&
         fhLinearity(.fhList[[fhCurrent()]]) != "fixed")
      {
        if(debug) message(prefix, "fixing linearity")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         linearity = "fixed", analyze = TRUE)
      }
      else 
        if(input$linearity == "OFF" &&
           fhLinearity(.fhList[[fhCurrent()]]) != "variable") { 
          if(debug) message(prefix, "modeling linearity")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           linearity = "variable", analyze = TRUE)
        }

      if(debug) message(prefix, "input$debris: ", input$debris)
      if(input$debris == "SC" &&
         fhDebris(.fhList[[fhCurrent()]]) != "SC")
      {
        if(debug) message(prefix, "switching to SC")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         debris = "SC", analyze = TRUE)
      }
      else 
        if(input$debris == "MC" &&
           fhDebris(.fhList[[fhCurrent()]]) != "MC") { 
          if(debug) message(prefix, "switching to MC")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           debris = "MC", analyze = TRUE)
        }
      if(debug){
        prefix <<- substring(prefix, 2)
        message(prefix, "returning from fhInitPlot")
      }
      
      .fhList[[fhCurrent()]]
    })

    fhCurrent <- reactive({
      if(debug) {
        message(prefix, "fhCurrent, starting at ", .fhI)
        prefix <<- paste(prefix, " ", sep = "")}
      if(input$nxt > nxtVal){
        if(debug) message(prefix, "moving forwards to ", .fhI + 1)
        nxtVal <<- input$nxt
        if(.fhI < length(.fhList))
          .fhI <<- .fhI + 1

        if(debug) message(prefix, "updating linval forward")
        if(fhLinearity(.fhList[[.fhI]]) == "fixed"){
          if(debug) message(prefix, "turning button ON")
          linVal <- "ON"
        } else {
          if(debug) message(prefix, "turning button OFF")
          linVal <- "OFF"
        }
        
        updateRadioButtons(session, "linearity",
                           selected = linVal)

      }

      if(input$prev > prevVal){
        if(debug){
          message(prefix, "moving backwards to ", .fhI - 1)
        }
        prevVal <<- input$prev
        if(.fhI > 1)
          .fhI <<- .fhI - 1

        if(debug) message(prefix, "updating linval backwards")
        if(fhLinearity(.fhList[[.fhI]]) == "fixed"){
          if(debug) message(prefix, "turning button ON")
          linVal <- "ON"
        } else {
          if(debug) message(prefix, "turning button OFF")
          linVal <- "OFF"
        }
        updateRadioButtons(session, "linearity",
                           selected = linVal)
      }

      if(fhDebris(.fhList[[.fhI]]) == "SC" ){
        if(debug) message(prefix, "Switching to SC")
        debVal <- "SC"
      } else {
        if(debug) message(prefix, "Switching to MC")
        debVal <- "MC"
      }
      
      updateRadioButtons(session, "debris",
                         selected = debVal)
      
      if(debug){
              prefix <<- substring(prefix, 2)
              message(prefix, "returning from fhCurrent")
      }
      .fhI
    })
    
    observe({
      if(input$exit > 0){
        stopApp()
      }
      
    })

    output$init <- renderPlot({
      if(debug){
        message(prefix, "renderPlot")
        prefix <<- paste(prefix, " ", sep = "")
      }
      plot(fhInitPlot(), init = TRUE, nls = TRUE, comps = TRUE)
      if(debug){
        prefix <<- substring(prefix, 2)
        message(prefix, "returning from renderPlot")
      }
    })

    output$flowNumber <- renderText({
      if(debug) message(prefix, "flowNumber")
      paste("File", tags$b(fhCurrent()), "of", length(.fhList))
    })

  }
  runApp(shinyApp(ui = ui, server = server))
  return(.fhList)
}

selectPeaks <- function(fh, peakA, peakB){
  pA <- fhHistData(fh)[round(peakA, 0), c("xx", "intensity")]
  if(is.numeric(peakB))                 
    pB <- fhHistData(fh)[round(peakB, 0), c("xx", "intensity")]
  
  fh <- resetFlowHist(fh)

  if(is.numeric(peakB))
    newPeaks <- as.matrix(rbind(pA, pB))
  else
    newPeaks <- as.matrix(rbind(pA))
  
  colnames(newPeaks) <- c("mean", "height")

  fhPeaks(fh) <- newPeaks
  
  fh <- addComponents(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)

  return(fh)
}
