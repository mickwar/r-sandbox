library(shiny)
shinyUI(fluidPage(

    titlePanel('AMS 5 resource'),

    sidebarLayout(
        sidebarPanel(
            h1("Manuel input"),
            textInput(inputId = "x", label = "x: List of numbers separated by comma (,)"),
            textInput(inputId = "y", label = "y: List of numbers separated by comma (,)"),
            br(), br(),
            h1("Random input (normal variables)"),
            numericInput(inputId = "rnormN", label = "N", value = 100),
            numericInput(inputId = "rnormXavg", label = "x Avg", value = 0),
            numericInput(inputId = "rnormXsd", label = "x SD", value = 1, min = 0),
            numericInput(inputId = "rnormYavg", label = "y Avg", value = 0),
            numericInput(inputId = "rnormYsd", label = "y SD", value = 1, min = 0),
            numericInput(inputId = "rnormCorr", label = "r (correlation)", value = 0,
                min = -1, max = 1),
            actionButton(inputId = "rnormGo", "Draw N random normal variables for x and y")
            ),
        mainPanel(
            plotOutput(outputId = "scatter")
            )
        )
    ))
