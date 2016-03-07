library(shiny)
library(MASS)
ui = fluidPage(

    titlePanel('AMS 5 resource'),

    sidebarLayout(
        sidebarPanel(
            h1("Manual input"),
            textInput(inputId = "x", label = "x: List of numbers separated by comma (,)"),
            textInput(inputId = "y", label = "y: List of numbers separated by comma (,)"),
            br(), br(),
            h1("Random input (normal variables)"),
            numericInput(inputId = "rnormN", label = "N", value = 100),
            numericInput(inputId = "rnormXavg", label = "x Avg", value = 0),
            numericInput(inputId = "rnormXsd", label = "x SD", value = 1, min = 0,
                step = 0.1),
            numericInput(inputId = "rnormYavg", label = "y Avg", value = 0),
            numericInput(inputId = "rnormYsd", label = "y SD", value = 1, min = 0,
                step = 0.1),
            numericInput(inputId = "rnormCorr", label = "r (correlation)", value = 0,
                min = -1, max = 1, step = 0.01),
            actionButton(inputId = "rnormGo", "Draw N random normal variables for x and y")
            ),
        mainPanel(
            plotOutput(outputId = "scatter")
            plotOutput(outputId = "summary")
            )
        )
    )

server = function(input, output, session){

    sd = function(z) sqrt(sum((z-mean(z))^2) / length(z))
    cor = function(w, z) mean( (w - mean(w)) / sd(w) * (z - mean(z)) / sd(z) )
#       (mean(w*z) - mean(w)*mean(z)) / ( sd(w) * sd(z) )
#       mean( (w - mean(w)) * (z - mean(z))) / ( sd(w) * sd(z) )

    observe({
        input$rnormGo
        isolate({
            require(MASS)
            sigma = matrix(c(input$rnormXsd^2,
                input$rnormXsd * input$rnormYsd * input$rnormCorr,
                input$rnormXsd * input$rnormYsd * input$rnormCorr,
                input$rnormYsd^2), 2, 2)
            draws = mvrnorm(input$rnormN, c(input$rnormXavg, input$rnormYavg), sigma)
            updateTextInput(session, "x", value = 
                paste(as.character(round(draws[,1], 3)), collapse = ","))
            updateTextInput(session, "y", value = 
                paste(as.character(round(draws[,2], 3)), collapse = ","))
            })
        })

    output$scatter = renderPlot({

        ### Handle the input for x and y
        x = as.numeric(strsplit(gsub(" ", "", input$x), ",")[[1]])
        y = as.numeric(strsplit(gsub(" ", "", input$y), ",")[[1]])

        # Empty on both
        if (length(x) == 0 && length(y) == 0){
            x = 0
            y = 0
            }

        # If one list is totally empty, make the other all zeros (for bivariate plot)
        if (length(x) == 0 && length(y) > 0)
            x = rep(0, length(y))
        if (length(y) == 0 && length(x) > 0)
            y = rep(0, length(x))

        # Add NA's to make lengths of x and y match
        d = length(x) - length(y)
        if (d > 0)
            y = c(y, rep(NA, abs(d)))
        if (d < 0)
            x = c(x, rep(NA, abs(d)))

        # Regression lines
        mhat = cor(x, y) * sd(y) / sd(x)
        bhat = mean(y) - mhat*mean(x)

        layout(matrix(c(1,1,1,1,2,3),3,2,byrow=TRUE))
        plot(x, y, pch = 20, axes = FALSE)
        if (!is.na(mhat))
            abline(bhat, mhat, col = 'red')
        axis(1); axis(2)
        hist(x, freq = FALSE, col = 'gray')
        hist(y, freq = FALSE, col = 'gray')
        })


    }

shinyApp(ui = ui, server = server)
