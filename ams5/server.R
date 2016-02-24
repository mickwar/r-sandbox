library(shiny)
shinyServer(function(input, output, session){


#   observeEvent(input$rnormGo, {
#       updateTextInput(session, "x", value = 
#           paste(as.character(round(rnorm(input$rnormN, input$Xavg, input$Xsd), 3)),
#               collapse = ","))
#       updateTextInput(session, "y", value = 
#           paste(as.character(round(rnorm(input$rnormN, input$Yavg, input$Ysd), 3)),
#               collapse = ","))
#       })

    output$scatter = renderPlot({

        input$rnormGo

        isolate({
            updateTextInput(session, "x", value = 
                paste(as.character(round(rnorm(input$rnormN, input$Xavg, input$Xsd), 3)),
                    collapse = ","))
#           updateTextInput(session, "y", value = 
#               paste(as.character(round(rnorm(input$rnormN, input$Yavg, input$Ysd), 3)),
#                   collapse = ","))
            })

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

        layout(matrix(c(1,1,1,1,2,3),3,2,byrow=TRUE))
        plot(x, y, pch = 20, axes = FALSE)
        axis(1); axis(2)
        hist(x, freq = FALSE, col = 'gray')
        hist(y, freq = FALSE, col = 'gray')
        })
    })
