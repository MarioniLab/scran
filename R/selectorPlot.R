selectorPlot <- function(x, y, persist=FALSE, plot.width=5, plot.height=500, pch=16, ...)
# This generates an interactive plot where x/y coordinates can be selected
# for further examination within R.
#
# created by Aaron Lun
# written 8 June 2016
{
    N <- length(x)
    if (length(y)!=N) { stop("length of x and y vectors should be equal") }
    collected <- new.env()
    resetValues(collected, data.frame(x=x, y=y))

    # Internal functions, to avoid passing many arguments around.
    plotFun1 <- function(output) {
        generatePlot1(output, collected, pch=pch, ...)
    }
    plotFun2 <- function(output) {
        generatePlot2(output, collected, ...)
    }
    updateSelect <- function(input, output, setting) {
        brushed <- brushedPoints(collected$coords, xvar="x", yvar="y", input$plot1_brush, allRows=TRUE)
        collected$current.selected[brushed$selected_] <- setting
        plotFun1(output)
    }

    # Generating the page layout.
    ui <- fluidPage(
        fluidRow(
            column(width = plot.width,
                plotOutput("plot1", height = plot.height,
                    brush = brushOpts(id = "plot1_brush")
                )
            ),
            column(width = plot.width,
                plotOutput("plot2", height = plot.height)
            )
        ),
        actionButton("select", "Select"),
        actionButton("unselect", "Deselect"),
        actionButton("clear", "Clear selection"),
        actionButton("list_add", "Add to list"),
        actionButton("reset", "Reset all"),
        actionButton("finish", "Save list to R")
    )

    # Setting up the server actions.
    server <- function(input, output) {
        plotFun1(output)
        plotFun2(output)

        observeEvent(input$select, { updateSelect(input, output, TRUE) })
        observeEvent(input$unselect, { updateSelect(input, output, FALSE) })
        observeEvent(input$clear, {
            collected$current.selected[] <- FALSE
            plotFun1(output)
        })

        observeEvent(input$list_add, {
            n <- length(collected$all.selected)
            collected$all.selected[[n+1]] <- collected$current.selected
            collected$old.selected[collected$current.selected] <- TRUE
            collected$current.selected[] <- FALSE
            plotFun1(output)
            plotFun2(output)
        })

        observeEvent(input$reset, {
            resetValues(collected)
            plotFun1(output)
            plotFun2(output)
        })

        observeEvent(input$finish, {
            keep.val <- collected$all.selected
            if (!persist) { resetValues(collected) }
            stopApp(keep.val)
        })
    }

    shinyApp(ui, server)
}

# A battery of internal functions, taken out to reduce the length of the main function.

resetValues <- function(collected, coords=NULL) {
    if (is.null(coords)) { coords <- collected$coords }
    else { collected$coords <- coords }
    N <- nrow(coords)
    collected$old.selected <- logical(N)
    collected$current.selected <- logical(N)
    collected$all.selected <- list()
}

generatePlot1 <- function(output, collected, ...) {
    output$plot1 <- renderPlot({
        cols <- rep("grey", length(collected$current.selected))
        cols[collected$old.selected] <- "orange"
        cols[collected$current.selected] <- "red"
        x <- collected$coords$x
        y <- collected$coords$y
        plot(x, y, col=cols, ...)
    })
}

generatePlot2 <- function(output, collected, ...) {
    output$plot2 <- renderPlot({
        n <- length(collected$all.selected)
        shading <- rev(grey.colors(n))
        x <- collected$coords$x
        y <- collected$coords$y
        plot(x, y, type="n", ...)
        for (i in seq_len(n)) {
            current <- collected$all.selected[[i]]
            text(x[current], y[current], labels=i, col=shading[i])
        }
    })
}

