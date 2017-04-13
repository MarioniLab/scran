exploreData <- function(x, cell.data, gene.data, red.dim, run=TRUE) 
# This function creates an interactive shiny app to explore expression data.
#
# created by Karsten Bach
# written 11 April 2017
# last modified 11 April 2017
{
    # Check on gene and cell names.
    xr <- rownames(x)
    gd <- rownames(gene.data)
    if (is.null(xr) && is.null(gd)) { 
        xr <- gd <- paste0("Gene", seq_len(nrow(x)))
    } else if (!identical(xr, gd)) {
        stop("'x' and 'gene.data' should have identical rownames")
    } else if (is.null(xr)) { 
        xr <- gd        
    } else {
        gd <- xr
    }

    xc <- colnames(x)
    cd <- rownames(cell.data)
    if (is.null(xc) && is.null(cd)) { 
        xc <- cd <- paste0("Cell", seq_len(ncol(x)))
    } else if (!identical(xc, cd)) {
        stop("colnames of 'x' and rownames of 'cell.data' should be identical")
    } else if (is.null(xc)) { 
        xc <- cd        
    } else {
        cd <- xc
    }

    # Set names 
    rownames(x) <- xr
    rownames(gene.data) <- gd
    colnames(x) <- xc
    rownames(cell.data) <- cd
    covariates <- colnames(cell.data)
    plot.data <- data.frame(Dim1=red.dim[,1], Dim2=red.dim[,2])

    # Set up shiny UI
    ui <- fluidPage(
        titlePanel("Data exploration"),
        sidebarLayout(
            sidebarPanel(
                inputPanel(
                    selectInput("colorBy", label = "Color by", choices=covariates, selected=covariates[1]),
                    selectInput("groupBy", label = "Group by", choices=covariates, selected=covariates[1]) 
                ),
                plotOutput("distPlot1")
            ),
            mainPanel(
                tabsetPanel(tabPanel("Plots",
                    splitLayout(plotOutput("tSNE"), plotOutput("tSNE2")),
                    dataTableOutput("table"))
                )
            )
        )
    )

    server <- function(input, output) {
        # Load the gene level data
        output$table <- renderDataTable({
            datatable(gene.data, filter="top", selection=list(mode="single", selected=1))
        })

        # tSNE plot colored by covariates
        output$tSNE <- renderPlot({
            plot.data.copy <- data.frame(plot.data, Covariate=cell.data[,input$colorBy])
            tsnPlot <- ggplot(plot.data.copy, aes_string(x="Dim1", y="Dim2", color="Covariate")) +
                geom_point(size=1.5) +
                labs(color=input$colorBy) +
                theme_void()
            tsnPlot
        })

        # tSNE plot colored by gene expression
        output$tSNE2 <- renderPlot({
            s <- input$table_rows_selected
            if (!is.null(s)) {
                gene <- rownames(gene.data)[s]
                plot.data.copy <- data.frame(plot.data, Expression=x[gene,])
                plot.data.copy <- plot.data.copy[base::order(plot.data.copy$Expression),]
                tsnPlot <- ggplot(plot.data.copy, aes_string(x="Dim1", y="Dim2", color="Expression")) +
                    geom_point(size=1.5) +
                    scale_color_viridis() +
                    ggtitle(gene) +
                    theme_void()
                tsnPlot
            }
        })

        # Distribution of gene expression by covariate
        output$distPlot1 <- renderPlot({
            s <- input$table_rows_selected
            if (!is.null(s)) {
                gene <- rownames(gene.data)[s]
                plot.data.copy <- data.frame(Covariate=cell.data[,input$groupBy], Expression=x[gene,])
                plt <- ggplot(plot.data.copy, aes_string(x="Covariate", y="Expression")) +
                    geom_boxplot() +
                    geom_point(position="jitter", alpha=0.2, shape=19) +
                    ggtitle(gene) +
                    xlab(input$groupBy) +
                    theme_bw()
                plt
            }
        })
    }

    app <- shinyApp(ui, server)
    if (run) {
        return(runApp(app))
    } else {
        return(app)
    }
}

