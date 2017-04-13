\name{Explore Data}
\alias{exploreData}

\title{Shiny app for explorative data analysis}
\description{Generate an interactive Shiny app to explore gene expression patterns in single-cell data}

\usage{
exploreData(x, cell.data, gene.data, red.dim, run=TRUE)
}

\arguments{
\item{x}{A numeric matrix of expression values to be visualized.}
\item{cell.data}{A data frame of cell information, where each row corresponds to a column of \code{x}.}
\item{gene.data}{A data frame of gene information, where each row corresponds to a row of \code{x}.}
\item{red.dim}{A numeric matrix with two colums, specifying the reduced-dimension coordinates for each cell.}
\item{run}{A logical scalar specifying whether the app should be run immediately.}
}

\details{
This function will return a Shiny app object that can be run with \code{\link{runApp}}.
The app allows the user to interactively explore gene expression patterns in single-cell RNA-seq data. 
Explorative analysis is focused on comparing gene exression between different groups of cells, as defined by the covariates of \code{cell.data}.

Three plots are shown in the app: 
\itemize{
\item a scatterplot of cell locations based on the \code{red.dim} coordinates, colored by a selected covariate
\item a scatterplot of cell locations based on the \code{red.dim} coordinates, colored by expression of a selected gene 
\item boxplot(s) of expression values for a selected gene, grouped by a selected covariate.
}

Several options are available within the app:
\describe{
\item{\dQuote{Color by}:}{Covariate to be used for coloring the first scatter plot.}
\item{\dQuote{Group by}:}{Covariate with which expression values are grouped in the boxplots.}
}
In addition, the \code{gene.data} data frame is rendered as an interactive table using the JavaScript library \code{DataTable}.
This allows the user to subset/search the feature data and select a gene by clicking on the corresponding row.
}

\value{
If \code{run=FALSE}, a Shiny app object is returned, which can be run with \code{\link{runApp}}.
If \code{run=TRUE}, a Shiny app object is created and run.
}

\seealso{
\code{\link{runApp}},
}

\author{
Karsten Bach
}

\examples{
# Set up example data
example(newSCESet)
x <- exprs(example_sceset)
cell.data <- pData(example_sceset)
gene.data <- data.frame("meanExp"=rowMeans(x))
rownames(gene.data) <- rownames(x)

# Mocking up  some reduced dimensions.
library(Rtsne)
tsn <- Rtsne(t(x), perplexity=10)
red.dim <- tsn$Y[,1:2]

# Creating the app object.
app <- exploreData(x, cell.data, gene.data, red.dim, run=FALSE)
if (interactive()) { shiny::runApp(app) }

\dontrun{# Running directly from the function.
saved <- exploreData(x, cell.data, gene.data, red.dim)
}
}

