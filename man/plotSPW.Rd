% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSPW.R
\name{plotSPW}
\alias{plotSPW}
\title{Plot subpathway network graph}
\usage{
plotSPW(
  subpathwayID,
  layout = NULL,
  margin = 0,
  vertex.label.cex = 0.8,
  vertex.label.font = 0.8,
  vertex.label.dist = 1,
  vertex.size = 12,
  edge.arrow.width = 3,
  edge.label.cex = 0.6,
  vertex.label.color = "black",
  vertex.color = "#F08080",
  vertex.frame.color = "dimgray",
  edge.color = "grey70",
  edge.label.color = "dimgray",
  sub = NULL,
  main = NULL
)
}
\arguments{
\item{subpathwayID}{Subpathway id .A character vector.}

\item{layout}{A matrix of x-y coordinates with two dims.
Determine the placement of the nodes for drawing a graph.}

\item{margin}{A numeric. The value is usually between -0.5 and 0.5, which is able to zoom in or out a subpathway graph. The default is 0.}

\item{vertex.label.cex}{A numeric vector of node label size.}

\item{vertex.label.font}{A numeric vector of label font.}

\item{vertex.label.dist}{A numeric vector of label dist.}

\item{vertex.size}{A numeric vector of Node size. See plot.igraph.}

\item{edge.arrow.width}{Edge arrow width. The default is 3.}

\item{edge.label.cex}{Edge label size.}

\item{vertex.label.color}{A vector of node label colors. The default is black.}

\item{vertex.color}{A vector of node colors. The default is the KEGG node color.}

\item{vertex.frame.color}{A vector of node frame color. The default is dimgray.}

\item{edge.color}{A vector of edge color. The default is dimgray}

\item{edge.label.color}{A vector of edge label color. The default is dimgray.}

\item{sub}{A character string of subtitle.}

\item{main}{A character string of main title.}
}
\value{
a subpathway map
}
\description{
The function plotSPW can visualize subpathway network graph.
}
\examples{
# load depend package
library(igraph)
# plot network graph of the subpathway "00020_4"
plotSPW("00020_4")
}
