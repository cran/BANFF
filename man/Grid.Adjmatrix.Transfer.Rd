\name{Grid.Adjmatrix.Transfer}
\alias{Grid.Adjmatrix.Transfer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Connector for Transferring Grid Coordinates into Adjacency Matrix
%%  ~~function to do ... ~~
}
\description{
This function is a supplementary function for transferring grid coordinates into adjacency matrix. The criteria for defining connected and disconnected: for any two vectors, if euclidean distance is less then a certain value, we consider these two vectors as connected, otherwise, disconnected. The purpose of this function is to quickly transform image data, such as brain image data, into adjacency matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Grid.Adjmatrix.Transfer(grid, euclidean.dist=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{a \code{x} by \code{y} matrix or data.frame, representing \code{y} vectors with dimension of \code{x} row binded (\code{rbind}) together.}
  \item{euclidean.dist}{a number representing the maximum euclidean distance to be considered as connected. The default setting is 1.}
%%     ~~Describe \code{x} here~~
}


\value{
a \code{y} by \code{y} binary (0/1) adjacency matrix, where \code{y} is the row of the argument \code{grid.}
}

\note{
Although in real world, the grid coordinates are less than 3-dimension, however, this function is still applicable for vectors larger than 3
%%  ~~further notes~~
}
\examples{
###10 3-D coordinates
x <- matrix(rnorm(3*10), nrow = 10)
###a 10x10 adjacency matrix returned for 'net'
###Distrance larger than 1 is considerded as connected
net<-Grid.Adjmatrix.Transfer(grid=x, euclidean.dist=1)


}

%% ~Make other sections like Warning with \section{Warning }{....} ~


