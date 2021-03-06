% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.label.R
\docType{data}
\name{class.label}
\alias{class.label}
\title{Gene class labels dataset}
\format{A vector of length=#genes, the value is 1, 2, 3 as down-regulated/ null/ up-regulated class.}
\usage{
data(class.label)
}
\description{
The simulated class labels for each gene node, by merged community detection algorithm.
}
\details{
Based on the simulated network structure. We firstly apply the fast community detection algorithm and then gradually merge the communities based on their pair-wised between edge counts, untill finally get three classes. The largest is assigned null, and then the upper/ down regulated are randomly chosen.
}
\keyword{datasets}
