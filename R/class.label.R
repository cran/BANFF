#' Gene class labels dataset
#'
#' The simulated class labels for each gene node, by merged community detection algorithm.
#'
#' @details Based on the simulated network structure. We firstly apply the fast community detection algorithm and then gradually merge the communities based on their pair-wised between edge counts, untill finally get three classes. The largest is assigned null, and then the upper/ down regulated are randomly chosen.
#' @docType data
#' @keywords datasets
#' @name class.label
#' @aliases class.label
#' @usage data(class.label)
#' @format A vector of length=#genes, the value is 1, 2, 3 as down-regulated/ null/ up-regulated class.
NULL
