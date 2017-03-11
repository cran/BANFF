#' Gene test statistics dataset
#'
#' The simulated test statistics on each gene node. With 10\% randomly selected as missing.
#'
#' @details We first simulate the underlying true feature labels by 25\% down-regulate, 50\% null and 25\% up-regulate. The pre-defined distribution for each regulation group is normal: the down-regulate genes are distributed as a normal with mean=-0.5 and sd=0.2, the null genes are distributed as a normal with mean=0 and sd=0.2, the up-regulated genes are distributed as a normal with mean=0.5 and sd=0.2. For the second step, we randomly knock out 10\% genes as missing nodes.
#' @docType data
#' @keywords datasets
#' @name test.stat
#' @aliases test.stat
#' @usage data(test.stat)
#' @format A vector of length=#genes, the value is the observed test statistics for each gene.
NULL
