\name{HyperPara.Select}
\alias{HyperPara.Select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selecting Hyper Parameters by Bayesian Model Averaging 
%%  ~~function to do ... ~~
}
\description{
This function is aiming at selecting appropriate hyper parameters by Bayesian Model Averaging. This process of this function has already embedded in function \code{Networks.STD()} and \code{Networks.Fast()}. The purpose of this function is for having the selected hyper-parameters independently. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HyperPara.Select(net,pvalue,piall,rhoall,n=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{an "n" by "n" binary adjacent matrix (0/1) for the network configuration with n=length(pvalue)}
  \item{pvalue}{a vector of p-values obtained from large scale statistical hypothesis testing}
  \item{piall}{a vector of possible choices for "pi0" in an increasing order}

\item{rhoall}{a vector of possible choices of "rho0" and "rho1" in an increasing order}
\item{n}{a number of iterations you set for Bayesian Model Averaging. The default setting is 30, which is accord with the embedded inner process of in function \code{Networks.STD()} and \code{Networks.Fast()}. More iterations could be set if you have a desire to have more accurate parameter}

}

\value{
a list of selected parameter
\describe{
\item{pi0}{the selected pi0}
\item{rho0}{the selected rho0}
\item{rho1}{the selected rho1}
}
}
\references{
Yize Zhao, Jian Kang, Tianwei Yu (2014) A Bayesian nonparametric model for selecting gene and gene sub network, Annals of Applied Statistics, in press. 

Zhou Lan, Jian Kang, Tianwei Yu, Yize Zhao, BANFF: an R package for network identifications via Bayesian nonparametric mixture models, working paper.
}
\author{
Zhou Lan, Jian Kang, Tianwei Yu and Yize Zhao


Department of Biostatistics and Bioinformatics, Emory University}

