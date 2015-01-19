\name{summary.Networks.Fast}
\alias{summary.Networks.Fast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing Bayesian Network Discovery using a Hybrid Fast Algorithm

%%  ~~function to do ... ~~
}
\description{Summary method for class \code{"Networks.Fast"}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{summary}{Networks.Fast}(object, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"Networks.Fast"} resulting of a call to \code{\link{Networks.Fast}}.
  
%%     ~~Describe \code{x} here~~
}
\item{\dots}{other arguments}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Networks.Fast}}.}
\examples{
library(igraph)
####Gene Network discovery
##Generating Scale free Gene Network
library(igraph)
g <- barabasi.game(50, power=1, zero.appeal=1.5,directed = FALSE)
net=as(get.adjacency(g,attr=NULL),"matrix")
##Random assign selected genes and make the signal intension as gaussian mixture
newz=rep(c(1,0,0,1,0),10)
Simnorm=function(n){
weight = c(0.4, 0.6)
mu = c(8,6)
sigma = c(1,0.5)
z = sample(c(1,2),size=n, prob=weight,replace=TRUE)
r = rnorm(n,mean=mu[z],sd=sigma[z])
return(r)
}
testcov<-0
for(i in 1:50){
 if(newz[i]==0){
   testcov[i]<-rnorm(1,mean=0,sd=1)
 }else{
  testcov[i]<-Simnorm(1)
 }
}
pvalue=pnorm(-testcov)
total1=Networks.Fast(pvalue,net,iter=5,nburns=2,
v=20,algorithms="DPM",DPparallel=FALSE,
piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(1, 2, 5, 10, 15)
)

######summary the object of class "Networks.Fast"
summary(total1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
