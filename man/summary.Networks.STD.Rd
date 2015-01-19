\name{summary.Networks.STD}
\alias{summary.Networks.STD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing Bayesian Network Discovery using a Standard MCMC Algorithm

%%  ~~function to do ... ~~
}
\description{Summary method for class \code{"Networks.STD"}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{summary}{Networks.STD}(object, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"Networks.STD"} resulting of a call to \code{\link{Networks.STD}}.
 
%%     ~~Describe \code{x} here~~
}
 \item{\dots}{other arguments}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Networks.STD}}.}
\examples{
library(igraph)
###Creating the network of 10X10 image
g <- graph.lattice(length=10,dim=2)
net=as(get.adjacency(g,attr=NULL),"matrix")##this is the input of argument 'net'
##Assign the signal elements with signal intention 
##as normal distribution N(1,0.2). While noise is set as N(0,0.2)
newz=rep(0,100)
for (i in 3:7)
{
 newz[(i*10+3):(i*10+7)]=1
}
testcov<-0
for(i in 1:100){
 if(newz[i]==0){
   testcov[i]<-rnorm(1,mean=0,sd=0.2)

 }else{
  testcov[i]<-rnorm(1,mean=1,sd=0.2)

 }
}
##The profile of the image
image(matrix(testcov,10,10),col=gray(seq(0,1,length=255)))
##Transform the signals into pvalue form and begin identification
pvalue=pnorm(-testcov)
total=Networks.STD(pvalue,net,iter=3,nburns=1,
piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(0.5,1,5,10,15))

######summary the object of class "Networks.STD"
summary(total)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
