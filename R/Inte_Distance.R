Inte_Distance=function(i,mclust)
{
  if (mclust$parameter$variance$modelName=="E"){mclust$parameter$variance$sigmasq=rep(mclust$parameter$variance$sigmasq,length(mclust$parameter$mean))}
  distance=stats::dnorm(mclust$parameter$mean[i],mclust$parameter$mean[i],sd=sqrt(mclust$parameter$variance$sigmasq[i]+mclust$parameter$variance$sigmasq[i]))-stats::dnorm(mclust$parameter$mean[i],mclust$parameter$mean[i+1],sd=sqrt(mclust$parameter$variance$sigmasq[i]+mclust$parameter$variance$sigmasq[i+1]))  -stats::dnorm(mclust$parameter$mean[i+1],mclust$parameter$mean[i],sd=sqrt(mclust$parameter$variance$sigmasq[i+1]+mclust$parameter$variance$sigmasq[i]))  +stats::dnorm(mclust$parameter$mean[i+1],mclust$parameter$mean[i+1],sd=sqrt(mclust$parameter$variance$sigmasq[i+1]+mclust$parameter$variance$sigmasq[i+1]))
  return(distance)
}
