#####Transferring large statistical pvalue into testing statistics
Transfer=function(pvalue)
{
  rstat=-(stats::qnorm(pvalue))
  return(rstat)
}
