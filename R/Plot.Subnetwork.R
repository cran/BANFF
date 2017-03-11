######Sub-network Plot
Plot.Subnetwork=function(net,trace)
{
  ztall=sapply(1:ncol(trace),function(kk) return(mean(trace[,kk])))
  eids=which(ztall>0.5)
  g <- igraph::graph.adjacency(net,mode="undirected")
  
  graphics::plot(g,layout=igraph::layout.kamada.kawai,mark.groups=eids,vertex.size=1,vertex.label=NA)
  
}
