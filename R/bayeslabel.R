#' @title {TransZ}
#' @description Compute a special transition matrix of Z: symmetric, down-regulated genes are not able to jump directly to up-regulated genes; a probability of a down-regulated gene jumping to null gene equals the probability of a up-regulated gene jumping to a null gene.
#' @param p The probability of a down-regulated gene jumping to null gene
#' @return A 3*3 transition matrix for Z
#' @keywords internal
#' @export
TransZ=function(p){
  m=matrix(c(p,(1-p)/2,0,(1-p),p,(1-p),0,(1-p)/2,p),nrow=3)
  return(m)
}


#' @title {biGaussianNull}
#' @description Compute a prior null distribution without given any prior biological knowledge. Methods used to estimate this prior include: "biGaussian","biGaussianWith0Cutoff","biGaussianMean0","biGaussianModeSplit" etc.
#' @param rstat The observed test statistics.
#' @param null.quantile Data quantiles, pre-fixed so that we could reduce the impact of extreme values for estimating prior null distribution. Default value is c(0.25, 0.75).
#' @param method A char. The methods we used to do estimation: either "biGaussian"-- EM algorithm for mixtures of two univariate normals, or "biGaussianWith0Cutoff"-- assume all negative test statistics forms one normal and all positive test statistics forms the other one normal. And proportion parameter is proportional to a number of observations each class, or "biGaussianMean0"-- null is formed by two half normals. "biGaussianModeSplit"-- split data from median value, then flip each part to the other side to estimate normal distribution.
#' @return A list with element
#'  \item{alpha}{cutoff value (if given)}
#'  \item{mu}{mean positon for two normals}
#'  \item{sd}{sd value for two normals}
#'  \item{prop}{proportion for each normal}
#' @keywords internal
#' @export
biGaussianNull=function(rstat,null.quantile=c(0.25, 0.75),method=c("biGaussianMean0","biGaussianModeSplit")){
  if(sum(is.na(rstat))>0){
    rstat=rstat[!is.na(rstat)]
  }
  rr=stats::quantile(rstat,prob=null.quantile)
  rr.center=rstat[which(rstat>rr[1] & rstat<rr[2])]
  methods=match.arg(method)
  if(methods=="biGaussianMean0"){
    rr.pos=rr.center[rr.center>0]
    rr.neg=rr.center[rr.center<=0]
    sd=sqrt(c(stats::var(rr.pos),stats::var(rr.neg)))
    res=list(alpha=0,mu=c(0,0), sd=sd, prop=sd/sum(sd))
  }else if(methods=="biGaussianModeSplit"){
    cutoff=stats::quantile(rr,0.5)
    rr.pos=rr.center[rr.center>cutoff]
    rr.pos=c((2*cutoff-rr.pos),rr.pos)
    rr.neg=rr.center[rr.center<=cutoff]
    rr.neg=c(rr.neg,(2*cutoff-rr.neg))
    res=list(alpha=cutoff,mu=c(cutoff,cutoff), sd=sqrt(c(stats::var(rr.pos),stats::var(rr.neg))), prop=c(length(rr.neg),length(rr.pos))/(length(rr.neg)+length(rr.pos)))
  }
  return(res)
}


#' @title {myKLdivergenece}
#' @description Compute Kullback Leibler divergence between biologic null distribution vs. proposed null distribution
#' @param nulldens A list, mu: mean location for each of normal density; sd: sd value for each of normal density; proportion: proportion for each normal density.
#' @param dens A list, mu: mean location for each of normal density; sd: sd value for each of normal density; proportion: proportion for each normal density.
#' @param integral A two element vector, the first element is integration lower bound; the second element is integration upper bound.
#' @param precision A number. The numerical precision for integration calculating. Default=0.001
#' @param method The way proposal density is calculated from wether it is "MixNorm" or "biGaussianMean0".
#' @return A number, KL distance between nulldens and dens
#' @keywords internal
#' @export
myKLdivergenece=function(nulldens,dens,integral=c(-6,6),precision=0.001,method=c("MixNorm","biGaussianMean0")){
  x=seq(integral[1],integral[2],by=precision)
  p=dens$prop/sum(dens$prop)
  y=sapply(1:length(p),function(pp){
    p[pp]*stats::dnorm(x,mean=dens$mu[pp],sd=dens$sd[pp])
  })
  y=rowSums(y)

  if(method=="biGaussianMean0"){
    ynull=x
    ll=which(x>=nulldens$alpha)
    ynull[ll]=stats::dnorm(x[ll],mean=nulldens$mu[2],sd=nulldens$sd[2])
    ll=which(x<nulldens$alpha)
    ynull[ll]=stats::dnorm(x[ll],mean=nulldens$mu[1],sd=nulldens$sd[1])
  }else{
    ynull=sapply(1:length(nulldens$prop),function(pp){
      nulldens$prop[pp]*stats::dnorm(x,mean=nulldens$mu[pp],sd=nulldens$sd[pp])
    })
    ynull=rowSums(ynull)
  }
  #   plot(x,ynull,type="l")
  KL=flexmix::KLdiv(cbind(ynull,y),eps =1e-9)
  return(KL[1,2])
}



#' @title {DPdensCluster}
#' @description First, apply DP-package to split data into several small normal groups based on their test statistics, then merge until we attain total three groups based on K-L distance between proposed null density vs. prior null density from biological background knowledge. Finally, return each group's mixture normal parameters as a list.
#' @param para The parameter list
#' @return The updated parameter list
#' @keywords internal
#' @export
DPdensCluster=function(para,twogrps=FALSE){
  #### initialized L-1 L0 L1 as sum(L_k)=fit$state$ncluster
  # set.seed(para$DPM$seed)

  fit=DPpackage::DPdensity(para$rStat[!is.na(para$rStat)],prior=para$DPM$prior,mcmc=para$DPM$mcmc,status=TRUE)
  utils::flush.console()

  #### sort sum(L_k)
  nclust=fit$state$ncluster
  prop=table(fit$state$ss)/sum(table(fit$state$ss))
  mu=fit$state$muclus[1:nclust]
  sigma=fit$state$sigmaclus[1:nclust]
  index=order(mu,decreasing=FALSE)
  mu=mu[index]
  sigma=sigma[index]
  prop=prop[index]

  KLdist=c()
  for(id in c(1:nclust)){
    # print(id)
    KLdist=c(KLdist,BANFF::myKLdivergenece(para$priorNull,list(mu=mu[id],sd=sigma[id],prop=prop[id]),integral=para$DPM$KLrange,precision=para$DPM$KLprecision,para$DPM$KLNullmethod))
  }
  nullclassID=which(KLdist==min(KLdist,na.rm=TRUE))

  if(max(nullclassID)==1){nullclassID=2}
  if(min(nullclassID)==nclust){nullclassID=nclust-1}
  KLD=KLdist[nullclassID][1]
  KLdiff=KLD

  while(nclust>3 & KLdiff>0){
    KLdist=c()
    # print(nclust)
    for(id in c(min(nullclassID)-1,min(max(nullclassID)+1,length(mu)))){
      # print(id)
      KLdist=c(KLdist,BANFF::myKLdivergenece(para$priorNull,list(mu=mu[c(nullclassID,id)],sd=sigma[c(nullclassID,id)],prop=prop[c(nullclassID,id)]),integral=para$DPM$KLrange,precision=para$DPM$KLprecision,para$DPM$KLNullmethod))
    }
    ll=which(KLdist==min(KLdist,na.rm=TRUE))
    KLdiff=KLD-KLdist[ll]
    if(KLdiff>0){
      nullclassID=c(nullclassID,c(min(nullclassID)-1,max(nullclassID)+1)[ll])
      nullclassID=unique(nullclassID)
      KLD=KLdist[ll]
      nclust=nclust-1
    }
  }

  if(twogrps){
    idlist=list(c(1:(min(nullclassID)-1)),sort(c(nullclassID,(max(nullclassID)+1):length(mu))))
  }else{
    idlist=list(c(1:(min(nullclassID)-1)),sort(nullclassID),c((max(nullclassID)+1):length(mu)))
  }

  if(max(table(unlist(idlist)))==1){
    ### empty para$mulist:
    para$mulist=lapply(1:para$numZ,function(k){list()})
    para$sdlist=lapply(1:para$numZ,function(k){list()})
    para$proplist=lapply(1:para$numZ,function(k){list()})

    para$zLab=rep(2,length(para$rStat))
    for(k in 1:length(idlist)){
      para$mulist[[k]]=mu[idlist[[k]]]
      para$sdlist[[k]]=sigma[idlist[[k]]]
      para$proplist[[k]]=prop[idlist[[k]]]/sum(prop[idlist[[k]]])
      ll=index[idlist[[k]]]
      para$zLab[which(is.element(fit$state$ss,ll))]=k
    }
  }
  return(para)
}

#' @title {z.update.Mclustprior}
#' @description z.update.Mclustprior
#' @param para The parameter list
#' @param oldz NULL
#' @param rho NULL
#' @param pi NULL
#' @return voted z
#' @keywords internal
#' @export
z.update.Mclustprior=function(para,oldz,rho,pi){
  zTrack=c()
  for(iter in 1:para$hyperPara$niter.z){
    newz=vector(length=length(oldz))
    for(node in 1:para$numNodes){
      nbr=which(para$net[node,]>0)
      loglike=sapply(1:length(para$zSet),function(z){
        log(pi[z])+rho[z]*sum(oldz[nbr]==z)+log(stats::dnorm(para$rStat[node],mean=para$hyperPara$mean.z[z],sd=para$hyperPara$sd.z[z]))
      })
      if(sum(is.infinite(loglike))>0){
        loglike[which(is.infinite(loglike))]=para$hyperPara$replaceInf
      }
      sampleprob=sapply(1:length(loglike),function(z){
        vv=loglike-loglike[z]
        return(1/sum(exp(vv)))
      })
      newz[node]=sample(para$zSet,1,prob=sampleprob)
      # print(node)
    }
    zTrack=cbind(zTrack,newz)
    oldz=newz
  }

  zvote=apply(zTrack,1,function(x){
    a=length(which(x==1))
    aa=length(which(x==2))
    aaa=length(which(x==3))
    mm=max(a,aa,aaa)
    return(ifelse(mm==a, 1, ifelse(mm==aa,2,3)))
  })
  return(zvote)
}

#' @title {llikDMH}
#' @description llikDMH
#' @param para NULL
#' @param z NULL
#' @param rho NULL
#' @param pi NULL
#' @return llik Value
#' @keywords internal
#' @export
llikDMH=function(para,z,rho,pi){
  res=sapply(1:para$numNodes,function(node){
    nbr=which(para$net[node,]>0)
    return(log(pi[z[node]])+rho[z[node]]*sum(z[nbr]==z[node])+stats::dnorm(para$rStat[node],mean=para$hyperPara$mean.z[z[node]],sd=para$hyperPara$sd.z[z[node]]))
  })
  return(sum(res))
}


#' @title {hyperParaDMH}
#' @description Apply Double Metropolis-Hastings sampler (DMH) to find hyperparameters: rho & pi
#' @param para The parameter list
#' @param plot Logical. T: plot trace plot; F: not plot output.
#' @return a list:
#' \item{rhoTrack}{saved rho value each iteration}
#' \item{piTrack}{saved pi value each iteration}
#' \item{zTrack}{saved z label value each iteration}
#' \item{rho}{posterior mean calculated based on second half of the iterations}
#' \item{zstat}{posterior z value calculated based on major vote method.}
#' @keywords internal
#' @export
hyperParaDMH=function(para,plot=TRUE){
  rho=para$hyperPara$rhostat
  pi=para$hyperPara$pistat
  z=para$zLab
  rhoTrack=rho
  piTrack=pi
  zTrack=z
  rho.L=length(rho)
  pi.L=length(pi)
  for(iter in 1:para$hyperPara$niter){
    #### revise initial values for rho & pi
    repeat{
      newrho=truncnorm::rtruncnorm(rho.L,mean=rho,sd=para$hyperPara$rhosd,a=para$hyperPara$rhoLowB,b=para$hyperPara$rhoUpB)
      newrho[rho.L]=0
      newpi=truncnorm::rtruncnorm(pi.L,mean=pi,sd=para$hyperPara$pisd,a=para$hyperPara$piLowB,b=para$hyperPara$piUpB)
      newpi[2]=1-sum(newpi[-2])
      if(min(newrho[-c(2,rho.L)])>newrho[2] & newpi[2]>0.5) break;
    }
    ### updating z
    newz=BANFF::z.update.Mclustprior(para,z,newrho,newpi)
    AccRate=BANFF::llikDMH(para,newz,rho,pi)+BANFF::llikDMH(para,z,newrho,newpi)-BANFF::llikDMH(para,z,rho,pi)-BANFF::llikDMH(para,newz,newrho,newpi)
    if(is.nan(AccRate)){
      accrate=1
    }else{
      accrate=min(AccRate,0)
    }
    if(log(stats::runif(1))<accrate){
      ### accept this new state.
      rho=newrho
      pi=newpi
      z=newz
    }
    rhoTrack=cbind(rhoTrack,rho)
    piTrack=cbind(piTrack,pi)
    zTrack=cbind(zTrack,z)
  }
  if(plot){
    graphics::par(mfrow=c(2,2))
    for(rr in 1:nrow(rhoTrack)){
      plot(seq(1:ncol(rhoTrack)),rhoTrack[rr,],type="l",main=paste0("DMH for rho",rr),xlab="iterations",ylab=paste0("rho",rr))
    }
    #     dev.off()
    graphics::par(mfrow=c(1,3))
    for(rr in 1:nrow(piTrack)){
      plot(seq(1:ncol(piTrack)),piTrack[rr,],type="l",main=paste0("DMH for pi",rr),xlab="iterations",ylab=paste0("pi",rr))
    }
    #     dev.off()
  }
  pi=vector(length=pi.L)
  pi[-2]=rowMeans(matrix(piTrack[-2,round(0.5*para$hyperPara$niter):para$hyperPara$niter],nrow=(pi.L-1)))
  pi[2]=1-sum(pi[-2])
  rho=rowMeans(rhoTrack[,round(0.5*para$hyperPara$niter):para$hyperPara$niter])
  zsw=apply(zTrack,1,function(x){
    l=length(which(x==1))
    ll=length(which(x==2))
    lll=length(which(x==3))
    mm=which(ll==max(ll))
    return(mm)
  })
  res=list(rhoTrack=rhoTrack,piTrack=piTrack,zTrack=zTrack,pi=pi,rho=rho,zstat=zsw)
  return(res)
}


#' @title {GCut.z}
#' @description Split network into three where nodes located within each group has the same class label (z value).
#' @param para The parameter list
#' @return The updated parameter list
#' @keywords internal
#' @export
GCut.z=function(para){
  GCut=list()
  for(z in para$zSet){
    ll=which(para$zLab==z)
    G=igraph::graph.adjacency(para$net[ll,ll],mode="undirected")
    G=igraph::set.vertex.attribute(G,"node",value=ll)
    GCut[[z]]=G
  }
  para$GCut=GCut
  return(para)
}


#' @title {SWCut}
#' @description Swendsen-Wang graph cut techniques. The edges located within each subgraph (same z) has a positive probability to stay "on" or "off". The edges with "off" labels are closed so that it is further cut to smaller subgraphs.
#' @param para parameter list
#' @param rhoScale a 4-element vector, used when rho values need to be scaled for a propor graph cut. Default=c(1,1,1,0)
#' @return  A list with three elements, representing z=1, z=2, or z=3. Each element is also a list, contains all subgraph with same z value.
#' @keywords internal
#' @export
SWCut=function(para){
  ProbEOn=1-exp(-para$rho[-length(para$rho)])
  gList=list()
  gl=1

  locvec=c()

  for(k in 1:length(para$GCut)){
    subg=para$GCut[[k]]
    #     plot(subg)
    E=igraph::E(subg)
    V=igraph::V(subg)
    vnode=igraph::get.vertex.attribute(subg,"node")
    if(length(vnode)==0){
      cat("Empty subgraph!","\n")
    }
    myz=para$zLab[vnode]
    if(length(E)>0){
      Eid=stats::rbinom(length(E),1,ProbEOn[myz[1]])
      newsubg=igraph::delete.edges(subg,(E[which(Eid==0)]))
      newsubg=igraph::clusters(newsubg)
      gList[gl:(gl+newsubg$no-1)]=lapply(1:newsubg$no,function(kk){
        igraph::induced.subgraph(graph=subg,vids=V[which(newsubg$membership==kk)])
      })

      gl=gl+newsubg$no
    }else{
      nnode=length(V)
      gList[gl:(gl+nnode-1)]=lapply(1:nnode,function(kk){
        mygraph=igraph::graph.empty(0, directed=FALSE)
        mygraph=igraph::add.vertices(mygraph,nv=1)
        mygraph=igraph::set.vertex.attribute(mygraph,"node",value=vnode[kk])
        return(mygraph)
      })
      gl=gl+nnode
    }
    locvec=c(locvec,vnode)
  }
  #   e=sapply(gList,function(x){length(E(x))})
  return(gList)
}


#' @title {logdensMixNorm}
#' @description logdensMixNorm
#' @param dat NULL
#' @param mu NULL
#' @param sd NULL
#' @param prop NULL
#' @return log density value
#' @keywords internal
#' @export
logdensMixNorm=function(dat,mu,sd,prop){
  res=sapply(1:length(mu),function(l){
    prop[l]*stats::dnorm(dat,mean=mu[l],sd=sd[l])
  })
  res=as.matrix(res)
  if(ncol(res)!=length(mu)){
    res=t(res)
  }
  x=rowSums(res)
  loc=which(x==0)
  if(length(loc)>0){
    #     cat("inf log density:,",loc,"\n")
    #     x=x[-loc]
    return(-999999)
  }else{
    return(sum(log(x)))
  }
}


#' @title {z.update.fastSW}
#' @description Given density specification. Update regulation type (z value) within each subgraph cut by Swendsen-Wang graph cut techniques.
#' @param para The parameter list
#' @param swcutList a list, the direct output from SWCut.
#' @return An updated parameter list.
#' @keywords internal
#' @export
z.update.fastSW=function(para,swcutList){
  for(gph in 1:length(swcutList)){
    if((min(table(para$zLab))<para$min.node) | (length(unique(para$zLab))<para$numZ)){
      next;
    }
    # cat(gph,"\n")
    currt.gph=swcutList[[gph]]
    node.gph=igraph::get.vertex.attribute(currt.gph,"node")
    # plot(currt.gph)
    nnode=length(node.gph)
    oldz=para$zLab[node.gph][1]
    nbr=lapply(1:length(node.gph),function(kk){which(para$net[node.gph[kk],]>0)})
    prob.z=sapply(1:length(para$zSet),function(z){
      mk=sapply(nbr,function(nbrs){sum(para$zLab[nbrs]==z)})
      return(log(para$pivec[z])*length(node.gph)+(para$rho[z]*sum(mk))+BANFF::logdensMixNorm(para$rStat[node.gph],mu=para$mulist[[z]],sd=para$sdlist[[z]],prop=para$proplist[[z]]))
    })
    if(sum(is.infinite(prob.z))>0){
      loc=which(is.infinite(prob.z))
      prob.z[loc]=1
    }
    sampleprob=sapply(1:length(prob.z),function(z){
      vv=prob.z-prob.z[z]
      return(1/sum(exp(vv)))
    })
    newz=sample(size=1,x=seq(1,length(prob.z)),prob=sampleprob)

    AccRate=exp(prob.z[newz]-prob.z[oldz])*para$TransZ[newz,oldz]/para$TransZ[oldz,newz]
    if(is.nan(AccRate)){
      accrate=0
    }else{
      accrate=min(AccRate,1)
    }

    if(stats::runif(1)<accrate){
      #       acclist[gph]=TRUE
      ### accept this new state.
      if((newz!=oldz)){
        if((para$mk[oldz]-nnode)>para$min.node){
          para$zLab[node.gph]=newz
          para$mk[c(oldz,newz)]=para$mk[c(oldz,newz)]+c(-1,+1)*nnode
        }
      }
    }
    #   return(list(para=para,ztrack=zz,problist=problist,acclist=acclist))
  }
  return(para)
}


#' @title {DPdensitySubset}
#' @description DPdensitySubset
#' @param ll NULL
#' @param subdat NULL
#' @param subprior NULL
#' @param submcmc NULL
#' @param substatus NULL
#' @return para
#' @keywords internal
#' @export
DPdensitySubset=function(ll,subdat,subprior,submcmc,substatus){
  fit=DPpackage::DPdensity(subdat[ll],prior=subprior,mcmc=submcmc,status=substatus)
  utils::flush.console()
  nclust=fit$state$ncluster
  prop=table(fit$state$ss)/sum(table(fit$state$ss))
  mu=fit$state$muclus[1:nclust]
  sigma=fit$state$sigmaclus[1:nclust]
  index=order(mu,decreasing=FALSE)
  oldg=fit$state$ss
  newg=oldg
  k=1
  gloc=list()
  for(g in index){
    lll=which(oldg==g)
    newg[lll]=k
    gloc[[k]]=ll[lll]
    k=k+1
  }
  return(list(mu=mu[index],sigma=sqrt(sigma)[index],prop=prop[index],gvec=newg,gloc=gloc))
}


#' @title {DensityDPOne}
#' @description Given regulation type (z value). Update density specification based on DPM fitting.
#' @param para parameter list.
#' @return An updated parameter list.
#' @keywords internal
#' @export
DensityDPOne=function(para){
  for(z in para$zSet){
    ll=which(para$zLab==z)
    if(length(ll)>para$min.node){
      myprior=para$HODCPara$prior[[z]]
      # myprior$m1=para$mulist[[z]]
      # nn=length(para$sdlist[[z]])
      # myprior$psiinv1=diag(para$sdlist[[z]],nrow=nn,ncol=nn)
      res=BANFF::DPdensitySubset(ll,subdat=para$rStat,subprior=myprior,submcmc=para$HODCPara$mcmc,substatus=TRUE)
      # utils::flush.console()
      para$mulist[[z]]=res$mu
      para$sdlist[[z]]=sqrt(res$sigma)
      para$proplist[[z]]=res$prop
      para$gvec[ll]=res$gvec
      para$zgloc[[z]]=res$gloc
    }else{
      para$mulist[[z]]=mean(para$rStat[ll])
      para$sdlist[[z]]=sqrt(stats::var(para$rStat[ll]))
      para$proplist[[z]]=1
      para$gvec[ll]=rep(1,length(ll))
      para$zgloc[[z]]=list(ll)
    }
  }
  return(para)
}


#' @title {r.knnImpute}
#' @description Impute missing node using their nearest neighbor test statistics.
#' @param para parameter list.
#' @return  An updated parameter list.
#' @keywords internal
#' @export
r.knnImpute=function(para){
  res=sapply(1:length(para$misLoc),function(k){
    nbr=para$nbrMis[[k]]
    rr=para$rStat[nbr]
    return(mean(rr))
  })
  para$rStat[para$misLoc]=res
  return(para)
}


#' @title {r.BayesImpute}
#' @description Impute missing node by sampling from posterior density.
#' @param para parameter list.
#' @return An updated parameter list.
#' @keywords internal
#' @export
r.BayesImpute=function(para){
  zloc=lapply(para$zgloc,function(x){unlist(x)})
  miszloc=sapply(para$misLoc,function(k){
    a=sapply(1:length(zloc),function(z){
      is.element(k,zloc[[z]])
    })
    which(a)
  })
  nn=length(para$misLoc)
  miszgloc=sapply(1:nn,function(k){
    tt=sapply(para$zgloc[[miszloc[k]]],function(x){is.element(para$misLoc[k],x)})
    return(which(tt))
  })
  miszgloc=cbind(miszloc,miszgloc)
  #   print(mizgloc[1,])
  mm=sapply(1:nn,function(k){c(para$mulist[[miszgloc[k,1]]][[miszgloc[k,2]]],para$sdlist[[miszgloc[k,1]]][[miszgloc[k,2]]])})
  rr=stats::rnorm(nn,mean=mm[1,],sd=mm[2,])
  #   print(rr[1])
  para$rStat[para$misLoc]=rr
  #   print(para$rStat[773])
  return(para)
}


#' @title {Bayesian nonparametric feature selection over large-scale networks with missing values}
#' @aliases BANFF2
#' @name BANFF2
#' @usage
#' BANFF2(net,test.stat,pvalue.stat=FALSE,candidate.z.set=c(-1,0,1),
#' seed.main=1024,na.action=c("NN","Bayes","na.remove"),niter.densupd=5,niter=10,
#' paras=list(tau=c(2,10,2),alpha=NULL,gamma=NULL,xi=NULL, beta=rep(10,3),
#' rho=c(1.003,0.479,0.988,0.000),pivec=c(0.15,0.7,0.15),densAcc=0.001,
#' null.quantile=c(0.25, 0.75),null.method="biGaussianModeSplit",
#' transitionMatrix.Z.11=0.6,miss.stat=2,min.node=5),
#' para.DPM=NULL,para.HODC=NULL,para.DMH=NULL)
#' @description Main function. Two steps: Given density specification, update selection indicator z by Swendsen- Wang; Given selection indicator z, update density specification by DPM fitting.
#' @param net The adjacent matrix with 0/1 indicating "connected" or "not directly connected
#' @param test.stat The observed test statistics. Missing values are represented as NAs. If they are pvalues, then the pvalue.stat should be T;
#' @param pvalue.stat Logical. Wether test.stat is generated as pvalues or not. Default F.
#' @param candidate.z.set Default is of three regulation type. Defalut=c(-1,0,1), 1=down-regulated, 2=not differentially expressed, 3=up-regulated.
#' @param seed.main Set seed before iteration for generating reproducible results. Default=1024.
#' @param na.action The method used to impute missing values. Can be "NN", "Bayes", or "na.remove".
#' @param niter.densupd The total number of iterations for updating density. Default=5
#' @param niter The total number of iterations for study. Default=10.
#' @param paras A list contains hyper-parameters and other parameters used for preparations.
#' \itemize{
#' \item niter.densupd The iteration is from 1 to the maximum steps when we update density specification by DPM. Default=20.
#' \item tau A three-element vector, default=c(2,10,2);
#' \item alpha A three-element vector. Default=NULL.
#' \item gamma A three-element vector. Default=NULL.
#' \item xi A three-element vector. Default=NULL.
#' \item beta A three-element vector. Default=rep(10,3).
#' \item rho A four-element vector. Default=c(1.003,0.479,0.988,0.000), indicating local smoothness for Potts prior. Note: the default value is calculated based on data(net) strucutre by DMH.
#' \item pivec A three-element vector. Default=c(0.15,0.7,0.15). Contains prior knowledge globally about selection indicator z.
#' \item densityAcc A number, need to specify precision for K-L integration when to use the numerical approximation. Default=0.001.
#' \item null.quantile A two element vector representing lower quantile and upper quantile for calculating prior null density if not given by biologists. Default=c(0.25, 0.75).
#' \item null.method A char. The method we used to estimate null density: "biGaussian"-- EM algorithm for mixtures of two univariate normals; "biGaussianWith0Cutoff"-- assume all negative test statistics forms one normal and all positive test statistics forms the other one normal. And proportion parameter is proportional to a number of observations each class; "biGaussianMean0"-- null is formed by two half normals. "biGaussianModeSplit"-- split data from median value, then flip each part to the other side to estimate normal distribution.
#' \item transitionMatrix.Z.11 [1,1] element in transition matrix for z. Default=0.6.
#' \item miss.stat impute NAs in test.test when apply Double Metropolis-Hastings sampler (DMH) to find hyperparameters: rho & pi.
#' \item min.node The minimum number of nodes in each group.
#' }
#' @param para.DPM A list object contains, if NULL, default value is used:
#' \itemize{
#' \item niter default=10
#' \item nsample default=10
#' \item KLrange default=c(-6,6), usually we consider wider range than c(floor(min(test.stat,na.rm=TRUE)),ceiling(max(test.stat,na.rm=TRUE)))
#' \item KLprecision default=0.001
#' \item KLNullmethod default="biGaussianMean0"
#' \item mcmc a list, default=list(nburn=10000,nsave=100,nskip=0,ndisplay=10000)
#' \item prior a list, default=list(alpha=3,m1=rep(0,1),psiinv1=diag(0.5,1),nu1=4,tau1=1,tau2=100)
#'}
#' @param para.HODC A list object contains, if NULL, default value is used:
#' \itemize{
#' \item nsample default=10
#' \item KLrange default=c(-6,6), usually we consider wider range than c(floor(min(test.stat,na.rm=TRUE)),ceiling(max(test.stat,na.rm=TRUE)))
#' \item KLprecision default=0.001
#' \item KLNullmethod default="biGaussianMean0",
#' \item mcmc a list, default=list(nburn=1000,nsave=100,nskip=0,ndisplay=1000)
#' \item prior a list, defaut is a list object where each of the element specify the prior used when fitting each density for class labels z. For each of the class, default parameters are the same, a list contains: alpha=3,m2=rep(0,1),s2=diag(100000,1),psiinv2=diag(temp.sdlist[1],1),nu1=4,nu2=4,tau1=1,tau2=100
#'}
#' @param para.DMH If rho & pivec is not given, DMH is used for pre-calculating rho & pivec. Default is a list object contains:
#' \itemize{
#' \item niter default=1000
#' \item pistat default=c(0.25,0.5,0.25)
#' \item pisd default=rep(0.03,3)
#' \item rhostat default=c(1,0.5,1,0)
#' \item rhosd default=rep(0.03,4)
#' \item rhoLowB default=c(0,0,0,0)
#' \item rhoUpB default=c(1.5,1.5,1.5,1.5)
#' \item piLowB default=c(0,0,0)
#' \item piUpB default=c(1,1,1)
#' \item niter.z default=1
#' \item replaceInf default=-99999
#' \item DMHplot default=FALSE
#' }
#'
#' @return A list:
#' \item{initialValue}{initial parameter list}
#' \item{zTrack}{trace for z}
#' \item{FinalValue}{final parameter list}
#' \item{iters}{total iterations}
#' \item{rmisTrack}{(if NAs in test.statistics) trace for test.statistics imputation. (only for those with NAs)}
#' @details The fully Bayesian updating algorithm is executed as below:
#' \itemize{
#' \item Input data r and graph G=<V,E>
#' \item Update z|theta via Swendsen-Wang
#' \item Update theta|z via DPM Fitting
#' }
#' @examples
#' \dontrun{
#' ## The simulation settings based on real gene network (takes time)
#' data("net")
#' data("test.stat")
#' res=BANFF2(net,test.stat,niter=300,na.action="NN")
#' res=BANFF2(net,pnorm(test.stat),pvalue.stat=TRUE,candidate.z.set=c(0,1),na.action="NN",
#' niter=300,
#' paras=list(tau=c(2,10),alpha=NULL,gamma=NULL,xi=NULL, beta=rep(10,2),rho=c(1,0.5,0),
#' pivec=c(0.2,0.8),densAcc=0.001,null.quantile=c(0.25, 1),
#' null.method="biGaussianModeSplit",transitionMatrix.Z.11=0.6,miss.stat=2,min.node=5))
#'
#' ## A toy example
#' simdata=SimulatedDataGenerator(nnode=100,missing=TRUE,missrate=0.1,dist="norm",
#' plot=TRUE,nbin=c(20,20,10),rng=1024)
#' res=BANFF2(net=simdata$net,test.stat=simdata$testcov,niter=100,na.action="NN")
#' classLabelEst=SummaryClassLabel(simdata$net,simdata$testcov,res$zTrack,
#' method="MajorVote",nburn=10)
#' print(table(classLabelEst))
#' }
#' @export
BANFF2=function(net,test.stat,pvalue.stat=FALSE,candidate.z.set=c(-1,0,1),seed.main=1024,na.action=c("NN","Bayes","na.remove"),niter.densupd=5,niter=10,paras=list(tau=c(2,10,2),alpha=NULL,gamma=NULL,xi=NULL, beta=rep(10,3),rho=c(1.003,0.479,0.988,0.000),pivec=c(0.15,0.70,0.15),densAcc=0.001,null.quantile=c(0.25, 0.75),null.method="biGaussianModeSplit",transitionMatrix.Z.11=0.6,miss.stat=2,min.node=5),para.DPM=NULL,para.HODC=NULL,para.DMH=NULL){

  if(!requireNamespace("mclust", quietly = TRUE)) {
    stop("Please load package: mclust as required ")
  }
  if(!requireNamespace("DPpackage", quietly = TRUE)) {
    stop("Please load package: DPpackage as required ")
  }

  set.seed(seed.main)
  # library(DPpackage,verbose = F,quietly = T)
  # library(mclust,verbose = F,quietly = T)
  if(pvalue.stat){
    test.stat=-stats::qnorm(test.stat)
  }

  #### if NA.action=NA.remove:
  if(sum(is.na(test.stat))>0){
    method.na.impute=match.arg(na.action)
    if(method.na.impute=="na.remove"){
      na.loc=which(is.na(test.stat))
      net=net[-na.loc,-na.loc]
      test.stat=test.stat[-na.loc]
    }
  }else{
    method.na.impute="NoImpute"
  }
  #### create a parameter list:
  numz=length(candidate.z.set)
  paraList=list(zSet=seq(1,numz),numZ=numz,pivec=paras$pivec,rho=paras$rho, rStat=test.stat, misLoc=which(is.na(test.stat)),numNodes=length(test.stat), zLab=NULL, tau=paras$tau, beta=paras$beta, alpha=rep(3,numz), gamma=NULL, xi=NULL, net=net, densityAccu=paras$densAcc,null.quantile=paras$null.quantile,null.method=paras$null.method,TransZ11=paras$transitionMatrix.Z.11,densIter=list(to=niter.densupd,niter=niter),min.node=paras$min.node)
  if(length(paraList$misLoc)>0){
    Imputeflag=TRUE
  }else{
    Imputeflag=FALSE
  }
  paraList$priorNull=biGaussianNull(paraList$rStat,paraList$null.quantile,method=paraList$null.method)
  # temp=mixtools::normalmixEM(paraList$rStat[!is.na(paraList$rStat)],k=numz)
  # paraList$mulist=as.list(temp$mu)
  # paraList$sdlist=as.list(temp$sigma)
  temp=mclust::Mclust(paraList$rStat[!is.na(paraList$rStat)],G=3,modelName="E")
  paraList$mulist=as.list(temp$parameter$mean)
  temp.sdlist=temp$parameter$variance$sigmasq
  if(length(temp.sdlist)!=3){
    temp.sdlist=rep(sqrt(temp.sdlist),3)
  }
  paraList$sdlist=as.list(temp.sdlist)
  # paraList$zLab=temp$classification
  paraList$proplist=list(1,1,1)

  if(is.null(paras$alpha)){
    paraList$alpha=paraList$beta/unlist(paraList$sdlist)+1
  }else{
    paraList$alpha=paras$alpha
  }
  if(is.null(paras$gamma)){
    paraList$gamma=unlist(paraList$mulist)
  }else{
    paraList$gamma=paras$gamma
  }
  if(is.null(paras$xi)){
    paraList$xi=unlist(paraList$sdlist)*sqrt(2)
  }else{
    paraList$xi=paras$xi
  }
  paraList$TransZ=TransZ(paraList$TransZ11)

  # mcmc.default=list(nburn=1000,nsave=100,nskip=0,ndisplay=1000)
  mcmc.default=list(nburn=10000,nsave=1000,nskip=0,ndisplay=10000)
  # prior.default=list(alpha=3,m2=rep(0,1),s2=diag(100000,1),psiinv2=diag(100000,1),nu1=4,nu2=4,tau1=1,tau2=100)
  prior.default=list(alpha=3,m1=0,psiinv1=0.5,nu1=4,tau1=1,tau2=100)
  para.DPM.default=list(niter=10,nsample=10,KLrange=c(floor(min(test.stat,na.rm=TRUE)),ceiling(max(test.stat,na.rm=TRUE))),KLprecision=0.001,KLNullmethod="biGaussianMean0")
  vars=setdiff(names(mcmc.default),names(para.DPM$mcmc))
  para.DPM$mcmc[vars]=lapply(vars,function(x){mcmc.default[[x]]})
  vars=setdiff(names(prior.default),names(para.DPM$prior))
  para.DPM$prior[vars]=lapply(vars,function(x){prior.default[[x]]})
  vars=setdiff(names(para.DPM.default),names(para.DPM))
  para.DPM[vars]=lapply(vars,function(x){para.DPM.default[[x]]})
  paraList$DPM=para.DPM


  paraList=DPdensCluster(paraList,twogrps=FALSE)
  utils::flush.console()
  paraList$zLab[!is.na(paraList$rStat)]=temp$classification
  paraList$zLab[is.na(paraList$rStat)]=rep(2,length(paraList$misLoc))
  paraList$mk=table(paraList$zLab)[paraList$zSet]

  mcmc.default=list(nburn=1000,nsave=100,nskip=0,ndisplay=1000)
  vars=setdiff(names(mcmc.default),names(para.HODC$mcmc))
  para.HODC$mcmc[vars]=lapply(vars,function(x){mcmc.default[[x]]})
  prior.default=list(alpha=3,m2=rep(0,1),s2=diag(100000,1),psiinv2=diag(100000,1),nu1=4,nu2=4,tau1=1,tau2=100)
  para.HODC$prior=lapply(1:numz,function(x){
    vars=setdiff(names(prior.default),names(para.HODC$prior[[x]]))
    para.HODC$prior[[x]][vars]=lapply(vars,function(x){prior.default[[x]]})
    return(para.HODC$prior[[x]])})
  vars=setdiff(names(para.DPM.default),names(para.HODC))
  para.HODC[vars]=lapply(vars,function(x){para.DPM.default[[x]]})
  paraList$HODCPara=para.HODC

  # paraList$zgloc=lapply(1:numz,function(x){list()})


  if(is.null(paraList$rho) | is.null(paraList$pivec)){
    para.DMH.default=list(niter=1000,pistat=c(0.25,0.5,0.25),pisd=rep(0.03,3),rhostat=c(1,0.5,1,0),rhosd=rep(0.03,4),rhoLowB=c(0,0,0,0),rhoUpB=c(1.5,1.5,1.5,1.5),piLowB=c(0,0,0),piUpB=c(1,1,1),niter.z=1,replaceInf=-99999,DMHplot=FALSE,mean.z=temp$parameters$mean,sd.z=temp.sdlist)

    vars=setdiff(names(para.DMH.default),names(para.DMH))
    para.DMH[vars]=lapply(vars,function(x){para.DMH.default[[x]]})
    mypara=paraList
    mypara$hyperPara=para.DMH
    mypara$rStat[is.na(paraList$rStat)]=rep(mean(mypara$rStat,na.rm=TRUE),sum(is.na(mypara$rStat)))
    hyperRes=hyperParaDMH(mypara,plot=mypara$hyperPara$DMHplot)
    paraList$rho=hyperRes$rho
    paraList$pivec=hyperRes$pi
  }

  paraList$rho=c(paraList$rho[2],paraList$rho)
  paraList$pivec=c(paraList$pivec,paraList$pivec[1])
  print("Completed creating the parameter list!")

  #### Main Loop for Fast S-W Updates
  simRes=list()
  simRes$initialValue=paraList
  zTrack=paraList$zLab
  iter=0
  if(Imputeflag){
    paraList$rStat[paraList$misLoc]=0
    rmisTrack=paraList$rStat[paraList$misLoc]
    if(method.na.impute=="NN"){
      paraList$nbrMis=lapply(paraList$misLoc,function(k){
        which(paraList$net[k,]>0)
      })
      names(paraList$nbrMis)=paraList$misLoc
    }
  }

  ptm=proc.time()
  #### Main Loop + density updates:
  kk=1
  while(iter<paraList$densIter$to){
    ### updating z
    paraList=GCut.z(paraList)
    currtSWCut=SWCut(paraList)
    paraList=z.update.fastSW(paraList,currtSWCut)
    ### updating mu/sd/etc
    paraList=DensityDPOne(paraList)
    # utils::flush.console()
    if(Imputeflag & method.na.impute=="NN"){
      paraList=r.knnImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }else if(Imputeflag & method.na.impute=="Bayes"){
      paraList=r.BayesImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }
    ## keep track of z changes:
    zTrack=cbind(zTrack,paraList$zLab)
    kk=kk+1
    iter=iter+1
  }
  print("Stop updating density, only updating the class indicators")
  #### Main Loop - density updates:
  while(iter<paraList$densIter$niter){
    ### updating z
    paraList=GCut.z(paraList)
    currtSWCut=SWCut(paraList)
    paraList=z.update.fastSW(paraList,currtSWCut)
    if(Imputeflag & method.na.impute=="NN"){
      paraList=r.knnImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }else if(Imputeflag & method.na.impute=="Bayes"){
      paraList=r.BayesImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }
    ## keep track of z changes:
    zTrack=cbind(zTrack,paraList$zLab)
    iter=iter+1
    # print(iter)
  }
  print("Algorithm finished!")
  #### Results
  colnames(zTrack)=NULL
  simRes$zTrack=zTrack
  simRes$FinalValue=paraList
  simRes$iters=iter
  if(Imputeflag){
    simRes$rmisTrack=rmisTrack
  }
  return(simRes)
}
