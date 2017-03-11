#' @title {Local FDR method for summarizing class labels}
#' @aliases SummaryClassLabel
#' @name SummaryClassLabel
#' @description The finalized class label can be either by method "Major vote", or "under local fdr control".
#' @usage SummaryClassLabel(net,test.stat,zTrack,method=c("MajorVote","LocfdrControl"),
#' plot=FALSE,nburn=50,nskip=2,locfdr.alpha=0.2)
#' @param net The adjacent matrix with 0/1 indicating "connected" or "not directly connected
#' @param test.stat The observed test statistics. Missing values are represented as NAs.
#' @param zTrack The trace for z. It is a num_of_genes by num_of_iterations matrix.
#' @param method A char. Either "MajorVote" or "LocfdrControl".
#' @param plot Logical. Default=FALSE, means the plot should be drawn.
#' @param nburn Default=50.
#' @param nskip Default=2
#' @param locfdr.alpha Default=0.2
#' @return A vector, where each value is the summarized class label in c(1,2,3) representing "down-regulate", "null", and "up-regulate".
#' @examples
#' \dontrun{
#' data(net)
#' data(test.stat)
#' res=FeatureClassFinder(net,test.stat,niter=100,na.action="NN")
#' classLabelEst=SummaryClassLabel(net,test.stat,res$zTrack,method="MajorVote")
#' print(table(classLabelEst))
#' }
#' @export
SummaryClassLabel=function(net,test.stat,zTrack,method=c("MajorVote","LocfdrControl"),plot=FALSE,nburn=50,nskip=2,locfdr.alpha=0.2){
  methods=match.arg(method)
  ll=seq(nburn+1,ncol(zTrack),by=nskip)
  zs=zTrack[,ll]
  ngenes=nrow(zs)
  z=vector(length=nrow(zs))
  nn=as.character(seq(1,3))
  if(methods=="LocfdrControl"){
    zpostdist=apply(zs,1,function(xx){
      t=table(xx)
      t=t[nn]
      if(sum(is.na(t))>0){
        t[is.na(t)]=0
      }
      names(t)=nn
      t=t/sum(t)
      t
    })
    zpostdist=t(zpostdist)
    p=zpostdist[,2]
    l0=which(p>locfdr.alpha)
    z[l0]=2
    l1=which(p<=locfdr.alpha)
    if(length(l1)>0){
      negMinusPos=zpostdist[l1,1]-zpostdist[l1,3]
      negind=l1[which(negMinusPos>0)]
      posind=l1[which(negMinusPos<=0)]
      if(length(negind)>0){
        z[negind]=1
      }
      if(length(posind)>0){
        z[posind]=3
      }
    }
  }else{
    z=apply(zs,1,function(xx){
      t=table(xx)
      t=t[nn]
      if(sum(is.na(t))>0){
        t[is.na(t)]=0
      }
      names(t)=nn
      t=t/sum(t)
      which(t==max(t))[1]
    })
  }

  if(plot){
    siggenes=which(z==1|z==3)
    nbrs=lapply(siggenes,function(x){which(net[x,]>0)})
    nbrs=unique(unlist(nbrs))
    allgenes=union(siggenes,nbrs)
    signet=net[allgenes,allgenes]
    siggraph=igraph::graph.adjacency(signet, mode="undirected")

    res=igraph::fastgreedy.community(siggraph)
    memberID=unique(res$membership)
    memberID=sort(memberID)

    ### color pallette:
    ncolor=100
    colfunc=grDevices::colorRampPalette(c("blue", "white","red"))
    mycolor=colfunc(ncolor)
    colvec=seq(floor(range(test.stat,na.rm=TRUE)[1]),ceiling(range(test.stat,na.rm=TRUE)[2]),length=100)
    usecolor=sapply(1:ngenes,function(k){
      if(is.na(test.stat[k])){
        "khaki2"
      }else{
        mycolor[which((test.stat[k]-colvec)<0)[1]]
      }
    })

    ### shape & size pallette
    mysize=c(20,5)
    names(mysize)=c("selected genes","neighbors")
    myshape=c(15:17)
    names(myshape)=c("down-regulated genes","null genes","up-regulated genes")

    for(id in memberID){
      ### plot
      allnode=allgenes[which(res$membership==id)]
      subz=z[allnode]
      sig=allnode[which(subz==1|subz==3)]
      nbr=setdiff(allnode,sig)
      names(allnode)=ifelse(allnode%in% sig,"selected genes","neighbors")

      zz=ifelse(allnode%in% nbr,"neighbors","selected genes")
      zzz=ifelse(subz==1,"down-regulated genes",ifelse(subz==2,"null genes","up-regulated genes"))

      g=network::network(net[allnode,allnode],directed=FALSE)
      if(unique(subz)[1]=="down-regulated genes" & length(unique(subz))==1){
        ggg=GGally::ggnet2(g,color=mycolor[subz],alpha = 0.75,node.size=zz,size.palette =mysize,shape=zzz,shape.palette=myshape,edge.size = 0.2,edge.alpha = 0.5, edge.color = "grey80")
      }else if(unique(subz)[1]=="up-regulated genes" & length(unique(subz))==1){
        ggg=GGally::ggnet2(g,color=mycolor[subz],node.size=zz,size.palette =mysize,shape=zzz,shape.palette=myshape,edge.size = 0.2, edge.alpha = 0.5,edge.color = "grey80")
      }else{
        ggg=GGally::ggnet2(g,color=usecolor[allnode],color.palette = mycolor,alpha = 0.75,node.size=zz,size.palette =mysize,shape=zzz,shape.palette=myshape,edge.alpha = 0.5,edge.size = 0.2, edge.color = "grey80")
      }
      plot(ggg)
    }
  }

  z=z-2
  return(z)
}


