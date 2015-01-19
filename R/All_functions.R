#####Transferring large statistical pvalue into testing statistics
Transfer=function(pvalue)
{
  rstat=-(qnorm(pvalue))
  return(rstat)
}

#####Defining parameters 
Parameter_Define=function(initial_mu_var,pirhopair,choice)
{
  parameter=list()
  parameter$beta0=10
  parameter$beta1=10
  parameter$tau0=10
  parameter$tau1=2
  parameter$alpha0=parameter$beta0/initial_mu_var$upvar[1]+1
  parameter$alpha1=parameter$beta1/initial_mu_var$upvar[2]+1
  parameter$gamma0=initial_mu_var$upmu[1]
  parameter$gamma1=initial_mu_var$upmu[2]
  parameter$eta0=initial_mu_var$upvar[1]
  parameter$eta1=initial_mu_var$upvar[2]
  parameter$pi0=pirhopair$pi0[choice]
  parameter$pi1=1-pirhopair$pi0[choice]
  parameter$rho0=pirhopair$rho0[choice]
  parameter$rho1=pirhopair$rho1[choice]
  return(parameter)
}

#####Using Keams for Providing Initial Values
Kmeans=function(rstat)
{
  wholeindex=kmeans(rstat,2)$cluster
  mean1=mean(rstat[which(wholeindex==1)])
  mean2=mean(rstat[which(wholeindex==2)])
  if ((mean2)>(mean1)){
    wholeindex<-wholeindex-rep(1,length(rstat))
  }else{
    wholeindex[which(wholeindex==2)]<-0}
  return(wholeindex)
}

#####Defining initial mu and var
Initial_mu_var=function(rstat,wholeindex)
{
  upmu<-rep(0,2)
  upvar<-rep(0,2)
  initial_mu_var=list()
  initial_mu_var$upmu[1]=mean(rstat[which(wholeindex==0)])
  initial_mu_var$upmu[2]=mean(rstat[which(wholeindex==1)])
  initial_mu_var$upvar[1]=var(rstat[which(wholeindex==0)])
  initial_mu_var$upvar[2]=var(rstat[which(wholeindex==1)])
  return(initial_mu_var)
}

#####Generating selecitons of possible pairs of pi and rho
Pi_rho_gen=function(piall,rhoall)
{ 
  pirhopair=list()
  for (j in 1:length(piall))
  {for (k in (j+1):length(rhoall))
  {
    pirhopair$rho0=c(pirhopair$rho0,rep(rhoall[j],length(piall)))
    pirhopair$rho1=c(pirhopair$rho1,rep(rhoall[k],length(piall)))
  }
  }
  for (i in 1:(length(rhoall)*(length(rhoall)-1)/2))
  {pirhopair$pi0=c(pirhopair$pi0,piall[1:length(piall)])
  }
  return(pirhopair)
}


#####Generating Zi by different setting of pi and rho
Generating_Zi=function(net,pirhopair,wholeindex,n=30)
{
  z=wholeindex
  znew=c()
  
  for (j in 1:length(pirhopair$pi0))
  {
    for(i in 1: n)
     
    {
      pro1<-0
      pro0<-0
      for(num1 in 1:length(wholeindex)){
        idx = which(net[num1,]==1)
        pro0=sum(z[idx]==0)
        pro1=sum(z[idx]==1)
                
                
        log0<-log(pirhopair$pi0[j])+(2*pirhopair$rho0[j]*pro0)
        
        log1<-log(1-pirhopair$pi0[j])+(2*pirhopair$rho1[j]*pro1)
        
        
        
        p0<-1/(1+exp(log1-log0))
        p1<-1-p0
        
        z[num1]<-sample(c(0,1),1,prob=c(p0,p1))
      }
   }
    znew=rbind(znew,z)
  }
  return(znew)
}  


#####Using likelyhood for choosing pi and rho
Pi_rho_selecting=function(znew,rstat,piall,rhoall)
{
  mylog<-0
  like<-0
  for(i in 1:(length(piall)*sum(seq(length(rhoall))-1))){
    mu0=mean(rstat[which(znew[i,]==0)])
    mu1=mean(rstat[which(znew[i,]==1)])
    var0=var(rstat[which(znew[i,]==0)])
    var1=var(rstat[which(znew[i,]==1)])
    for(num1 in 1:length(rstat)){
            if(znew[i,num1]==0){
        mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
      }else{
        mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
      }
      
    }
    like=c(like,sum(mylog))
    likenew=like[-1]
  }
  maxlike=max(likenew[!is.na(likenew)])
  choice=which(likenew==maxlike)
  if (length(choice)!=0){choice=sample(choice,1)}
  return(choice)
}
####Gaussian Quadrature Part Integration
Innerfunc=function(zz,rstat,model,
                   ww2=0.8049141,ww1=0.08131284,ww3=0.8049141,ww4=0.08131284,
                   xx1=-1.65068,xx2=-0.5246476,xx3= 0.5246476,xx4= 1.65068){
  if(zz==1){
    ff1=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx1-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    ff2=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx2-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    ff3=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx3-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    ff4=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx4-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    
    outvalue=exp(lgamma(model$parameter$alpha1+0.5)-lgamma(model$parameter$alpha1))/sqrt(pi*model$parameter$beta1)
    
  }else{
    ff1=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx1-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    ff2=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx2-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    ff3=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx3-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    ff4=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx4-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    
    outvalue=exp(lgamma(model$parameter$alpha0+0.5)-lgamma(model$parameter$alpha0))/sqrt(pi*model$parameter$beta0)
    
  }
  ff<-c(ff1,ff2,ff3,ff4)
  integralvalue=outvalue*(c(ww1,ww2,ww3,ww4)%*%ff)
  return(integralvalue)
}

#####giibs for mu ars for var
Gibbsfortheta<-function(initalmu,initalvar,index,model)
{
  uprr=initalmu
  initialvarvar=initalvar
  
  if(index==0){
    
    
    myvar<-rigamma(1,model$parameter$alpha0,model$parameter$beta0)
    
    
    meanforupmu=(myvar*model$parameter$gamma0+model$parameter$eta0^2*uprr)/(myvar+model$parameter$eta0^2)
    varforupmu=myvar*model$parameter$eta0^2/(myvar+model$parameter$eta0^2)
    
    
    mymu<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
    
    
  }else{
    
    
    myvar<-rigamma(1,model$parameter$alpha1,model$parameter$beta1)
    
    
    
    meanforupmu=(myvar*model$parameter$gamma1+model$parameter$eta1^2*uprr)/(myvar+model$parameter$eta1^2)
    varforupmu=myvar*model$parameter$eta1^2/(myvar+model$parameter$eta1^2)
    
    mymu<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
    
    
  }
  
  
  return(c(mymu,myvar))
}



#####Iterations
####Define a list called model which would be included in iterations later
Model=function(wholeindex,net,rstat,parameter,initial_mu_var)
{
  model=list()
  model$wholeindex=wholeindex
  model$net=net
  model$rstat=rstat
  model$parameter=parameter
  model$initial_mu_var=initial_mu_var
  return(model)
}
######DMP1 Iteration functions
Step1_1_Update_gi_zi=function(model,num1)
{
  model$uniqueindex<-sort(unique(model$wholeindex))
  model$originalindex<-model$uniqueindex
  sampleprob<-rep(NA,length(model$uniqueindex)+2)
  
  
  idx = which(model$net[num1,]==1)
  
  pro0=sum((model$wholeindex[idx]<1))
  pro1=sum((model$wholeindex[idx]>0))
  ## calculate mk
  propro0=sum(model$wholeindex<1)
  propro1=sum(model$wholeindex>0)
  
  
  sumforg=0
  uu<-0
  for(gg in model$uniqueindex){
    uu<-uu+1
    sumforg[uu]<-sum(model$wholeindex==gg)-(model$wholeindex[num1]==gg)
    if(gg<1){
      sampleprob[uu+1]=log(model$parameter$pi0)+(2*model$parameter$rho0*pro0-(model$rstat[num1]-model$initial_mu_var$upmu[uu])^2*0.5/model$initial_mu_var$upvar[uu])-0.5*log(model$initial_mu_var$upvar[uu])+log(sumforg[uu])-log(model$parameter$tau0+propro0-1)
      
    }else{
      sampleprob[uu+1]=log(1-model$parameter$pi0)+(2*model$parameter$rho1*pro1-(model$rstat[num1]-model$initial_mu_var$upmu[uu])^2*0.5/model$initial_mu_var$upvar[uu])-0.5*log(model$initial_mu_var$upvar[uu])+log(sumforg[uu])-log(model$parameter$tau1+propro1-1)
    }
  }
  
  test1whole=model$wholeindex
  
  sampleprob[1]= log(model$parameter$pi0)+(2*model$parameter$rho0*pro0)+log(model$parameter$tau0)-log(model$parameter$tau0+propro0-1)+log(Innerfunc(0,model$rstat[num1],model=model))
  sampleprob[length(model$uniqueindex)+2]= log(1-model$parameter$pi0)+(2*model$parameter$rho1*pro1)+log(model$parameter$tau1)-log(model$parameter$tau1+propro1-1)+log(Innerfunc(1,model$rstat[num1],model=model))
  newsampleprob<-0
  for(uu in 1:length(sampleprob)){ newsampleprob[uu]=1/sum(exp(sampleprob-sampleprob[uu]))}
  newsampleprob[which(newsampleprob=="NaN")]=0
  model$origindex<-model$wholeindex[num1]
  model$wholeindex[num1]=sample((min(model$wholeindex)-1):(max(model$wholeindex)+1),1,prob=newsampleprob)
  model$sampleindex=model$wholeindex[num1]  
  return(model)
}

## check whether one cluster is disappear
Step1_2_Check=function(model,num1)
{
  model$myuniqueindex<-sort(unique(model$wholeindex))
  if(sum(model$wholeindex==model$origindex)==0){
    
    deletindex=which(model$originalindex==model$origindex)
    
    model$initial_mu_var$upmu=model$initial_mu_var$upmu[-deletindex]
    model$initial_mu_var$upvar=model$initial_mu_var$upvar[-deletindex]
    if(model$origindex<1){
      model$wholeindex[which(model$wholeindex<model$origindex)]=model$wholeindex[which(model$wholeindex<model$origindex)]+1
    }else{
      
      model$wholeindex[which(model$wholeindex>model$origindex)]=model$wholeindex[which(model$wholeindex>model$origindex)]-1
      
    }
    model$uniqueindex<-sort(unique(model$wholeindex))
    if((model$sampleindex %in% model$originalindex)==0){
      if(model$sampleindex<1){model$uniqueindex=model$uniqueindex[-1]}else{model$uniqueindex=model$uniqueindex[-length(model$uniqueindex)]}
      
    }
    
    
  }
  
  
  test1mu=model$initial_mu_var$upmu
  test1var=model$initial_mu_var$upvar
  test1unique=model$uniqueindex
  return(model)
}


## update mu and variance part
Step_2_Update_mu_var=function(model,num1)
{
  uu=0
  for(gg in (model$uniqueindex)){
    uu=uu+1
    if(gg <1){
      meanforupmu=(model$initial_mu_var$upvar[uu]*model$parameter$gamma0+model$parameter$eta0^2*sum(model$rstat[which(model$wholeindex==gg)]))/(model$initial_mu_var$upvar[uu]+model$parameter$eta0^2*sum(model$wholeindex==gg))
      varforupmu=model$initial_mu_var$upvar[uu]*model$parameter$eta0^2/(model$initial_mu_var$upvar[uu]+model$parameter$eta0^2*sum(model$wholeindex==gg))
      
      model$initial_mu_var$upmu[uu]<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
      all0<-model$parameter$alpha0+sum(model$wholeindex==gg)*0.5
      laa0<-sum((model$rstat[which(model$wholeindex==gg)]-model$initial_mu_var$upmu[uu])^2)*0.5+model$parameter$beta0
      model$initial_mu_var$upvar[uu]<-rigamma(1,all0,laa0)
      
    }else{
      meanforupmu=(model$initial_mu_var$upvar[uu]*model$parameter$gamma1+model$parameter$eta1^2*sum(model$rstat[which(model$wholeindex==gg)]))/(model$initial_mu_var$upvar[uu]+model$parameter$eta1^2*sum(model$wholeindex==gg))
      varforupmu=model$initial_mu_var$upvar[uu]*model$parameter$eta1^2/(model$initial_mu_var$upvar[uu]+model$parameter$eta1^2*sum(model$wholeindex==gg))
      
      model$initial_mu_var$upmu[uu]<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
      all1<-model$parameter$alpha1+sum(model$wholeindex==gg)*0.5
      laa1<-sum((model$rstat[which(model$wholeindex==gg)]-model$initial_mu_var$upmu[uu])^2)*0.5+model$parameter$beta1
      model$initial_mu_var$upvar[uu]<-rigamma(1,all1,laa1)
      
    }
    
  }
  
  if((model$sampleindex %in% model$originalindex)==0){
    if(model$sampleindex==min(model$myuniqueindex)){
      
      results=Gibbsfortheta(model$rstat[num1],model$parameter$beta0/(model$parameter$alpha0+0.5),0,model=model)
      newupmu=results[1]
      newupvar=results[2]
      model$initial_mu_var$upmu<-c(newupmu,model$initial_mu_var$upmu)
      model$initial_mu_var$upvar<-c(newupvar,model$initial_mu_var$upvar)
      
    }else if(model$sampleindex==max(model$myuniqueindex)){
      results=Gibbsfortheta(model$rstat[num1],model$parameter$beta1/(model$parameter$alpha1+0.5),1,model=model)
      newupmu=results[1]
      newupvar=results[2]
      model$initial_mu_var$upmu<-c(model$initial_mu_var$upmu,newupmu)
      model$initial_mu_var$upvar<-c(model$initial_mu_var$upvar,newupvar)
      
    }
  }
  return(model)
}

## switch the label
Step_3_switch_the_label=function(model,num1)
{
  model$uniqueindex<-sort(unique(model$wholeindex))
  test2whole=model$wholeindex
  test2mu=model$initial_mu_var$upmu
  test2var=model$initial_mu_var$upvar
  
  if(length(model$uniqueindex)>1){
    orderindex=rank(model$initial_mu_var$upmu)
    ssupmu=0
    ssupvar=0
    sswholeindex<-rep(NA,length(model$rstat))
    for(kk in 1:length(model$uniqueindex)){
      ssupmu[kk]=model$initial_mu_var$upmu[which(orderindex==kk)]
      ssupvar[kk]=model$initial_mu_var$upvar[which(orderindex==kk)]
      dudu=which(orderindex==kk)
      sswholeindex[which(model$wholeindex==model$uniqueindex[dudu])]=model$uniqueindex[kk]
      
      
    }
    model$initial_mu_var$upmu=ssupmu
    model$initial_mu_var$upvar=ssupvar
    model$wholeindex=sswholeindex
  }
  return(model)
}

Inte_Distance=function(i,mclust)
{
  if (mclust$parameter$variance$modelName=="E"){mclust$parameter$variance$sigmasq=rep(mclust$parameter$variance$sigmasq,length(mclust$parameter$mean))}
  distance=dnorm(mclust$parameter$mean[i],mclust$parameter$mean[i],sd=sqrt(mclust$parameter$variance$sigmasq[i]+mclust$parameter$variance$sigmasq[i]))-dnorm(mclust$parameter$mean[i],mclust$parameter$mean[i+1],sd=sqrt(mclust$parameter$variance$sigmasq[i]+mclust$parameter$variance$sigmasq[i+1]))  -dnorm(mclust$parameter$mean[i+1],mclust$parameter$mean[i],sd=sqrt(mclust$parameter$variance$sigmasq[i+1]+mclust$parameter$variance$sigmasq[i]))  +dnorm(mclust$parameter$mean[i+1],mclust$parameter$mean[i+1],sd=sqrt(mclust$parameter$variance$sigmasq[i+1]+mclust$parameter$variance$sigmasq[i+1]))
return(distance)
}

Inte_Distance_DPdensity=function(i,mclust)
{
  
  distance=dnorm(mclust$parameter$mean[i],mclust$parameter$mean[i],sd=sqrt(mclust$parameter$variance$sigmasq[i]+mclust$parameter$variance$sigmasq[i]))
  -dnorm(mclust$parameter$mean[i],mclust$parameter$mean[i+1],sd=sqrt(mclust$parameter$variance$sigmasq[i]+mclust$parameter$variance$sigmasq[i+1]))
  -dnorm(mclust$parameter$mean[i+1],mclust$parameter$mean[i],sd=sqrt(mclust$parameter$variance$sigmasq[i+1]+mclust$parameter$variance$sigmasq[i]))
  +dnorm(mclust$parameter$mean[i+1],mclust$parameter$mean[i+1],sd=sqrt(mclust$parameter$variance$sigmasq[i+1]+mclust$parameter$variance$sigmasq[i+1]))
  return(distance)
}


######HODCMclust
HODCMclust=function(mclust,rstat)
{
  while (length(mclust$parameter$variance$sigmasq)!=length(unique(mclust$parameter$variance$sigmasq))){
    for (i in 1:length(mclust$parameter$mean)){
      for (j in 1:length(mclust$parameter$mean)){
        if(i!=j){if(mclust$parameter$variance$sigmasq[i]==mclust$parameter$variance$sigmasq[j]){mclust$parameter$variance$sigmasq[j]=mclust$parameter$variance$sigmasq[j]+0.0001}}
      }
    }
  }
  
  while (length(mclust$parameter$pro)!=length(unique(mclust$parameter$pro))){
    for (i in 1:length(mclust$parameter$mean)){
      for (j in 1:length(mclust$parameter$mean)){
        if(i!=j){if(mclust$parameter$pro[i]==mclust$parameter$pro[j]){mclust$parameter$pro[j]=mclust$parameter$pro[j]+0.0001}}
        
      }
    }
  }
  
  if (length(mclust$parameter$mean)==1) {print("warning: the input is not appropriate for mclust since only one cluster was detected by the function Mclust" )}
  if (length(mclust$parameter$mean)==1) break
  ###Step1 find the min distance
  if (length(mclust$parameter$mean)==2) {
    hodcmclust=list()
    hodcmclust$mean=unique(mclust$parameter$mean)
    hodcmclust$pro=unique(mclust$parameter$pro)
    hodcmclust$variance=unique(mclust$parameter$variance$sigmasq[!is.na(mclust$parameter$variance$sigmasq)])
  }else{
    repeat{
   
      distance_all=0
      mclust$parameter$mean=unique(mclust$parameter$mean)
      mclust$parameter$pro=unique(mclust$parameter$pro)
      mclust$parameter$variance$sigmasq=unique(mclust$parameter$variance$sigmasq[!is.na(mclust$parameter$variance$sigmasq)])
      for (i in 1:(length(mclust$parameter$mean)-1))
      {
        distance=Inte_Distance(i,mclust)
        distance_all=c(distance_all,distance)
      }
      lmin=which(distance_all[-1]==min(distance_all[-1]))
      
      if (length(lmin)!=1){lmin=sample(lmin,1)}
      
      for (l in 1:(length(mclust$parameter$mean)-1))
      {
        if (l<lmin){mclust$parameter$mean[l]=mclust$parameter$mean[l]
                    mclust$parameter$pro[l]=mclust$parameter$pro[l]
        
        }else if (l==lmin){mclust$parameter$mean[l]=mclust$parameter$mean[l]*mclust$parameter$pro[l]/(mclust$parameter$pro[l]+mclust$parameter$pro[l+1])+mclust$parameter$mean[l+1]*mclust$parameter$pro[l+1]/(mclust$parameter$pro[l]+mclust$parameter$pro[l+1])
                          mclust$parameter$pro[l]=mclust$parameter$pro[l]+ mclust$parameter$pro[l+1]
                          k=lmin
                           repeat{
                             k=k+1                           
                             if (length(mclust$classification[which(mclust$classification==k)])!=0){mclust$classification[which(mclust$classification==k)]=lmin
                                                                                                    break}
                           }
                          
                          if (mclust$parameter$variance$modelName!="E"){mclust$parameter$variance$sigmasq[l]=var(rstat[which(mclust$classification==l)])}
        }else if (l>lmin){mclust$parameter$mean[l]=mclust$parameter$mean[l+1]
                         mclust$parameter$variance$sigmasq[l]=mclust$parameter$variance$sigmasq[l+1]
                         mclust$parameter$pro[l]=mclust$parameter$pro[l+1]
        }
        
      }
      
      if (lmin==(length(mclust$parameter$mean)-1)){mclust$parameter$mean=mclust$parameter$mean[-length(mclust$parameter$mean)]
                                                   mclust$parameter$pro=mclust$parameter$pro[-length(mclust$parameter$pro)]
                                                   if (mclust$parameter$variance$modelName!="E"){mclust$parameter$variance$sigmasq=mclust$parameter$variance$sigmasq[-length(mclust$parameter$variance$sigmasq)] }}
      if (length(unique(mclust$parameter$mean))==2) break
    }
    hodcmclust=list()
    index=sort(unique(mclust$classification))
    hodcmclust$mean[1]=mean(rstat[which(mclust$classification==index[1])])
    hodcmclust$mean[2]=mean(rstat[which(mclust$classification==index[2])])
    hodcmclust$variance[1]=var(rstat[which(mclust$classification==index[1])])
    hodcmclust$variance[2]=var(rstat[which(mclust$classification==index[2])])
    hodcmclust$pro[1]=length(rstat[which(mclust$classification==index[1])])/length(rstat)
    hodcmclust$pro[2]=length(rstat[which(mclust$classification==index[2])])/length(rstat)
    hodcmclust$classification=mclust$classification}
  return(hodcmclust) 
}


#####Iteration3_Mclust
Iteration3_Mclust<-function(iter,wholeindex,hodcmclust,net,pirhopair,choice,rstat,show.steps,showlikelihood,likelihood.frequency){
  z<-wholeindex
  total<-matrix(rep(0,length(wholeindex)*iter),ncol=length(wholeindex),nrow=iter)
  
  for(jj in 1: iter){
    
    if(jj%%show.steps==0){
      cat("iter: ",jj,"\n")
      flush.console()
    }
    for(num1 in 1:length(wholeindex)){
      ztemp=c()
      pro1<-0
      pro0<-0
      idx = which(net[num1,]==1)
      pro0=sum((z[idx]==0))
      pro1=sum((z[idx]==1))
      if (length(hodcmclust$variance)==1){hodcmclust$variance=c(hodcmclust$variance,hodcmclust$variance)}
     
        log0<-log(pirhopair$pi0[choice])+2*pirhopair$rho0[choice]*pro0-(rstat[num1]-hodcmclust$mean[1])^2/(2*hodcmclust$variance[1])-log(sqrt(2*pi*hodcmclust$variance[1]))
        log1<-log(1-pirhopair$pi0[choice])+2*pirhopair$rho1[choice]*pro1-(rstat[num1]-hodcmclust$mean[2])^2/(2*hodcmclust$variance[2])-log(sqrt(2*pi*hodcmclust$variance[2]))
        
        
        p0<-1/(1+exp(log1-log0))
        p1<-1-p0
        
      if (is.na(p0) | is.na(p1)) {z[num1]=sample(c(1,0),1,prob=c(0.5,0.5))###in case the mu is null
      }else{
      z[num1]<-sample(c(0,1),1,prob=c(p0,p1))}
      
      
        total[jj,num1]<-z[num1]
      
      #ratio[jj,num1]<-total[jj,num1]/jj
      
    }
    
    if(showlikelihood==TRUE){
      if(jj%%likelihood.frequency==0){
        mylog<-0
        mu0=mean(rstat[which(total[jj,]<=0)])
        mu1=mean(rstat[which(total[jj,]>=1)])
        var0=var(rstat[which(total[jj,]<=0)])
        var1=var(rstat[which(total[jj,]>=1)])
        for(num1 in 1:length(rstat)){
          
          if(total[jj,num1]<=0){
            mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
          }else{
            mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
          }
          
        }
        cat("Now for the step:" ,jj, "the log-likelihood value is" ,sum(mylog) , "\n")
        flush.console()
      }
    }
  }
  
  return(total)
}


#####DPdensity
#####DPdensity
DPdensitycluster<-function(v,rstat,DPM.mcmc,DPM.prior){
  results<-list()
  
 
    
    mcmc <- DPM.mcmc
    #prior1 <- list(alpha=0.5,m1=0.2,nu1=100,psiinv1=0.1,k0=100)
    #prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
     #              psiinv2=solve(diag(0.5,1)),
      #             nu1=4,nu2=4,tau1=1,tau2=100)
    fit<-DPdensity(rstat,prior=DPM.prior,mcmc=mcmc,status=TRUE)
    
  for(oo in 1:v){
    nburn <-0
    nsave <-1
    nskip <-0
    ndisplay <-10
    mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
    #prior1 <- list(alpha=0.5,m1=0.2,nu1=100,psiinv1=0.1,k0=100)
    prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)
    fit<-DPdensity(rstat,prior=prior2,mcmc=mcmc,state=fit$state,status=FALSE)
    
    ##calculate cluster
    ss<-fit$state$ss
    num<-fit$state$ncluster
    prop<-table(fit$state$ss)/num
    mu<-fit$state$muclus[1:num]
    sigma<-fit$state$sigmaclus[1:num]
    index<-order(mu,decreasing=T)
    newm<-mu[index]
    news<-sigma[index] 
    newp<-prop[index]
    
    
    if(num==2 | num==1){results[[oo]]<-ss-1}else{
      #hh1<-1/2/sqrt(pi*news[1])+1/2/sqrt(pi*news[2])-2*dnorm(newm[1],newm[2],sd=sqrt(news[1]+news[2]))
      
      #hh2<-1/2/sqrt(pi*news[2])+1/2/sqrt(pi*news[3])-2*dnorm(newm[2],newm[3],sd=sqrt(news[2]+news[3]))
      
      clmatrix<-matrix(rep(0,num*num),ncol=num,nrow=num)
      for(i in 1:num){
        for(j in 1:num){
          clmatrix[i,j]<-dnorm(newm[i],newm[j],sd=sqrt(news[i]+news[j]))
          
        }
      }
      
      
      temp1<-0
      temp2<-0
      com1<-0
      com2<-0
      res1<-0
      res2<-0
      final1<-0
      final2<-0
      latent<-0
      latent[1]<-0
      final1<-c(1,-1,rep(0,num-2))
      final2<-c(0,-1,1,rep(0,num-3))
      
      for(iii in 1:(num-2)){
        dec1<-t(final1)%*%clmatrix%*%final1
        dec2<-t(final2)%*%clmatrix%*%final2
        
        if(dec1>dec2){latent[iii+1]<-1
        }else{
          latent[iii+1]<-0
        }
        if(latent[iii+1]==1){
          ddd<-max(which(latent==0))
          temp1<-rep(0,num)
          temp1[1:ddd]<-1
          com1<-temp1*newp
          com1<-scale(com1,center=F,scale=sum(com1))
          temp1<-rep(0,num)
          temp1[(ddd+1):(iii+1)]<--1
          res1<-temp1*newp
          res1<-scale(res1,center=F,scale=sum(abs(res1)))
          final1<-com1+res1
          temp2<-rep(0,num)
          temp2[(ddd+1):(iii+1)]<--1
          res2<-temp2*newp
          res2<-scale(res2,center=F,scale=sum(abs(res2)))
          com2<-rep(0,num)
          com2[iii+2]<-1
          final2<-com2+res2
        }else{
          temp1<-rep(0,num)
          temp1[1:iii]<-1
          com1<-temp1*newp
          com1<-scale(com1,center=F,scale=sum(com1))
          res1<-rep(0,num)
          res1[iii+1]<--1
          final1<-com1+res1
          final2<-rep(0,num)
          final2[iii+1]<--1
          final2[iii+2]<-1
        }
        
        
      }
      latent[num]<-1
      
      uuu<-max(which(latent==0))
      ss[ss %in% index[(uuu+1):num]]<-0
      ss[ss %in% index[1:uuu]]<-1
      results[[oo]]<-ss
    }
  }
  dpdensitycluster<-do.call(rbind,results)
  return(dpdensitycluster)
}

Iteration3_DPdensity_Par<-function(iter,wholeindex,dpdensitycluster,net,pirhopair,choice,rstat,v,show.steps,n.cores){
  z<-wholeindex
  total<-matrix(rep(0,length(wholeindex)*iter),ncol=length(wholeindex),nrow=iter)
  pro1<-0
  pro0<-0
  kk<-0
  mu0=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==0)])))
  mu1=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==1)])))
  var0=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==0)])))
  var1=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==1)])))
  registerDoParallel( cores=n.cores)
  
  
  ztemp<-foreach (kk = 1:v) %dopar% {
    mu0new=mu0[kk]
    mu1new=mu1[kk]
    if (is.na(var0[kk])){var0[kk]=var1[kk]
    }else if (is.na(var1[kk])){var1[kk]=var0[kk]}
    var0new=var0[kk]
    var1new=var1[kk]
    
    
    for(jj in 1: iter){
      
        if(jj%%show.steps==0){
          cat("iter: ",jj,"\n")
          flush.console()
        }
      for(num1 in 1:length(wholeindex)){
        ztemp=c()
        
        idx = which(net[num1,]==1)
        pro0=sum((z[idx]==0))
        pro1=sum((z[idx]==1))
        
        log0<-log(pirhopair$pi0[choice])+2*pirhopair$rho0[choice]*pro0-(rstat[num1]-mu0new)^2/(2*var0new)-log(sqrt(2*pi*var0new))
        log1<-log(1-pirhopair$pi0[choice])+2*pirhopair$rho1[choice]*pro1-(rstat[num1]-mu1new)^2/(2*var1new)-log(sqrt(2*pi*var1new))
        
        
        p0<-1/(1+exp(log1-log0))
        p1<-1-p0
        
        if (is.na(p0) | is.na(p1)) {z[num1]=sample(c(1,0),1,prob=c(0.5,0.5))###in case the mu is null
        }else{
          z[num1]<-sample(c(0,1),1,prob=c(p0,p1))}
        
        total[jj,num1]<-z[num1]
        
      }#ratio[jj,num1]<-total[jj,num1]/jj
      
      
    }
    matrix=total
    
  }
  results=t(sapply(1:iter, function(it) return(ztemp[[sample(1:v,1)]][it,])))
  stopImplicitCluster()
  return(results)
}

Iteration3_DPdensity<-function(iter,wholeindex,dpdensitycluster,net,pirhopair,choice,rstat,v,show.steps,trace){
  z<-wholeindex
  total<-matrix(rep(0,length(wholeindex)*iter),ncol=length(wholeindex),nrow=iter)
  pro1<-0
  pro0<-0
  
  mu0=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==0)])))
  mu1=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==1)])))
  var0=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==0)])))
  var1=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==1)])))
  for(jj in 1: iter){
   
      if(jj%%show.steps==0){
        cat("iter: ",jj,"\n")
        flush.console()
      }
    for(num1 in 1:length(wholeindex)){
      ztemp=c()
      
      idx = which(net[num1,]==1)
      pro0=sum((z[idx]==0))
      pro1=sum((z[idx]==1))
      
      for (kk in 1:v){
        
        if (is.na(var0[kk])){var0[kk]=var1[kk]
        }else if (is.na(var1[kk])){var1[kk]=var0[kk]}
        
        log0<-log(pirhopair$pi0[choice])+2*pirhopair$rho0[choice]*pro0-(rstat[num1]-mu0[kk])^2/(2*var0[kk])-log(sqrt(2*pi*var0[kk]))
        log1<-log(1-pirhopair$pi0[choice])+2*pirhopair$rho1[choice]*pro1-(rstat[num1]-mu1[kk])^2/(2*var1[kk])-log(sqrt(2*pi*var1[kk]))
        
        p0<-1/(1+exp(log1-log0))
        p1<-1-p0
        if (is.na(p0) | is.na(p1)) {ztemp=ztemp
        }else{
          ztemp<-c(ztemp,sample(c(0,1),1,prob=c(p0,p1)))
        }}
      if (is.null(ztemp)) {ztemp=sample(c(1,0),1,prob=c(0.5,0.5))###in case the mu is null
      }else{
        z[num1]<-sample(ztemp,1,prob=c(rep(1/length(ztemp),length(ztemp))))}
      
      
      total[jj,num1]<-z[num1]
      
      #ratio[jj,num1]<-total[jj,num1]/jj
      
    }
  }
  
  return(total)
}


####True Positive Rate
TPR=function(num,Trace)
{
  new1=Trace[,num]
  ratiosum=length(new1[which(new1>=1)])/(length(new1))
  
  finalratio=mean(ratiosum)
}


