##############################################
#####      calculate basic parameter     #####
##############################################
calcul<-function(time,status,group,t0,tau){

  if(is.null(tau)==TRUE){
    dl<-data.frame(time,status,group)
    lasttime<-c(max(dl$time[dl$group==0]),max(dl$time[dl$group==1]))
    lastevent<-c(max(dl$status[dl$time==lasttime[1] & dl$group==0]),max(dl$status[dl$time==lasttime[2] & dl$group==1]))
    if(all(lastevent==0)) {
      tau=min(lasttime)}
    if(any(lastevent==0) & any(lastevent!=0)){
      tau=max((lasttime[1]*(1-(lastevent[1]>0))),(lasttime[2]*(1-(lastevent[2]>0))))}
    if(all(lastevent!=0)){
      tau=max(lasttime)}
  }
  else{tau=tau}

  ng=table(group)
  lg=labels(ng)$group
  dd=table(time,status,group)
  ddd=dd

  if (dim(dd)[2]==2)
  {ddd=array(0,dim=c(dim(dd)[1],3,2))
  if (sum(status==0)==0){
    ddd[,2:3,]=dd
  }
  if (sum(status==1)==0)
  {warning("There are only competing risks")
    ddd[,1,]=dd[,1,]
    ddd[,3,]=dd[,2,]
  }
  if (sum(status==2)==0)
  {warning("There are no competing risks")
    ddd[,1:2,]=dd
  }
  }

  dd=ddd
  tt=sort(unique(time))
  tt1=c(tt[2:length(tt)],NA)
  nt<-length(which(tt1<=tau))
  deltat=tt1-tt
  dd=dd[1:nt,,]
  deltat=deltat[1:nt]
  tt=tt[1:nt]

  for (i in 1:2){
    dd1=dd[,,i]
    dd2=apply(dd1,1,sum)
    nrisk=ng[i]-cumsum(c(0,dd2[1:(nt-1)]))
    dev=dd1[,2]
    dcr=dd1[,3]
    dcens=dd1[,1]+dcr
    dall=dev+dcr
    si=(nrisk-dall)/nrisk
    s=cumprod(si)
    sminus=c(1,s[1:(length(s)-1)])
    fi=dev/nrisk*sminus
    f=cumsum(fi)
    fcri=dcr/nrisk*sminus
    fcr=cumsum(fcri)
    Ci=(nrisk-dcens)/nrisk
    C=cumprod(Ci)
    Cminus=c(1,C[1:(length(C)-1)])

    if (i==1)
    {nrisk1=nrisk
    si1<-si
    s1=s
    sminus1=sminus
    f1=f
    fcr1=fcr
    C1=C
    C1minus=Cminus
    dall1=dall
    dev1=dev}
    if (i==2)
    {nrisk2=nrisk
    si2<-si
    s2=s
    sminus2=sminus
    f2=f
    fcr2=fcr
    C2=C
    C2minus=Cminus
    dall2=dall
    dev2=dev}
  }
  data.all<-data.frame(nrisk1,nrisk2,s1,s2,sminus1,sminus2,f1,f2,dall1,dall2,dev1,dev2,tau)
  data.all
}


##############################################
#####              ABC test              #####
##############################################
linxu<-function(time,status,group,t0,tau){

  if(is.null(tau)==TRUE){
    dl<-data.frame(time,status,group)
    lasttime<-c(max(dl$time[dl$group==0]),max(dl$time[dl$group==1]))
    lastevent<-c(max(dl$status[dl$time==lasttime[1] & dl$group==0]),max(dl$status[dl$time==lasttime[2] & dl$group==1]))
    if(all(lastevent==0)) {
      tau=min(lasttime)}
    if(any(lastevent==0) & any(lastevent!=0)){
      tau=max((lasttime[1]*(1-(lastevent[1]>0))),(lasttime[2]*(1-(lastevent[2]>0))))}
    if(all(lastevent!=0)){
      tau=max(lasttime)}
  }
  else{tau=tau}
  rho=0.5
  tt=sort(unique(time))
  tt1=c(tt[2:length(tt)],NA)
  nt<-length(which(tt1<=tau))
  deltat=tt1-tt
  deltat=deltat[1:nt]
  tt=tt[1:nt]

  data.all<-calcul(time,status,group,t0=t0,tau=tau)

  if(any(is.na(data.all$s1)==TRUE)){
    data.all$s1[which(is.na(data.all$s1)==TRUE)]<-data.all$s1[min(which(is.na(data.all$s1)==TRUE))-1]
  }
  if(any(is.na(data.all$s2)==TRUE)){
    data.all$s2[which(is.na(data.all$s2)==TRUE)]<-data.all$s2[min(which(is.na(data.all$s2)==TRUE))-1]
  }
  if(any(is.na(data.all$sminus1)==TRUE)){
    data.all$sminus1[which(is.na(data.all$sminus1)==TRUE)]<-data.all$sminus1[min(which(is.na(data.all$sminus1)==TRUE))-1]
  }
  if(any(is.na(data.all$sminus2)==TRUE)){
    data.all$sminus2[which(is.na(data.all$sminus2)==TRUE)]<-data.all$sminus2[min(which(is.na(data.all$sminus2)==TRUE))-1]
  }
  if(any(is.na(data.all$f1)==TRUE)){
    data.all$f1[which(is.na(data.all$f1)==TRUE)]<-data.all$f1[min(which(is.na(data.all$f1)==TRUE))-1]
  }
  if(any(is.na(data.all$f2)==TRUE)){
    data.all$f2[which(is.na(data.all$f2)==TRUE)]<-data.all$f2[min(which(is.na(data.all$f2)==TRUE))-1]
  }

  I<-as.numeric(tt>t0)
  si=(abs(data.all$f1-data.all$f2))*deltat*I
  delta=sum(si)

  var1.1<-(data.all$sminus1*data.all$f1)^2*((data.all$nrisk1-1)/data.all$nrisk1^3)*(data.all$dall1)
  var1.2<-data.all$sminus1^2*(1-2*data.all$f1)*((data.all$nrisk1-1)/data.all$nrisk1^3)*data.all$dev1
  sigmai1<-var1.1+var1.2
  sigmai1[which(is.na(sigmai1)==TRUE)]<-0
  sigma1=cumsum(sigmai1)

  var2.1<-(data.all$sminus2*data.all$f2)^2*((data.all$nrisk2-1)/data.all$nrisk2^3)*(data.all$dall2)
  var2.2<-data.all$sminus2^2*(1-2*data.all$f2)*((data.all$nrisk2-1)/data.all$nrisk2^3)*data.all$dev2
  sigmai2<-var2.1+var2.2
  sigmai2[which(is.na(sigmai2)==TRUE)]<-0
  sigma2=cumsum(sigmai2)

  Edelta<-sum(sqrt((2/pi)*(sigma1+sigma2))*deltat*I)
  Vardelta1<-sum((1-2/pi)*(sigma1+sigma2)*(deltat)^2*I)
  Vardelta2<-0
  t<-c(tt[1],tt1[1:nt])
  I<-c(I[1],I)
  for(m in 1:(length(t)-2)){
    for(n in (m+1):(length(t)-1)){
      Vardelta2=Vardelta2+2*rho*I[m]*I[n]*(t[m+1]-t[m])*(t[n+1]-t[n])*(1-2/pi)*sqrt((sigma1[m]+sigma2[m])*(sigma1[n]+sigma2[n]))
    }
  }
  deltastar<-(delta-Edelta)/sqrt(Vardelta1+Vardelta2)
  Pvalue<-2*(1-pnorm(abs(deltastar)))
  result<-data.frame(t0,tau,delta,"var(delta)"=(Vardelta1+Vardelta2)[[1]],statistic=deltastar[[1]],Pvalue=Pvalue[[1]])
  result
}


ABC<-function(time,status,group,t0=0,tau=NULL){

  if(is.null(tau)==TRUE){
    dl<-data.frame(time,status,group)
    lasttime<-c(max(dl$time[dl$group==0]),max(dl$time[dl$group==1]))
    lastevent<-c(max(dl$status[dl$time==lasttime[1] & dl$group==0]),max(dl$status[dl$time==lasttime[2] & dl$group==1]))
    if(all(lastevent==0)) {
      tau=min(lasttime)}
    if(any(lastevent==0) & any(lastevent!=0)){
      tau=max((lasttime[1]*(1-(lastevent[1]>0))),(lasttime[2]*(1-(lastevent[2]>0))))}
    if(all(lastevent!=0)){
      tau=max(lasttime)}
  }
  else{tau=tau}

  if(t0>tau) stop("tau should be greater than t0")
  if(t0<0) stop("t0 should be equal to or greater than 0")
  if(tau<=unlist(by(time,status,min))[[2]]) stop("tau should be greater than the minimum of event time")
  if(tau>max(unlist(by(time,status,max)))) stop("tau should be equal to or smaller than the maximum of last time")

  ng=table(group)
  lg=labels(ng)$group
  dd=table(time,status,group)

  if (dim(dd)[3]!=2) stop("This test is only for two groups")
  if (dim(dd)[2]>3) stop("All competing risks should be grouped as code 2")
  if (dim(dd)[2]<2) stop("Either all observations are censored or \n there is only one type of event and no censor observations")

  linxu(time,status,group,t0,tau)
}


##############################################
#####       ABC test (permutation)       #####
##############################################
ABC.perm<-function(time,status,group,t0=0,tau=NULL,nperm=1000,seed=12345,bias=FALSE){

  if(is.null(tau)==TRUE){
    dl<-data.frame(time,status,group)
    lasttime<-c(max(dl$time[dl$group==0]),max(dl$time[dl$group==1]))
    lastevent<-c(max(dl$status[dl$time==lasttime[1] & dl$group==0]),max(dl$status[dl$time==lasttime[2] & dl$group==1]))
    if(all(lastevent==0)) {
      tau=min(lasttime)}
    if(any(lastevent==0) & any(lastevent!=0)){
      tau=max((lasttime[1]*(1-(lastevent[1]>0))),(lasttime[2]*(1-(lastevent[2]>0))))}
    if(all(lastevent!=0)){
      tau=max(lasttime)}
  }
  else{tau=tau}

  if(t0>tau) stop("tau should be greater than t0")
  if(t0<0) stop("t0 should be equal to or greater than 0")
  if(tau<=unlist(by(time,status,min))[[2]]) stop("tau should be greater than the minimum of event time")
  if(tau>max(unlist(by(time,status,max)))) stop("tau should be equal to or smaller than the maximum of last time")

  dd=table(time,status,group)

  if (dim(dd)[3]!=2) stop("This test is only for two groups")
  if (dim(dd)[2]>3) stop("All competing risks should be grouped as code 2")
  if (dim(dd)[2]<2) stop("Either all observations are censored or \n there is only one type of event and no censor observations")

  detal_detalstar<-function(data,i){
    d.d.d<-data[i,]
    time=d.d.d$time
    cens=d.d.d$status
    group=data$group

    linxu(time,cens,group,t0,tau)$delta
  }

  cifdata<-data.frame(time,status,group)
  set.seed(seed)
  boot.result<-boot(data=cifdata,statistic=detal_detalstar,R=nperm,sim="permutation")
  B<-boot.result$R

  delta<-boot.result$t0[[1]]
  variance<-apply(boot.result$t,2,var)
  pval<-(sum(abs(boot.result$t)>abs(delta)))/B

  if (bias==TRUE) {
    bias<-(mean(boot.result$t)-boot.result$t0)[[1]]
    result<-data.frame(t0,tau,delta,"var(delta)"=variance,bias,Pvalue=pval)}
  else {result<-data.frame(t0,tau,delta,"var(delta)"=variance,Pvalue=pval)}
  message("The permutation resampling times =", nperm, "\n\n")
  result
}


##############################################
#####        ABC test (two-stage)        #####
##############################################
twostage<-function(time,status,group,nboot,alpha,seed){
  ds<-data.frame(time,status,group)
  out.dif1<-CIFsm(ds,method="dif",pp=0,qq=0)
  PT1<-out.dif1$ave/out.dif1$avese
  PP1<-out.dif1$avepval
  a1<-(1-sqrt(1-alpha))/2
  PT.CT1<-abs(qnorm(a1))

  out.dif2<-linxu(time,status,group,t0=0,tau=NULL)
  PT2<-out.dif2$statistic
  PP2.1<-out.dif2$Pvalue

  PT11<-PT22<-c()
  data<-data.frame(time,status)
  group1<-c(rep(0,length(group[group==0])),rep(1,length(group[group==1])))
  set.seed(seed)
  for(i in 1:nboot){
    h<-sample(seq(1:length(time)),replace=T)
    dd<-data[h,]
    ds11<-data.frame(time=dd$time,status=dd$status,group=group1)
    out.dif11<-CIFsm(ds11,method="dif",pp=0,qq=0)
    PT11[i]<-out.dif11$ave/out.dif11$avese
    out.dif22<-linxu(ds11$time,ds11$status,ds11$group,t0=0,tau=NULL)
    PT22[i]<-out.dif22$statistic
  }
  PTT2<-PT22[abs(PT11)<PT.CT1]
  PT.CT2<-sort(PTT2)[ceiling(length(PTT2)*(1-a1))]
  abjust.P<-sum(abs(PTT2)>PT2)/length(PTT2)
  PP2.2<-sum(abs(PT22)>PT2)/length(PT22)

  if(PP1<=a1*2){
    pvalue<-PP1
  }
  if(PP1>a1*2){
    pvalue<-a1*2+abjust.P*(1-a1*2)
  }
  result<-data.frame(method=c("Li","ABC","Two-stage"),
                     Pvalue=c(PP1,PP2.1,pvalue))
  result
}

ABC.ts<-function(time,status,group,nboot=1000,alpha=0.05,seed=12345){

  dd=table(time,status,group)

  if (dim(dd)[3]!=2) stop("This test is only for two groups")
  if (dim(dd)[2]>3) stop("All competing risks should be grouped as code 2")
  if (dim(dd)[2]<2) stop("Either all observations are censored or \n there is only one type of event and no censor observations")

  message("The bootstrap resampling times =", nboot, "\n\n")
  twostage(time,status,group,nboot,alpha,seed)

}

##############################################
#####         ABC test (combined)        #####
##############################################
combined<-function(time,status,group,nboot,seed){
  ds<-data.frame(time,status,group)
  out.dif1<-CIFsm(ds,method="dif",pp=0,qq=0)
  Stat1<-out.dif1$ave/out.dif1$avese
  P1<-out.dif1$avepval
  out.dif2<-linxu(time,status,group,t0=0,tau=NULL)
  Stat2<-out.dif2$statistic
  P2<-out.dif2$Pvalue
  O.Stat<-max(abs(Stat1),abs(Stat2))

  N.Stat<-Stat11<-Stat22<-c()
  data<-data.frame(time,status)
  group1<-c(rep(1,length(group[group==1])),rep(0,length(group[group==0])))
  set.seed(seed)
  for(i in 1:nboot){
    h<-sample(seq(1:length(time)),replace=T)
    dd<-data[h,]
    ds11<-data.frame(time=dd$time,status=dd$status,group=group1)
    out.dif11<-CIFsm(ds11,method="dif",pp=0,qq=0)
    Stat11[i]<-out.dif11$ave/out.dif11$avese
    out.dif22<-linxu(ds11$time,ds11$status,ds11$group,t0=0,tau=NULL)
    Stat22[i]<-out.dif22$statistic
    N.Stat[i]<-max(abs(Stat11[i]),abs(Stat22[i]))
  }

  boot.P<-sum(N.Stat>=O.Stat,rm.na=T)/nboot

  result<-data.frame(method=c("Li","ABC","Combined"),
                     statistic=c(Stat1,Stat2,O.Stat),
                     Pvalue=c(P1,P2,boot.P))
  result
}


ABC.comb<-function(time,status,group,nboot=1000,seed=12345){

  dd=table(time,status,group)

  if (dim(dd)[3]!=2) stop("This test is only for two groups")
  if (dim(dd)[2]>3) stop("All competing risks should be grouped as code 2")
  if (dim(dd)[2]<2) stop("Either all observations are censored or \n there is only one type of event and no censor observations")

  message("The bootstrap resampling times =", nboot, "\n\n")
  combined(time,status,group,nboot,seed)
}


##############################################
#####  fix point test(Gaynor and Aalen)  #####
##############################################
Fix.tp<-function(time,status,group,timepoint){

  dl<-data.frame(time,status,group)
  lasttime<-c(max(dl$time[dl$group==0]),max(dl$time[dl$group==1]))
  lastevent<-c(max(dl$status[dl$time==lasttime[1] & dl$group==0]),max(dl$status[dl$time==lasttime[2] & dl$group==1]))
    if(all(lastevent==0)) {
      tau=min(lasttime)}
    if(any(lastevent==0) & any(lastevent!=0)){
      tau=max((lasttime[1]*(1-(lastevent[1]>0))),(lasttime[2]*(1-(lastevent[2]>0))))}
    if(all(lastevent!=0)){
      tau=max(lasttime)}

  ng=table(group)
  lg=labels(ng)$group
  dd=table(time,status,group)
  ddd=dd

  if (dim(dd)[2]==2)
  {ddd=array(0,dim=c(dim(dd)[1],3,2))
  if (sum(status==0)==0){
    ddd[,2:3,]=dd
  }
  if (sum(status==1)==0)
  { ddd[,1,]=dd[,1,]
  ddd[,3,]=dd[,2,]
  }
  if (sum(status==2)==0)
  {ddd[,1:2,]=dd}
  }

  dd=ddd

  tt=sort(unique(time))
  tt1=c(tt[2:length(tt)],NA)
  nt<-length(which(tt1<=tau))
  deltat=tt1-tt
  dd=dd[1:nt,,]
  deltat=deltat[1:nt]
  tt=tt[1:nt]

  data.all<-calcul(time,status,group,t0=0,tau)

  if(any(is.na(data.all$sminus1)==TRUE)){
    data.all$sminus1[which(is.na(data.all$sminus1)==TRUE)]<-data.all$sminus1[min(which(is.na(data.all$sminus1)==TRUE))-1]
  }

  if(any(is.na(data.all$sminus2)==TRUE)){
    data.all$sminus2[which(is.na(data.all$sminus2)==TRUE)]<-data.all$sminus2[min(which(is.na(data.all$sminus2)==TRUE))-1]
  }

  Var11<-c()
  Var1<-matrix(0,nrow(data.all),nrow(data.all))
  sumi1<-c(0,cumsum(data.all$dall1/(data.all$nrisk1*(data.all$nrisk1-data.all$dall1))))
  Var22<-c()
  Var2<-matrix(0,nrow(data.all),nrow(data.all))
  sumi2<-c(0,cumsum(data.all$dall2/(data.all$nrisk2*(data.all$nrisk2-data.all$dall2))))
  for(i in 1:nrow(data.all)){
    for(j in i:nrow(data.all)){
      if(i==j){
        deltakl1<-1
        deltakl2<-1
      }
      if(i!=j){
        deltakl1<-0
        deltakl2<-0
      }
      vv11<-((data.all$dev1[i]/data.all$nrisk1[i])*data.all$sminus1[i])*((data.all$dev1[j]/data.all$nrisk1[j])*data.all$sminus1[j])
      vv21<-((data.all$dev2[i]/data.all$nrisk2[i])*data.all$sminus2[i])*((data.all$dev2[j]/data.all$nrisk2[j])*data.all$sminus2[j])
      Var1[i,j]<-vv11*(deltakl1/data.all$dev1[i]-1/data.all$nrisk1[i]+sumi1[i])
      if(data.all$dev1[i]==0 | data.all$dev1[j]==0 ){Var1[i,j]<-0}
      Var2[i,j]<-vv21*(deltakl2/data.all$dev2[i]-1/data.all$nrisk2[i]+sumi2[i])
      if(data.all$dev2[i]==0 | data.all$dev2[j]==0 ){Var2[i,j]<-0}
    }
    Var11[i]<-sum(2*Var1[1:i,1:i])-sum(diag(Var1[1:i,1:i]))
    Var22[i]<-sum(2*Var2[1:i,1:i])-sum(diag(Var2[1:i,1:i]))
  }
  sigma1<-Var11
  sigma2<-Var22
  tt<-c(0,tt)[1:length(sigma1)]

  gaynor.var1<-sigma1[max(which(tt<=timepoint))]
  gaynor.var2<-sigma2[max(which(tt<=timepoint))]

  fit<-timepoints(cuminc(time,status,group,cencode=0),timepoint)
  CIF1<-fit$est[1]
  CIF2<-fit$est[2]
  Aalen.var1<-fit$var[1]
  Aalen.var2<-fit$var[2]

  line<-((CIF1-CIF2)^2)/(gaynor.var1+gaynor.var2)
  line.P<-1-pchisq(line,df=1)
  log<-((log(CIF1)-log(CIF2))^2)/(gaynor.var1/(CIF1^2)+gaynor.var2/(CIF2^2))
  log.P<-1-pchisq(log,df=1)
  cloglog<-((log(-log(CIF1))-log(-log(CIF2)))^2)/((gaynor.var1/(CIF1^2))/((log(CIF1))^2)+(gaynor.var2/(CIF2^2))/((log(CIF2))^2))
  cloglog.P<-1-pchisq(cloglog,df=1)
  v1<-(gaynor.var1/CIF1)/(4*(1-CIF1))
  v2<-(gaynor.var2/CIF2)/(4*(1-CIF2))
  arcsin<-((asin(sqrt(CIF1))-asin(sqrt(CIF2)))^2)/(v1+v2)
  arcsin.P<-1-pchisq(arcsin,df=1)
  logit<-((log(CIF1/(1-CIF1))-log(CIF2/(1-CIF2)))^2)/((gaynor.var1/(CIF1^2))/((1-CIF1)^2)+(gaynor.var2/(CIF2^2))/((1-CIF2)^2))
  logit.P<-1-pchisq(logit,df=1)
  Gaynor<-data.frame(method=c("Line","Log","Cloglog","Arcsin-square","Logit"),
                     est0=rep(CIF1,5),var0=rep(gaynor.var1,5),
                     est1=rep(CIF2,5),var1=rep(gaynor.var2,5),
                     statistic=c(line,log,cloglog,arcsin,logit),
                     Pvalue=c(line.P,log.P,cloglog.P,arcsin.P,logit.P))

  line<-((CIF1-CIF2)^2)/(Aalen.var1+Aalen.var2)
  line.P<-1-pchisq(line,df=1)
  log<-((log(CIF1)-log(CIF2))^2)/(Aalen.var1/(CIF1^2)+Aalen.var2/(CIF2^2))
  log.P<-1-pchisq(log,df=1)
  cloglog<-((log(-log(CIF1))-log(-log(CIF2)))^2)/((Aalen.var1/(CIF1^2))/((log(CIF1))^2)+(Aalen.var2/(CIF2^2))/((log(CIF2))^2))
  cloglog.P<-1-pchisq(cloglog,df=1)
  v1<-(Aalen.var1/CIF1)/(4*(1-CIF1))
  v2<-(Aalen.var2/CIF2)/(4*(1-CIF2))
  arcsin<-((asin(sqrt(CIF1))-asin(sqrt(CIF2)))^2)/(v1+v2)
  arcsin.P<-1-pchisq(arcsin,df=1)
  logit<-((log(CIF1/(1-CIF1))-log(CIF2/(1-CIF2)))^2)/((Aalen.var1/(CIF1^2))/((1-CIF1)^2)+(Aalen.var2/(CIF2^2))/((1-CIF2)^2))
  logit.P<-1-pchisq(logit,df=1)
  Aalen<-data.frame(method=c("Line","Log","Cloglog","Arcsin-square","Logit"),
                    est0=rep(CIF1,5),var0=rep(Aalen.var1,5),
                    est1=rep(CIF2,5),var1=rep(Aalen.var2,5),
                    statistic=c(line,log,cloglog,arcsin,logit),
                    Pvalue=c(line.P,log.P,cloglog.P,arcsin.P,logit.P))

  list(Gaynor=Gaynor,Aalen=Aalen)
}


fixpoint<-function(time,status,group,timepoint,type=1){

  if(sum(timepoint<=min(unlist(by(time,status,min))))==1) stop("timepoint should be greater than the minimum of event time")
  if(sum(timepoint>=max(unlist(by(time,status,max))))==1) stop("timepoint should be smaller than the maximum of last time")

  dd=table(time,status,group)

  if (dim(dd)[3]!=2) stop("This test is only for two groups")
  if (dim(dd)[2]>3) stop("All competing risks should be grouped as code 2")
  if (dim(dd)[2]<2) stop("Either all observations are censored or \n there is only one type of event and no censor observations")

  if (timepoint>unlist(by(time,status,max))[[2]] || timepoint<unlist(by(time,status,min))[[2]])
  {stop("timepoint should be greater than the minimum of event time and smaller than the maximum of last time")}

  if (sum(status==1)==0)
  {warning("There are only competing risks")}
  if (sum(status==2)==0)
  {warning("There are no competing risks")
  }

  message("The estimation of CIFs at the timepoint (timepoint) =", timepoint, "\n\n")
  if (type==1){result<-Fix.tp(time,status,group,timepoint)$Gaynor}
  if (type==2){result<-Fix.tp(time,status,group,timepoint)$Aalen}
  result
}


##############################################
#####         plot ABC function          #####
##############################################
ABC.plot<-function(time,status,group,tau=NULL,max.x=NULL,max.y=NULL,col=c(1,1,8),lwd=c(3,3),lty=c(1,2),
                   lab.x="Time", lab.y="CIF", cex.main=1.5, cex.lab=1.5, cex.axis=1.5){

  if(is.null(tau)==TRUE){
    dl<-data.frame(time,status,group)
    lasttime<-c(max(dl$time[dl$group==0]),max(dl$time[dl$group==1]))
    lastevent<-c(max(dl$status[dl$time==lasttime[1] & dl$group==0]),max(dl$status[dl$time==lasttime[2] & dl$group==1]))
    if(all(lastevent==0)) {
      tau=min(lasttime)}
    if(any(lastevent==0) & any(lastevent!=0)){
      tau=max((lasttime[1]*(1-(lastevent[1]>0))),(lasttime[2]*(1-(lastevent[2]>0))))}
    if(all(lastevent!=0)){
      tau=max(lasttime)}
  }
  else{tau=tau}

  if (is.null(max.y)==TRUE) {ymax=1}
  else  {ymax=max.y}
  if (is.null(max.x)==TRUE) {xmax=max(time)}
  else  {xmax=max.x}

  if (ymax>1) {stop("The cumulative functions should be equal to or smaller than 1")}
  if (ymax<0) {stop("The cumulative functions should be larger than 0")}
  if (xmax<0) {stop("The Time should be larger than 0")}

  cr1<-cuminc(time,status,group,cencode=0)

  fit1<-data.frame(t=cr1$'0 1'$time,e=cr1$'0 1'$est)
  fit2<-data.frame(t=cr1$'1 1'$time,e=cr1$'1 1'$est)
  x1<-fit1$t
  y1<-fit1$e
  x2<-fit2$t
  y2<-fit2$e

  plot(x1,y1,xlim=c(0,xmax),ylim=c(0,ymax),type="n",xlab=lab.x,ylab=lab.y,
       cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis)
  x11<-c(x1[x1<=tau],tau)
  y11<-c(y1[x1<=tau],y1[max(which(x1<=tau))])
  x22<-c(x2[x2<=tau],tau)
  y22<-c(y2[x2<=tau],y2[max(which(x2<=tau))])

  x<-c(x11,rev(x22))
  y<-c(y11,rev(y22))
  polygon(x,y,col=col[3],border =NA)
  lines(x1,y1,type="s",col=col[1],lwd=lwd[1],lty=lty[1])
  lines(x2,y2,type="s",col=col[2],lwd=lwd[2],lty=lty[2])

}


##############################################
#####           data generation          #####
##############################################
rpweibull<-function (n,A1,A2,t,scale){
  re <- rweibull(n,A1,scale)
  ret1<- re[which(re <= t)]

  ret2<-NULL
  while(length(ret1)<n){
    ind<-n-length(ret1)
    re <- rweibull(ind,A2,scale)
    success <- which(re > t)
    ret2<- re[success]
    ret1<-c(ret1,ret2)
  }
  ret1
}

set.seed(20200502)
obs.time1<-rpweibull(100,A1=3,A2=2,t=2.5,scale=2)
obs.status1<-rbinom(100,size=1,prob=0.7)
obs.status1<-ifelse(obs.status1==0,2,1)
time1<-obs.time1
status1<-obs.status1
group1<-rep(0,100)

set.seed(20200502+10000)
obs.time2<-rpweibull(100,A1=0.6,A2=2,t=2.5,scale=2)
obs.status2<-rbinom(100,size=1,prob=0.7)
obs.status2<-ifelse(obs.status2==0,2,1)
time2<-obs.time2
status2<-obs.status2
group2<-rep(1,100)

time<-c(time1,time2)
status<-c(status1,status2)
group<-c(group1,group2)
crossdata<-data.frame(time=round(time,1),status,group)
#save(crossdata,file="C:/Users/Administrator/Desktop/ComparisonCR/data/crossdata.rda")


