Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL ){
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}

dLegendre.model <-function( t, mu, tmin=NULL, tmax=NULL ){
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1]*0 + 1*mu[2];
  if (np.order>=2)
    L <- L + 0.5 * (6 * ti)* mu[3] ;
  if (np.order>=3)
    L <- L +0.5 * (15 * ti ^ 2 - 3)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)* mu[5];
  if (np.order>=5)
    L <- L + 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)*mu[6];
  if (np.order>=6)
    L <- L + (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *ti)* mu[7];
  if (np.order>=7)
    L <- L + (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)* mu[8];
  return(L);
}

Likehood_Legendre_ind <- function(times,para,E){
  
  sum( (E -Legendre.model(t=times,mu=para))^2)
  
}

Likehood_Legendre_cluster <- function(times,para,E){
  
  sum((apply(E, 2, mean)-Legendre.model(t=times,mu=para))^2)
  
} 

smooth.optim_ind <- function(times,para,effect_value){
  
  mean_par <- c()
  
  L <- try(optim(para,Likehood_Legendre_ind,E=effect_value,times=times,method = "BFGS"),TRUE)
  L <- try(optim(L$par,Likehood_Legendre_ind,E=effect_value,times=times,method = "BFGS"),TRUE)
  L <- try(optim(L$par,Likehood_Legendre_ind,E=effect_value,times=times,method = "BFGS"),TRUE)
  if(class(L) == "try-error")
    mean_par<-NA
  else{
    mean_par <- L$par
  mean_par
  }
  
}

smooth.optim_cluster <- function(times,para,effect_value){
  
  mean_par <- c()
  
  L <- optim(para,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  mean_par <- L$par
  mean_par
}

get_cluster <- function(cutree,cluster_num,plast_Value){
  
  allcluster <- list()    
  for (a in 1:cluster_num) {
    
    allcluster[[a]] <- as.data.frame(plast_Value)[which(cutree==c(1:cluster_num)[a]),]
    
  }
  allcluster
} 

varsel1 <- function(X,Y,tt,order=5){
  
  nobs = nrow(X)
  ndim = ncol(X)
  dfo = rep(order-1,ndim)
  index = rep(1:ndim,times=dfo)
  aa2 <- c()
  for(k in 1:ndim){
    aa1 <- c()
    for(j in 1:(order-1)){
      aa <- c()
      for(i in 1:length(tt)){
        aa <- c(aa,Legendre.modelve(tt[i],np.order = j,tmin = min(tt), tmax = max(tt)))
      }
      aa1 <- cbind(aa1,aa*X[,k])
    }
    aa2 <- cbind(aa2,aa1)
  }
  
  
  Xc = scale(aa2,center=TRUE,scale=TRUE)
  n = nrow(Xc)
  
  connect = matrix(0,nrow=ndim,ncol=ndim)
  coefest = matrix(0,nrow=sum(dfo),ncol=ndim)
  regfun = vector("list",length=ndim)
  for(i in 1:ndim)
  {
    inde <- index[-which(index==i)]
    XXc <- Xc[,-which(index==i)]
    
    yc <- Y[,i]-mean(Y[,i])
    
    out1 <- GrpLasso(X=XXc,y=yc,index=inde,lambda=30,crit="BIC")
    var.grp <- out1$var.select  # genes selected
    coef.grp <- out1$coef
    
    ### Adaptive Group Lasso
    index.adp <- index[is.element(index,var.grp)]
    W.adp = sapply(1:length(var.grp),function(j) sqrt(sum(coef.grp[index.adp==var.grp[j]]^2)))
    Xc.adp = Xc[,is.element(index,var.grp)]
    Xcs.adp = scale(Xc.adp,center=F,scale=rep(1/W.adp,times=dfo[var.grp]))
    init.adp = coef.grp/rep(W.adp,times=dfo[var.grp])
    lambda = lambdamax(Xcs.adp,yc,index=index.adp,coef.ini=init.adp,
                       penscale=sqrt,center=F,standardize=F,model=LinReg())*0.7^(1:18)
    out2 = GrpLasso(X=Xc.adp,y=yc,W=W.adp,index=index.adp,ini=coef.grp,
                    lambda=lambda,crit="BIC")
    var.adp = out2$var.select
    coef.adp = out2$coef
    connect[i,var.adp] <-  1
    coefest[is.element(index,var.adp),i] <-  coef.adp
    regfun[[i]] <-  sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    
    
    regfun[[i]] = sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    cat("var=",i,var.adp,"\n")
  }
  return(list(connect=connect,regfun=regfun,coefest=coefest))
}

Legendre.modelve <- function(t, np.order, tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- NA
  L <- 1;
  if (np.order >= 2)
    L <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L <- 0.5 * (15 * ti ^ 2 - 3)
  if (np.order >= 4)
    L <- 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)
  if (np.order >= 5)
    L <- 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L <- (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                       ti)
  if (np.order >= 7)
    L <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                       ti ^ 2 - 35) 
  return(L);
}

GrpLasso <- function(X,y,W=NULL,index,ini=rep(0,ncol(X)),lambda=NULL,crit=c("BIC","EBIC"),center=F){
  if(center==T){
    y = y-mean(y)
    X = scale(X,center=T,scale=F)
  }
  n = nrow(X)
  ind = unique(index)
  p = length(ind)
  dfo = sapply(1:p,function(j) sum(index==ind[j]))
  
  # fit model for a sequence of penalty parameters
  if(!is.null(W)){
    X = scale(X,center=F,scale=rep(1/W,times=dfo))
    ini = ini/rep(W,times=dfo)
  }
  
  # set up the candidates for penalty parameter
  if(is.null(lambda)){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:20)
  }else if(length(lambda)==1){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:lambda)
  }
  
  fit = grplasso(X,y,index=index,lambda=lambda,model=LinReg(),center=F,
                 standardize=F,coef.ini=ini,penscale=sqrt,
                 control=grpl.control(update.hess="lambda",tol=10^-8,trace=0))
  
  # calculate BIC/EBIC
  nlambda = length(lambda)
  rss = sapply(1:nlambda,function(j) sum((y-fit$fitted[,j])^2))
  var.select = sapply(1:nlambda,function(j) unique(index[abs(fit$coef[,j])>0]))
  dfo.lambda = sapply(1:nlambda,function(j) sum(dfo[is.element(ind,var.select[[j]])]))
  if(crit!="BIC" & crit!="EBIC"){
    cat("Error: Criterion not implemented. Reset to BIC!\n")
    crit = "BIC"
  }
  if(crit=="BIC"){
    bic = log(rss)+dfo.lambda*log(n)/n
  }else if(crit == "EBIC"){
    bic = log(rss)+dfo.lambda*log(n)/n+0.5*dfo.lambda*log(p)/n
  }
  
  # select model with smallest value of selection criterion
  id.ss = which.min(bic)
  var.ss = var.select[[id.ss]]
  fit.ss = fit$fitted[,id.ss]
  coef.ss = fit$coef[,id.ss]
  if(!is.null(W)){
    coef.ss = coef.ss*rep(W,times=dfo)
  }
  coef.ss = coef.ss[is.element(index,var.ss)]
  
  return(list(var.select=var.ss,coefficients=coef.ss,fitted=fit.ss,BIC=bic,
              lambda=lambda,fit.org=fit))
}
optim.parallel <- function(connect,effect,n.cores,proc,order,times,nstep){
  
  diag(connect) <- 1
  
  LL <- L_derivative(nt=times,nstep=nstep,order=order)
  
  nx <- dim(effect)[2]
  
  grp <- floor(nx/n.cores)
  grp.i <- c()
  if(n.cores==1){
    grp.i <- c(grp.i,rep(1,nx))
  }else{
    for(ii in 1:n.cores){
      if(ii==n.cores){
        grp.i <- c(grp.i,rep(ii,nx-grp*(ii-1)))
      }else{
        grp.i <- c(grp.i,rep(ii,grp))
      }
    }
  }
  
  grp.ii <- unique(grp.i)
  
  res.list <- mclapply(grp.ii, function(i)
  {
    y.c <- 	which(grp.i==i)
    A <- sapply(y.c, proc, connect=connect,effect=effect,LL=LL,nstep=nstep,order=order,times=times);
    return (unlist(A));
  }, mc.cores=n.cores )
  
  res1 <- do.call("c", res.list)
  res2 <- parallel.data.optim(res1,connect,times)
  return(res2)
}
L_derivative <- function(nt,nstep,order){
  
  stp <- (max(nt)-min(nt))/nstep
  res <- c()
  for(j in 1:nstep){
    
    tg1 <- Legendre.model1((j-1)*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg2 <- Legendre.model1(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg3 <- Legendre.model1(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg4 <- Legendre.model1(j*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tmp1 <- rbind(tg1,tg2,tg3,tg4)
    res <- rbind(res,tmp1)
  }
  res
}
Legendre.model11 <- function(t, np.order,tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- 1;
  if (np.order >= 2)
    L[2] <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L[3] <- 0.5 * (15 * ti ^ 2 - 3) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * 4 * ti ^ 3 - 60 * ti) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                         ti) 
  if (np.order >= 7)
    L[7] <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                          ti ^ 2 - 35)
  return(L);
}
Legendre.model1 <- function(t, np.order,tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- ti;
  if (np.order >= 2)
    L[2] <- 0.5 * (3 * ti^2 - 1) 
  if (np.order >= 3)
    L[3] <- 0.5 * (5 * ti ^ 3 - 3 * ti) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * ti ^ 4 - 30 * ti^2 + 3) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * ti ^ 5 - 70 * ti ^ 3 + 15*ti)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * ti ^ 6 - 315 * ti ^ 4 + 105 * ti^2 -5) 
  
  return(L);
}
ode.optim <- function(y.c,connect,effect,LL,nstep,order,times){
  self <- y.c 
  indexx <- which(connect[y.c,]==1)
  #para <- rep(0.00001,length(indexx)*(order-1))
  para <- rep(0.001,length(indexx)*(order-1))
  #res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
  #            LL=LL,method="BFGS",control=list(maxit=2000,trace=T))
  res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
               LL=LL,control=list(maxit=4000,trace=F))
  cat("Gene=",y.c," ",res$value,"\n")
  A <- ode.sovle.ind(NG=(effect),res$par,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,LL=LL,self=self)
  return(A)
}
fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){
  
  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL,self = self)
  sum((NG[,self]-rowSums(odes))^2) ##最小二???
  
}
ode.sovle.ind <- function(NG,fitpar,nconnect,nt,order,nstep,LL,self){
  stp <- (max(nt)-min(nt))/nstep
  index <- which(nconnect==1)
  
  ind.par <- matrix(fitpar[1:(length(index)*(order-1))],ncol=order-1,byrow=T)
  allrep <- matrix(rep(0,length(index)),nrow=1)
   allrep[which(index==self)] <- NG[1,self]
  nn <- 1
  for(j in 1:nstep){
    tg1 <- (rowSums(t(apply(ind.par,1,"*",LL[nn,])))*NG[j,index])
    tg2 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))*NG[j,index])
    tg3 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))*NG[j,index])
    tg4 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))*NG[j,index])
    tmp <- allrep[j,] +stp*(tg1+2*tg2+2*tg3+tg4)/6
    allrep <- rbind(allrep,tmp)
    nn <- nn + 4
  }
  
  
  self_name <- colnames(NG)[self]
  no_name <- which(colnames(allrep)==self_name)
  allrep[,no_name] <- allrep[,no_name]
  allrep
}

parallel.data.optim <- function(rd,nm,ntt){
  
  nrd <- matrix(rd,nrow=length(ntt))
  nn <- dim(nm)[1]
  ki <- 0
  allist <- list()
  for(i in 1:nn){
    iii <- (which(nm[i,]==1))
    iiil <- length(iii)
    tmp.d <- nrd[,(ki+1):(ki+iiil)]
    if(is.matrix(tmp.d)){
      colnames(tmp.d) <- iii
    }else{
      names(tmp.d) <- iii
    }
    
    allist[[i]] <- tmp.d
    ki <- ki + iiil
  }
  
  return(allist)
}
get_eageInform <- function(eage,cluster){
  
  eage_sati <- c()
  for (n in rownames(cluster)) {
    
    
    if(any(which(eage[,1]==n))){
      a1 <- eage[which(eage[,1]==n),]
      a1.0 <- 0
      a1.1 <- 0
      if(!is.matrix(a1)){
        if(a1[3]==0){
          a1.0 <- 1
          a1.1 <- 0}
        else{
          a1.0 <- 0
          a1.1 <- 1}}
      else{
        for (q in 1:dim(a1)[1]) {
          a2 <- as.matrix(a1)[q,]
          if(a2[3]==0){
            a1.0 <- a1.0+1
            a1.1 <- a1.1+0
          }
          else{
            a1.0 <- a1.0+0
            a1.1 <- a1.1+1
          }
        }
      }
      a3 <- c(a1.0,a1.1)
      
    }
    else{a3 <- c(0,0)}
    
    
    b1.0 <- 0
    b1.1 <- 0
    b1 <- eage[which(eage[,2]==n),]
    if(!is.matrix(b1)){
      if(b1[3]==0){
        b1.0 <- 1
        b1.1 <- 0
      }
      else{
        b1.0 <- 0
        b1.1 <- 1
      }
    }
    else{
      for (w in 1:dim(b1)[1]) {
        b2 <- b1[w,]
        if(b2[3]==0){
          b1.0 <- b1.0+1
          b1.1 <- b1.1+0
        }
        else{
          b1.0 <- b1.0+0
          b1.1 <- b1.1+1}
      }
    }
    b3 <- c(b1.0,b1.1)
    cc <- c(a3,b3)
    eage_sati <- rbind(eage_sati,cc)
    
  }
  eage_sati <- cbind(rownames(cluster),eage_sati)
  colnames(eage_sati) <- c("no","sour_con","sour_pro","tar_con","tar_pro")
  rownames(eage_sati) <- c(1:dim(cluster)[1])
  eage_sati
  
  
}

get_integrate <- function(out,t){
  s <- t[2] - t[1]
  c <- c()
  for (n in 1:length(out[1,])) {
    a <- out[,n]
    b <- c()
    for (i in 1:(length(out[,1])-1)) {
      b <- c(b,a[i] * s)
    }
    c <- cbind(c,b)
  }
  colnames(c) <- colnames(out)
  c/100
}
regasso_moudule <- function(connect1,gene,interaction){
  aaa <- colnames(connect1)
  diag(connect1) <- 0
  ng <- dim(gene)[1]
  allcor <- list()
  for(i in 1:ng){
    a1 <- gene[i,]
    nng1 <- (interaction[[i]])
    if(!is.matrix(nng1)){
      next
    }else{
      nng <- as.matrix(nng1[,-which(colnames(nng1)==aaa[i])])
      corr <- c()
      for(j in 1:dim(nng)[2]){
        corr <- c(corr,cor(as.numeric(a1),as.numeric(nng[,j])))
      }
      allcor[[i]] <- corr
    }
  }
  
  for (q in 1:length(allcor)) {
    connect1[q,which(connect1[q,]==1)] <- allcor[[q]]
  }
  
  return(connect1)
}


##解微分方程，求得SNP互作数据
optim_inter <- function(all_cluster_value,ind_Lpar11,norder){
  
  library(splines)
  library(orthogonalsplinebasis)
  library(MASS)
  library(grplasso)
  library(parallel)
  
  ind_Lpar1 <- ind_Lpar11
  clusterAA <- all_cluster_value
  
  
  ind_LparAA <- ind_Lpar1[as.numeric(rownames(clusterAA)),]
  rownames(ind_LparAA) <-  rownames(clusterAA) 
  
  ttt <- seq(1,57,length.out = (dim(ind_LparAA)[1]+5))
  cAA_mean <- apply(ind_LparAA,1, Legendre.model,t=ttt)
  cAA_d <- apply(ind_LparAA,1, dLegendre.model,t=ttt)
  
  #变量选择variable select
  clusterAA_50 <- varsel1(X=cAA_mean,Y=cAA_d,tt=ttt)
  colnames(clusterAA_50$connect) <- as.numeric(rownames(clusterAA))
  rownames(clusterAA_50$connect) <- as.numeric(rownames(clusterAA))
  
  
  
  clusterAA_list <- list()
  for (n in 1:dim(ind_LparAA)[1]) {
    
    ttemp <- list()
    ttemp[[1]] <-  rownames(clusterAA)[n]
    ttemp[[2]] <- names(which(clusterAA_50$connect[n,]==1))
    clusterAA_list[[n]] <- ttemp
  }
  
  #求解微分方程
  tryyAA <- optim.parallel(connect=clusterAA_50$connect,effect=t(all_cluster_value),
                           n.cores=1,proc=ode.optim,order=norder,times=seq(1,57,0.5),nstep=112)
  
  for (e in 1:length(tryyAA)) {
    
    colnames(tryyAA[[e]]) <- rownames(clusterAA)[as.numeric(colnames(tryyAA[[e]]))]
    
  } 
  
  
  
  ##计算相关系数
  inter_effect <- regasso_moudule(connect1=clusterAA_50$connect,gene=all_cluster_value,interaction=tryyAA)
  
  
  
  for (n in 1:length(clusterAA_list)) {
    
    
    nono <- which(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0])[order(as.numeric(names(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0]))))]==0)
    
    clusterAA_list[[n]][[3]] <- colSums(get_integrate(out=tryyAA[[n]],t=seq(1,57,0.5)))
    
    clusterAA_list[[n]][[3]][nono] <- 0
  } 
  
  
  
  ##assign(paste("optim_cluster_",which_cluster,"list",sep=""),clusterAA_list, envir = .GlobalEnv)
  optim_cluster <- list()
  optim_cluster[[1]] <- clusterAA_list
  optim_cluster[[2]] <- tryyAA
  optim_cluster
}


##?γ?cytoscape????ͼ?????ڵ????ߵ?????
out_Netdata_cytoscape <- function(optim_cluster,which_cluster,clusterAA){
  
  clusterAA_list <- optim_cluster[[1]]
  afterAA <- c()
  effect_selfAA <- c()
  for (i in 1:length(clusterAA_list)){
    dep <- clusterAA_list[[i]][[1]]
    ind <- clusterAA_list[[i]][[2]]
    effectAA <- clusterAA_list[[i]][[3]]
    effect_selfAA <- rbind(effect_selfAA,effectAA[which(names(effectAA)==dep)])
    effectAA <-effectAA[-which(names(effectAA)==dep)]
    one <- c()
    for (j in 1:length(ind)) {
      if(effectAA[j] >= 0){
        type <- 1
      }else{
        type <- 0
      }
      one <- rbind(one,c(ind[j],dep,abs(effectAA[j]),type))
    }
    afterAA <- rbind(afterAA,one)
    
  }
  eageAA <- list(source=afterAA[,1],target=afterAA[,2],colour=as.numeric(afterAA[,4]),
                 Id=as.numeric(afterAA[,2]),Weight=as.numeric(afterAA[,3]))
  eage<-cbind(afterAA[,1],afterAA[,2],as.numeric(afterAA[,4]),as.numeric(afterAA[,2]),as.numeric(afterAA[,3]))
  
  new <- eageAA
  new[[1]] <- paste("S",new[[1]],sep = "")
  new[[2]] <- paste("S",new[[2]],sep = "")
  write.csv(new,paste("M",which_cluster,".csv",sep = ""), row.names = F,fileEncoding = "UTF-8")
  
  get_eageInform(eage,cluster=clusterAA) 
  
}


out_Netdata_cytoscape1 <- function(snp=which(clus[[6]][,62]==24),sigqtl1,optim_cluster,which_cluster,clusterAA){
  
  samesnp<-which(snp%in%sigqtl1)
  snpnames<-paste("SNP",snp,sep="")
  snpnames[samesnp]<-paste("QTL",snp[samesnp],sep="")
  
  clusterAA_list <- optim_cluster[[1]]
  afterAA <- c()
  effect_selfAA <- c()
  for (i in 1:length(clusterAA_list)){
    dep <- clusterAA_list[[i]][[1]]
    ind <- clusterAA_list[[i]][[2]]
    effectAA <- clusterAA_list[[i]][[3]]
    effect_selfAA <- rbind(effect_selfAA,effectAA[which(names(effectAA)==dep)])
    effectAA <-effectAA[-which(names(effectAA)==dep)]
    one <- c()
    for (j in 1:length(ind)) {
      if(effectAA[j] >= 0){
        type <- 1
      }else{
        type <- 0
      }
      one <- rbind(one,c(ind[j],dep,abs(effectAA[j]),type))
    }
    afterAA <- rbind(afterAA,one)
    
  }
  eageAA <- list(source=afterAA[,1],target=afterAA[,2],colour=as.numeric(afterAA[,4]),
                 Id=snpnames[as.numeric(afterAA[,2])],Weight=as.numeric(afterAA[,3]))
  eage<-cbind(afterAA[,1],afterAA[,2],as.numeric(afterAA[,4]),as.numeric(afterAA[,2]),as.numeric(afterAA[,3]))
  
  new <- eageAA
  new[[1]] <- paste("S",new[[1]],sep = "")
  new[[2]] <- paste("S",new[[2]],sep = "")
  write.csv(new,paste("M",which_cluster,".csv",sep = ""), quote=FALSE,row.names = F,fileEncoding = "UTF-8")
  
  get_eageInform(eage,cluster=clusterAA) 
  
}

