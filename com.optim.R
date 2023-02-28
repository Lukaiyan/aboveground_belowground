com.est1 <- function(dat1,dat2,dat3,dat4,interval=c(1,1),parin2){
  
  y1 <- as.matrix( dat1$pheno)
  y2 <- as.matrix( dat2$pheno)
  y3 <- as.matrix( dat3$pheno)
  y4 <- as.matrix( dat4$pheno)
  
  times <- dat1$sample_times
  geno_table <- dat1$snps
  nm <- dim(geno_table)[1]
  n1 <- interval[1]
  n2 <- interval[2]
  if(n2 >nm)
    n2 <- nm
  res <- matrix(NA,nrow=length(c(n1:n2)),ncol=120)
  for(i in n1:n2){
    SNP <- geno_table[i,]
    missing <- which(is.na(SNP))
    if ( length(missing) > 0)
    {
      SNP1 <- SNP[-(missing)]
      y11 <- y1[ -(missing), ]
      y22 <- y2[ -(missing), ]
      y33 <- y3[ -(missing), ]
      y44 <- y4[ -(missing), ]
    }else{
      SNP1 <- SNP
      y11 <- y1
      y22 <- y2
      y33 <- y3
      y44 <- y4
    }
    
    ndat1 <- dat1
    ndat2 <- dat2
    ndat3 <- dat3
    ndat4 <- dat4
    
    ndat1$pheno <- y11
    ndat2$pheno <- y22
    ndat3$pheno <- y33
    ndat4$pheno <- y44
    
    #parin2<-allpar1[9,-c(1:6)]
    h01 <- try(com.H0(ndat1,ndat2,ndat3,ndat4,parin2),TRUE)
    if (class(h01) == "try-error") 
      h01 <- NA
    h02 <- try(com.H1(y11,y22,y33,y44,SNP1,init.par=h01,times),TRUE)
    if (class(h02) == "try-error") 
      h02 <- NA
    LR <- 2*(h01[1]-h02[1])
    if(is.na(h01)||is.na(h02)){
      allpar <- c(LR,rep(NA,25))
    }else{
      allpar <- c(LR,h02)
    }
    

    cat("snp", i, "=", allpar[1:16], "\n");
    res[(i-(n1-1)),(1:length(allpar))] <- allpar
  }
  return(res)
}

com.H1 <- function(y11,y22,y33,y44,SNP1,init.par,times){
  
  index <- table(SNP1)
  snp.type <- as.numeric(names(index))
  
  g.par <- c()
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    SNP.n <- which(SNP1==snp.type[i])
    SNP.p <- as.numeric(c(colMeans(cbind(y11,y22,y33,y44)[SNP.n,],na.rm=T)))
    r1 <- optim(init.par[16:47],ind.mle,s.y=SNP.p,s.t=times,x1=SNP.p[1],x2=SNP.p[16],x3=SNP.p[31],x4=SNP.p[46],
                method="BFGS",control=list(maxit=32000))
    par <- r1$par
  
    g.par <- c(g.par,par)
    SNP.index[[i]] <- SNP.n
  }
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-5;
  max_err <- 1;
  
  #parin <- c(init.par[2:16],init.par[7:16])
  parin <- c(init.par[2:15],g.par)
  while(loop_k<max_iter && max_err>epsi){
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[-(1:14)])
      AA <- mle.fun(nnpar,y11=y11,y22=y22,y33=y33,y44=y44,times=times,SNP.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:14],mle.covar1,method = "BFGS",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- mle.fun(nnpar,y11=y11,y22=y22,y33=y33,y44=y44,times=times,SNP.index,snp.type)
      AA
    }
    r1.g <- optim(c(parin[-(1:14)]),mle1.g,method = "BFGS",control=list(maxit=32000))
    #cat("r1.g:",unlist(r1.g$par), "\n");
    
    newpar <- c(new.covar1,r1.g$par)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  return(c(r1.g$value,newpar))
}

mle.fun <- function(par,y11=y11,y22=y22,y33=y33,y44=y44,times=times,
                    SNP.index=SNP.index,snp.type=snp.type){
  
  Y1 <- cbind(y11,y22,y33,y44)
  n  <- length(y11[,1])
  
  len.cov <- 14
  par.covar <- par[1:len.cov]
  
  sigma <- SAD3.get_mat(par.covar,times, 4)
  len.gen <- 32
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    mu.g <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    yy1 <- Y1[SNP.index[[i]],]
    
    nyy1 <- yy1
    Myy1 <- as.numeric(colMeans(nyy1,na.rm=T))
    mu <- ind.get_mu(mu.g, times,x1=Myy1[1],x2=Myy1[16],x3=Myy1[31],x4=Myy1[46])

    fy1 <- dmvnorm( nyy1, mu, sigma)
    #fy1[which(fy1<=.Machine$double.eps)] <- .Machine$double.eps
    A1 <- c(A1,-sum(log(fy1)))
    len <- len + len.gen
  }
  A <- sum(A1)
  #cat("LL=",A,"\n")
  return (A);
}
