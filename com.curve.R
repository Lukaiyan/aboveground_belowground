com.H0 <- function(dat1,dat2,dat3,dat4,parin2){
  
  allpheno <- cbind(dat1$pheno,dat2$pheno,dat3$pheno,dat4$pheno)
  times <- dat1$sample_times

  #parin2<-c(0.5,3.0366,0.5,17,0.5,6.1777,0.5,7.7994,0.5,0.5,0.5,0.5,0.4,0.5,
  #          fit_result$par)
  
  mpheno <- as.numeric(colMeans(allpheno))
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin2);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin2[15:46])
      AA <- curve.mle(nnpar,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[16],x3=mpheno[31],x4=mpheno[46])
      AA
    }
    r1.covar <- optim(parin2[1:14],mle.covar1,method = "BFGS",control=list(maxit=2000))
    #r1.covar <- optim(parin2[1:14],mle.covar1,method = "Nelder-Mead",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[16],x3=mpheno[31],x4=mpheno[46])
      AA
    }
    r1 <- optim(c(parin2[15:46]),mle.1,method = "BFGS",control=list(maxit=32000)) 
    #r1 <- optim(c(parin2[16:47]),mle.1,method = "Nelder-Mead",control=list(maxit=2000)) 
    new1 <- r1$par
    
    #r2<- optim(c(new.covar1,new1),curve.mle,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[2],
               #method = "BFGS",control=list(maxit=32000,trace=T))    
    #cat("new1:",unlist( new1), "\n");
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin2,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[16],x3=mpheno[31],x4=mpheno[46])
  return(c(LL,parin2))
}



curve.mle <-function( par,y,time.std,x1,x2,x3,x4)
{
  len.cov <- 14
  par.covar <- par[1:len.cov]
  n  <- length(y[,1])
  #sig.inv3 <- SAD3.get_inv_mat(par.covar,time.std, 2)
  sigma <- SAD3.get_mat(par.covar,time.std, 4)#solve( sig.inv3 )
  
  curve.par <- par[(len.cov+1):(len.cov+ 32)]
  mu <- ind.get_mu(curve.par,time.std,x1,x2,x3,x4)
  
  yy <- y
  fy <- dmvnorm(yy,mu,sigma)
  #fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps
  A <- -sum(log(fy))
  #cat("LL=",A,"\n")
  return(A)
}