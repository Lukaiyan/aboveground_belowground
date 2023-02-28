com.parin <- function(dat,ret,file="sigqtl.csv",n=1){
  
  sigqtl <- read.csv(file)
  loci <- sigqtl[,1][n]
  loci.chr <- sigqtl[,2][n]
  
  parin <- ret[loci,1:80]
  parin <- c(parin,loci.chr,loci)
  
  parin
}


com.sim.par <- function(npar,n){
 
  par <- list(covar = as.numeric(npar[3:16]),
              c.par = as.numeric(npar[17:80]),
              times = seq(1,57,by=4),
              sample_p = 1,
              sample_n = n,
              chr=npar[81],
              loci=as.numeric(npar[82]))
  

  return(par)
}


com.simulate<-function( par,dat1,dat2,dat3,dat4 )
{
  sim_dat<-list(
    model        = "double",
    name         = "univariate",
    chr          = 1,
    sample_times = par$times,
    sample_p       = par$sample_p,
    pheno_file   = "simu.double.pheno",
    geno_file    = "simu.double.geno",
    raw_snps     = NULL,
    snps         = NULL,
    pheno_y      = NULL)
  
  class( sim_dat) <- "double.dat"
  
  #--generate snp data
  sim_dat$snps <- t(com_simu_geno(par))
  #-- generate SNP info
  snps.info <- c(1:par$sample_p)
  sim_dat$snps.info<- snps.info
  snps.names <- paste("snp", par$sample_p, sep="")
  sim_dat$snps.names<- snps.names
  
  #--generate phenotype data
  sim_dat$pheno_y <- com_simu_pheno(par,snp=sim_dat$snps,dat1,dat2,dat3,dat4);
  
  ##output to console
  cat("** Simulation is successful and data object is returned.\n");
  
  return(sim_dat);
}


com_simu_pheno<-function( par, snp,dat1,dat2,dat3,dat4)
{
  covar <- as.numeric(par$covar)
  
  SAD3 <- SAD3.get_mat(covar,par$times,4)

  mu <- list()
  y1 <- as.matrix(dat1$pheno)
  y2 <- as.matrix(dat2$pheno)
  y3 <- as.matrix(dat3$pheno)
  y4 <- as.matrix(dat4$pheno)
  
  geno_table <- dat1$snps
  
  SNP <- geno_table[par$loci,]
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
  index <- table(SNP1)
  snp.type <- as.numeric(names(index))
  
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    SNP.n <- which(SNP1==snp.type[i])
    SNP.p <- as.numeric(c(colMeans(cbind(y11,y22,y33,y44)[SNP.n,],na.rm=T)))
    x1<-SNP.p[1]
    x2<-SNP.p[16]
    x3<-SNP.p[31]
    x4<-SNP.p[46]
    mu[[i]] <- ind.get_mu(as.numeric(par$c.par[(32*(i-1)+1):(32*i)]),par$times,x1,x2,x3,x4)
    
  }

  simu_Y <- matrix(NA,par$sample_n,length(par$times)*4)
  for (i in 1:par$sample_n){
      ny <- rmvnorm(1,mu[[snp[i]+1]], SAD3 )
      while(any(ny<0)){
        ny <- rmvnorm(1,mu[[snp[i]+1]], SAD3)
      }
      simu_Y[i,]<-ny
    
  }

  ##output to console
  cat("** Phenotypical data is simulated successfully.\n")
  
  return(simu_Y);
}


com_simu_geno<-function( par ){
  c <- qnorm(1/4)
  
  all_SNPs <- c()
  for (i in 1:par$sample_n) 
  {
    all_SNPs <- rbind(all_SNPs, rnorm(par$sample_p, mean = 0, sd = 1))
  }	
  for (i in 1:par$sample_n) 
  {
    for (j in 1:par$sample_p) 
    {
      tmp <- all_SNPs[i,j];
      #all_SNPs[i,j] <- -1*(tmp<c)+(tmp>-c);#three genotype
      all_SNPs[i,j] <- -1*(tmp<c)#two genotype
    }
  }
  colnames(all_SNPs)<-paste("snp", c(1:par$sample_p),sep="");
  rownames(all_SNPs)<-paste("id", c(1:par$sample_n),sep="");
  
  ##output to console
  cat("** Genotypical data is simulated successfully.\n");
  
  #--QQ:2, Qq:1, qq:0
  return (all_SNPs+1);
}



delta<-function(dat1,dat2,dat3,dat4,par=parin,spar=allpar,loci=7490){
  m1<-matrix(NA,ncol=10,nrow=8);
  m2<-matrix(NA,ncol=10,nrow=8);
  m3<-matrix(NA,ncol=10,nrow=8);
  #m4<-matrix(NA,ncol=100,nrow=8);m5<-matrix(NA,ncol=100,nrow=8);
  for(i in 1:10){
    t<-seq(1,57,runif(1))
    #t<-sort(sample(seq(1,57,by=0.01),100,replace=TRUE))
    for(j in 1:4){
      result<-evalue(dat1,dat2,dat3,dat4,npar=as.numeric(par[17:80]),loci,t,ind=j)
      sim_result<-evalue1(npar=spar,t,ind=j)
      
        for(k in 1:2){
          phe<-result[[k]]$pheno
          phe1<-result[[k]]$pheno1
          phe2<-result[[k]]$pheno2
          #phe3<-result[[k]]$pheno3
          #phe4<-result[[k]]$pheno4
          
          sim_phe<-sim_result[[k]]$pheno
          sim_phe1<-sim_result[[k]]$pheno1
          sim_phe2<-sim_result[[k]]$pheno2
          #sim_phe3<-sim_result[[k]]$pheno3
          #sim_phe4<-sim_result[[k]]$pheno4
          
          m1[2*(j-1)+k,i]<-max(sqrt((phe-sim_phe)^2))
          m2[2*(j-1)+k,i]<-max(sqrt((phe1-sim_phe1)^2))
          m3[2*(j-1)+k,i]<-max(sqrt((phe2-sim_phe2)^2))
          #m4[2*(j-1)+k,i]<-max(sqrt((phe3-sim_phe3)^2))
          #m5[2*(j-1)+k,i]<-max(sqrt((phe4-sim_phe4)^2))
        }
      
      }
  }
  
  mm<-c()
  for(k in 1:8){
    mm<-c(mm,c(min(m1[k,],na.rm = T),min(m2[k,],na.rm = T),min(m3[k,],na.rm = T)))
  }
  return(mm)
}



evalue<-function(dat1,dat2,dat3,dat4,npar,loci,t,ind){
  
  y1 <- as.matrix(dat1$pheno)
  y2 <- as.matrix(dat2$pheno)
  y3 <- as.matrix(dat3$pheno)
  y4 <- as.matrix(dat4$pheno)
  geno_table <- dat1$snps
  
  SNP <- geno_table[loci,]
  
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
  
  index <- table(SNP)
  snp.type <- as.numeric(names(index))
  
  result<-list()
  for(i in 1:length(snp.type)){
    SNP.n <- which(SNP1==snp.type[i])
    SNP.p <- as.numeric(c(colMeans(cbind(y11,y22,y33,y44)[SNP.n,],na.rm=T)))
   
    c.par<-npar[(32*(i-1)+1):(32*i)]
    fity<-ind.get_mu(par=c.par,times=t,
                      x1=SNP.p[1],x2=SNP.p[16],x3=SNP.p[31],x4=SNP.p[46])
    fity1<-ind.get_mu1(par=c.par,times=t,
                      x1=SNP.p[1],x2=SNP.p[16],x3=SNP.p[31],x4=SNP.p[46])
    fit.tn <- data.frame(pheno=fity[((ind-1)*length(t)+1):(ind*length(t))],
                         pheno1=fity1[((ind-1)*length(t)+1):(ind*length(t))],
                         pheno2=fity[((ind-1)*length(t)+1):(ind*length(t))]-fity1[((ind-1)*length(t)+1):(ind*length(t))])
    result[[i]]<-fit.tn
  }
  
  return(result)
}



evalue1<-function(npar,t,ind){

  result<-list()
  for(i in 1:2){
    fity_all<-c()
    fity_all1<-c()
    for(iter in 1:dim(npar)[1]){
      x1=as.numeric(npar[iter,64+(i-1)*4+1])
      x2=as.numeric(npar[iter,64+(i-1)*4+2])
      x3=as.numeric(npar[iter,64+(i-1)*4+3])
      x4=as.numeric(npar[iter,64+(i-1)*4+4])
      c.par<-npar[iter,(32*(i-1)+1):(32*i)]
      tmpfity<-ind.get_mu(par=c.par,times=t,x1,x2,x3,x4)
      tmpfity1<-ind.get_mu1(par=c.par,times=t,x1,x2,x3,x4)
      if((!any(is.nan(c(tmpfity,tmpfity1))))&&any(c(tmpfity,tmpfity1)!=-Inf)){
        fity_all<-rbind(fity_all,tmpfity)
        fity_all1<-rbind(fity_all1,tmpfity1)
      }
    }
    fity<-colMeans(fity_all)
    fity1<-colMeans(fity_all1)
    fit.tn <- data.frame(pheno=fity[((ind-1)*length(t)+1):(ind*length(t))],
                         pheno1=fity1[((ind-1)*length(t)+1):(ind*length(t))],
                         pheno2=fity[((ind-1)*length(t)+1):(ind*length(t))]-fity1[((ind-1)*length(t)+1):(ind*length(t))])
    result[[i]]<-fit.tn
  }
  
  return(result)
}




sim.test<-function(dat1,dat2,dat3,dat4,parin,h2,n,max_iter){
  
  allpar<-c()
  for(iter in 1:max_iter){
    par_obj<- com.sim.par(npar=parin,n)
    var.g<-c(0.01,0.02,0.01,0.01)
    #var.g<-c(0.1,0.2,0.1,0.1)
    var.e <- (1/h2-1)*var.g
    par_obj1<-par_obj
    par_obj1$covar[c(2,4,6,8)] <- sqrt(var.e)
    
    dat_obj <- com.simulate(par=par_obj1,dat1,dat2,dat3,dat4)
    dat_obj1<-dat_obj
    dat_obj1$pheno<-dat_obj$pheno_y[,c(1:15)]
    dat_obj2<-dat_obj
    dat_obj2$pheno<-dat_obj$pheno_y[,c(16:30)]
    dat_obj3<-dat_obj
    dat_obj3$pheno<-dat_obj$pheno_y[,c(31:45)]
    dat_obj4<-dat_obj
    dat_obj4$pheno<-dat_obj$pheno_y[,c(46:60)]
    
    init_pars<-c(49,0.1,102,0.06,33,0.13,19,0.11,
                 0.007,-0.059,0.012,0.372,-0.035,0.094,
                 0.801,0.346,-0.588,0.296,-0.162,0.221, 
                 0.143,-0.33,-0.077,-0.016,-0.072,-0.43,
                 -0.03,0.179, 0.143,0.056,-0.099,-0.09)
    fit_result1<-try(get_value(dat_obj1,dat_obj2,dat_obj3,dat_obj4,init_pars),TRUE)
      for(j in 1:8){
        for(k in 1:8){
          for(l in 1:8){
            parin2<-c(0.2,sd(colMeans(dat_obj1$pheno)),0.1+j*0.1,sd(colMeans(dat_obj2$pheno)),0.1+k*0.1,sd(colMeans(dat_obj3$pheno)),0.1+l*0.1,sd(colMeans(dat_obj4$pheno)),
                      0.5,0.5,0.5,0.5,0.5,0.5,fit_result1$par)
            #           0.5,0.5,0.5,0.5,0.5,0.5,parin2[-c(1:14)])
            ret_obj<-try(com.est1(dat_obj1,dat_obj2,dat_obj3,dat_obj4,interval=c(1,1),parin2),TRUE)
            
            if(!is.na(ret_obj[1])) break
          }
          if(!is.na(ret_obj[1])) break
        }
        if(!is.na(ret_obj[1])) break
      }
    
    if(is.na(ret_obj[1])||any(ret_obj[1:16]<0)){
      ret_obj<-c()
    }else{
      y1 <- as.matrix(dat_obj1$pheno)
      y2 <- as.matrix(dat_obj2$pheno)
      y3 <- as.matrix(dat_obj3$pheno)
      y4 <- as.matrix(dat_obj4$pheno)
      SNP <- dat_obj$snps
      
      index <- table(SNP)
      snp.type <- as.numeric(names(index))
      
      xi<-c()
      for(i in 1:length(snp.type)){
        SNP.n <- which(SNP==snp.type[i])
        SNP.p <- as.numeric(c(colMeans(cbind(y1,y2,y3,y4)[SNP.n,],na.rm=T)))
        
        xi<- c(xi,c(SNP.p[1],SNP.p[16],SNP.p[31],SNP.p[46]))
      }
      
      ret_obj<-c(ret_obj[17:80],xi)
      allpar<-rbind(allpar,ret_obj)
    }
  }
  
  return(allpar)
}




sim.test<-function(dat1,dat2,dat3,dat4,parin,h2,n,max_iter){
  
  allpar<-c()
  for(iter in 1:max_iter){
    par_obj<- com.sim.par(npar=parin,n)
    var.g<-c(0.05,0.05,0.05,0.05)
    var.e <- (1/h2-1)*var.g
    par_obj1<-par_obj
    par_obj1$covar[c(2,4,6,8)] <- sqrt(var.e)
    
    dat_obj <- com.simulate(par=par_obj1,dat1,dat2,dat3,dat4)
    dat_obj1<-dat_obj
    dat_obj1$pheno<-dat_obj$pheno_y[,c(1:15)]
    dat_obj2<-dat_obj
    dat_obj2$pheno<-dat_obj$pheno_y[,c(16:30)]
    dat_obj3<-dat_obj
    dat_obj3$pheno<-dat_obj$pheno_y[,c(31:45)]
    dat_obj4<-dat_obj
    dat_obj4$pheno<-dat_obj$pheno_y[,c(46:60)]
    
    
    init_pars<-c(49,0.1,102,0.06,33,0.13,19,0.11,
                 0.007,-0.059,0.012,0.372,-0.035,0.094,
                 0.801,0.346,-0.588,0.296,-0.162,0.221, 
                 0.143,-0.33,-0.077,-0.016,-0.072,-0.43,
                 -0.03,0.179, 0.143,0.056,-0.099,-0.09)
    fit_result1<-get_value(dat_obj1,dat_obj2,dat_obj3,dat_obj4,init_pars)
    parin2<-c(0.5,sd(colMeans(dat_obj1$pheno)),0.5,sd(colMeans(dat_obj2$pheno)),0.5,sd(colMeans(dat_obj3$pheno)),0.5,sd(colMeans(dat_obj4$pheno)),
              0.5,0.5,0.5,0.5,0.5,0.5,fit_result1$par)
    ret_obj<-try(com.est1(dat_obj1,dat_obj2,dat_obj3,dat_obj4,interval=c(1,1),parin2),TRUE)
    
    if(is.na(ret_obj[1])||any(ret_obj[1:16]<0)){
      ret_obj<-c()
    }else{
      y1 <- as.matrix(dat_obj1$pheno)
      y2 <- as.matrix(dat_obj2$pheno)
      y3 <- as.matrix(dat_obj3$pheno)
      y4 <- as.matrix(dat_obj4$pheno)
      SNP <- dat_obj$snps
      
      index <- table(SNP)
      snp.type <- as.numeric(names(index))
      
      xi<-c()
      for(i in 1:length(snp.type)){
        SNP.n <- which(SNP==snp.type[i])
        SNP.p <- as.numeric(c(colMeans(cbind(y1,y2,y3,y4)[SNP.n,],na.rm=T)))
        
        xi<- c(xi,c(SNP.p[1],SNP.p[16],SNP.p[31],SNP.p[46]))
      }
      
      ret_obj<-c(ret_obj[17:80],xi)
      allpar<-rbind(allpar,ret_obj)
    }
    }
    
  return(allpar)
}


sim.test<-function(dat1,dat2,dat3,dat4,parin,h2,n,max_iter){
  
  allpar<-c()
  for(iter in 1:max_iter){
    par_obj<- com.sim.par(npar=parin,n)
    var.g<-c(0.01,0.02,0.01,0.01)
    var.e <- (1/h2-1)*var.g
    par_obj1<-par_obj
    par_obj1$covar[c(2,4,6,8)] <- sqrt(var.e)
    
    dat_obj <- com.simulate(par=par_obj1,dat1,dat2,dat3,dat4)
    dat_obj1<-dat_obj
    dat_obj1$pheno<-dat_obj$pheno_y[,c(1:15)]
    dat_obj2<-dat_obj
    dat_obj2$pheno<-dat_obj$pheno_y[,c(16:30)]
    dat_obj3<-dat_obj
    dat_obj3$pheno<-dat_obj$pheno_y[,c(31:45)]
    dat_obj4<-dat_obj
    dat_obj4$pheno<-dat_obj$pheno_y[,c(46:60)]
    
    for(j in 1:8){
      for(k in 1:8){
        for(l in 1:8){
          parin2<-c(0.2,sd(colMeans(dat_obj1$pheno)),0.1+j*0.1,sd(colMeans(dat_obj2$pheno)),0.1+k*0.1,sd(colMeans(dat_obj3$pheno)),0.1+l*0.1,sd(colMeans(dat_obj4$pheno)),
          #          0.5,0.5,0.5,0.5,0.5,0.5,fit_result1$par)
                     0.5,0.5,0.5,0.5,0.5,0.5,parin2[-c(1:14)])
          
          ret_obj<-try(com.est1(dat1=dat_obj1,dat2=dat_obj2,dat3=dat_obj3,dat4=dat_obj4,interval=c(1,1),parin2),TRUE)[17:80]
          
          if(!is.na(ret_obj[1])) break
        }
        if(!is.na(ret_obj[1])) break
      }
      if(!is.na(ret_obj[1])) break
    }
    if(is.na(ret_obj[1])||any(ret_obj[1:16]<0)){
      ret_obj<-c()
    }else{
      y1 <- as.matrix(dat_obj1$pheno)
      y2 <- as.matrix(dat_obj2$pheno)
      y3 <- as.matrix(dat_obj3$pheno)
      y4 <- as.matrix(dat_obj4$pheno)
      SNP <- dat_obj$snps
      
      index <- table(SNP)
      snp.type <- as.numeric(names(index))
      
      xi<-c()
      for(i in 1:length(snp.type)){
        SNP.n <- which(SNP==snp.type[i])
        SNP.p <- as.numeric(c(colMeans(cbind(y1,y2,y3,y4)[SNP.n,],na.rm=T)))
        
        xi<- c(xi,c(SNP.p[1],SNP.p[16],SNP.p[31],SNP.p[46]))
      }
      
      ret_obj<-c(ret_obj,xi)
      allpar<-rbind(allpar,ret_obj)
    }
  }
  
  return(allpar)
}


sim.test<-function(dat1,dat2,dat3,dat4,parin,h2,n,max_iter){
  
  allpar<-c()
  for(iter in 1:max_iter){
    par_obj<- com.sim.par(npar=parin,n)
    var.g<-c(0.01,0.02,0.01,0.01)
    var.e <- (1/h2-1)*var.g
    par_obj1<-par_obj
    par_obj1$covar[c(2,4,6,8)] <- sqrt(var.e)
    
    dat_obj <- com.simulate(par=par_obj1,dat1,dat2,dat3,dat4)
    dat_obj1<-dat_obj
    dat_obj1$pheno<-dat_obj$pheno_y[,c(1:15)]
    dat_obj2<-dat_obj
    dat_obj2$pheno<-dat_obj$pheno_y[,c(16:30)]
    dat_obj3<-dat_obj
    dat_obj3$pheno<-dat_obj$pheno_y[,c(31:45)]
    dat_obj4<-dat_obj
    dat_obj4$pheno<-dat_obj$pheno_y[,c(46:60)]
    
    for(j in 1:8){
      for(k in 1:8){
        for(l in 1:8){
          parin2<-c(0.2,sd(colMeans(dat_obj1$pheno)),0.1+j*0.1,sd(colMeans(dat_obj2$pheno)),0.1+k*0.1,sd(colMeans(dat_obj3$pheno)),0.1+l*0.1,sd(colMeans(dat_obj4$pheno)),
                    #          0.5,0.5,0.5,0.5,0.5,0.5,fit_result1$par)
                    0.5,0.5,0.5,0.5,0.5,0.5,parin2[-c(1:14)])
          
          ret_obj<-try(com.est1(dat1=dat_obj1,dat2=dat_obj2,dat3=dat_obj3,dat4=dat_obj4,interval=c(1,1),parin2),TRUE)[17:80]
          
          if(!is.na(ret_obj[1])) break
        }
        if(!is.na(ret_obj[1])) break
      }
      if(!is.na(ret_obj[1])) break
    }
    if(is.na(ret_obj[1])){
      ret_obj<-c()
    }else{
      y1 <- as.matrix(dat_obj1$pheno)
      y2 <- as.matrix(dat_obj2$pheno)
      y3 <- as.matrix(dat_obj3$pheno)
      y4 <- as.matrix(dat_obj4$pheno)
      SNP <- dat_obj$snps
      
      index <- table(SNP)
      snp.type <- as.numeric(names(index))
      
      xi<-c()
      for(i in 1:length(snp.type)){
        SNP.n <- which(SNP==snp.type[i])
        SNP.p <- as.numeric(c(colMeans(cbind(y1,y2,y3,y4)[SNP.n,],na.rm=T)))
        
        xi<- c(xi,c(SNP.p[1],SNP.p[16],SNP.p[31],SNP.p[46]))
      }
      
      ret_obj<-c(ret_obj,xi)
      allpar<-rbind(allpar,ret_obj)
    }
  }
  
  return(allpar)
}

