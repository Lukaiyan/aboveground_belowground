SD<-function(dat1,dat2,dat3,dat4,interval=c(1,10)){
  
  y1 <- as.matrix( dat1$pheno)
  y2 <- as.matrix( dat2$pheno)
  y3 <- as.matrix( dat3$pheno)
  y4 <- as.matrix( dat4$pheno)
  
  times <- dat1$sample_times
  geno_table <- dat1$snps
  nm <- dim(geno_table)[1]
  n1 <- interval[1]
  n2 <- interval[2]
  
  SD<-c()
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
    
    index <- table(SNP1)
    snp.type <- as.numeric(names(index))
    
  
  SNP.phe <- c()
  for(i in 1:length(snp.type)){
    SNP.n <- which(SNP1==snp.type[i])
    SNP.p <- as.numeric(c(colMeans(cbind(y11,y22,y33,y44)[SNP.n,],na.rm=T)))
    
    SNP.phe<-rbind(SNP.phe,SNP.p)
  }
  SD<-rbind(SD,apply(SNP.phe,2,sd))
  }
  SD[which(is.na(SD))]<-0
  return(SD)
}



LgdP <- expression( 1,tt,
                    ( 3* tt^2 - 1 )/2 , 
                    ( 5 * tt^3 - 3* tt )/2, 
                    ( 35 * tt^4 - 30 * tt^2 + 3)/8,
                    ( 63 * tt^5 - 70 * tt^3 + 15 * tt )/8,
                    ( 231 * tt^6 - 315 * tt^4 + 105 * tt^2 - 5)/16  )

Legendre<-function(par,r,times){
  tnum = length(times) 
  f<-c()
  for(t in 1:tnum ){
    tt <- -1 + 2*(times[t] - times[1])/(times[tnum] - times[1])
    ff<-0
    for(i in 1:r){
      ff<-ff+par[i]*eval(LgdP[i])
    }
    f<-c(f,ff)
  }
  return(f)
}

s.mle<-function(par,data,r,times){
  y <- Legendre(par,r,times)
  yi <- data
  res <- sum((yi-y)^2)
  return(res)
}

get_Legendre_par<-function(initial_f_par,data,r,times){
  
  a <- optim(initial_f_par,s.mle,data=data,r=r,times=times,method = "Nelder-Mead")
  curve_par_i<-a$par
}

get_initial_par <- function(effect,k,r){
  
  init_cluster <- kmeans(effect,k,iter.max = 100)
  init_curve_para<-c()
  init_sd_para<-c()
  for(c in 1:k){
    tmp<-effect[which(init_cluster$cluster==c),]
    initial_f_par <- rep(0.1,r)
    a1<-get_Legendre_par(initial_f_par,colMeans(tmp)[1:length(times)],r,times)
    a2<-get_Legendre_par(initial_f_par,colMeans(tmp)[(length(times)+1):(2*length(times))],r,times)
    a3<-get_Legendre_par(initial_f_par,colMeans(tmp)[(2*length(times)+1):(3*length(times))],r,times)
    a4<-get_Legendre_par(initial_f_par,colMeans(tmp)[(3*length(times)+1):(4*length(times))],r,times)
    init_curve_para<-c(init_curve_para,c(a1,a2,a3,a4))
    
    cusd.s <- mean(apply(tmp[,1:length(times)],2,sd))
    cusd.t <- mean(apply(tmp[,(length(times)+1):(2*length(times))],2,sd))
    cusd.l <- mean(apply(tmp[,(2*length(times)+1):(3*length(times))],2,sd))
    cusd.lm <- mean(apply(tmp[,(3*length(times)+1):(4*length(times))],2,sd))
    init_sd_para<-c(init_sd_para,c(0.5,cusd.s,0.5,cusd.t,0.5,cusd.l,0.5,cusd.lm,0.5,0.5,0.5,0.5,0.5,0.5))
  }  
  
  init_pro <- table(init_cluster$cluster)/nrow(effect)
  
  return_object <- list(init_sd_para,init_curve_para,init_pro)
  names(return_object)<-c("init_sd_par","init_curve_par","init_pro")
  return(return_object)
}

get_cluster <- function(data=effect,k,r,init_sd_par,init_curve_par,init_pro){
  Delta <- 1000; iter <- 0; itermax <- 1000;
  
  mle <- function(par,data,prob){
    par1<-par[1:(14*k)]
    par2<-par[-c(1:(14*k))]
    tmp_S<- 0
    for (c in 1:k) {
      covM<-SAD3.get_mat(par1[(1+14*(c-1)):(14*c)],times, 4)
      cpar<-par2[(1+4*r*(c-1)):(4*r*c)]
      mcurve<-c(Legendre(cpar[1:r],r,times),Legendre(cpar[(r+1):(2*r)],r,times),
                Legendre(cpar[(2*r+1):(3*r)],r,times),Legendre(cpar[(3*r+1):(4*r)],r,times))
      tmp_S<-tmp_S+dmvnorm(data,mcurve,covM)*prob[c]
    }
    
    LL <- -sum(log(tmp_S))
    return(LL)
  }
  
  while ( Delta > 0.1 && iter <= itermax ) {
    # initiation
    if(iter == 0){
      init_sd_para <- init_sd_par
      init_curve_para <- init_curve_par 
      pro <- init_pro
      
    }
    #E step, calculate the posterior probability
    old_par <- c(init_sd_para,init_curve_para)
    LL_mem <- mle(old_par,data,pro)
    
    tmp_M<-matrix(NA,nrow = dim(data)[1],ncol = k)
    for (c in 1:k) {
      covM<-SAD3.get_mat(init_sd_para[(1+14*(c-1)):(14*c)],times, 4)
      cpar<-init_curve_para[(1+4*r*(c-1)):(4*r*c)]
      mcurve<-c(Legendre(cpar[1:r],r,times),Legendre(cpar[(r+1):(2*r)],r,times),
                Legendre(cpar[(2*r+1):(3*r)],r,times),Legendre(cpar[(3*r+1):(4*r)],r,times))
      tmp_M[,c]<-dmvnorm(data,mcurve,covM)*pro[c]
    }
    
    omega <- tmp_M/rowSums(tmp_M)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    new_par <- optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead")
    L_Value <- new_par$value
    
    init_sd_para <- new_par$par[1:(14*k)]
    init_curve_para <- new_par$par[-c(1:(14*k))]
    Delta <- abs(L_Value-LL_mem)
    cat('\n',"iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  BIC <- 2*L_Value+log(nrow(data))*length(new_par$par)
  cluster <- apply(omega,1,which.max)
  clustered_df <- cbind(row.names(data),data,cluster)
  
  return_object <- list(init_sd_para,init_curve_para,pro,LL_mem,BIC,clustered_df)
  names(return_object)<-c("sd_par", "curve_par", "pro", "LL", "BIC","clustered_data")
  
  return(return_object)
}
