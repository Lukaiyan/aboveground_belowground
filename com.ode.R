
require("deSolve")

s.mle <- function(s.par,s.y,s.t){
  A <- sum((s.y - com.get_mu(s.par,s.t))^2 )
  A
}

ind.mle <- function(s.par,s.y,s.t,x1,x2,x3,x4){
  A <- sum((s.y - ind.get_mu(s.par,s.t,x1,x2,x3,x4))^2 )
  # cat("LL=",A,"\n")
  A
}

com.get_mu <- function(par, times, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8],
      l11= par[9],
      l12= par[10],
      l13= par[11],
      l14= par[12],
      l15= par[13],
      l16= par[14],
      l21= par[15],
      l22= par[16],
      l23= par[17],
      l24= par[18],
      l25= par[19],
      l26= par[20],
      l31= par[21],
      l32= par[22],
      l33= par[23],
      l34= par[24],
      l35= par[25],
      l36= par[26],
      l41= par[27],
      l42= par[28],
      l43= par[29],
      l44= par[30],
      l45= par[31],
      l46= par[32]);
  }
  
  state0 <- c(X1=0.4437582,X2=3.734670,X3=0.1264504,X4=0.37056179);
  y <- COMP.f( par0, state0, times );
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}


COMP.f <-function( parameters, state, times,max_state ){
  


  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX1 <- r1*X1*(1-(X1/K1))+l11*X2^l12+l13*X3^l14+l15*X4^l16
            dX2 <- r2*X2*(1-(X2/K2))+l21*X1^l22+l23*X3^l24+l25*X4^l26
            dX3 <- r3*X3*(1-(X3/K3))+l31*X1^l32+l33*X2^l34+l35*X4^l36
            dX4 <- r4*X4*(1-(X4/K4))+l41*X1^l42+l43*X2^l44+l45*X3^l46
            list(c(dX1, dX2, dX3 ,dX4))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

com.get_mu1 <- function(par, times, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8]);
  }
  
  state0 <- c(X1=0.4437582,X2=3.734670,X3=0.1264504,X4=0.37056179);
  y <- COMP.f1( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}


COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX1 <- r1*X1*(1-(X1/K1))
            dX2 <- r2*X2*(1-(X2/K2))
            dX3 <- r3*X3*(1-(X3/K3))
            dX4 <- r4*X4*(1-(X4/K4))
            
            list(c(dX1, dX2, dX3 ,dX4))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}



com.get_mu2 <- function(par, times, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8],
      l11= par[9],
      l12= par[10],
      l21= par[15],
      l22= par[16],
      l31= par[21],
      l32= par[22],
      l41= par[27],
      l42= par[28]);
  }
  
  state0 <- c(X1=0.4437582,X2=3.734670,X3=0.1264504,X4=0.37056179);
  y <- COMP.f2( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}


COMP.f2 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX1 <- r1*X1*(1-(X1/K1))+l11*X2^l12
            dX2 <- r2*X2*(1-(X2/K2))+l21*X1^l22
            dX3 <- r3*X3*(1-(X3/K3))+l31*X1^l32
            dX4 <- r4*X4*(1-(X4/K4))+l41*X1^l42
            
            list(c(dX1, dX2, dX3 ,dX4))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


com.get_mu3 <- function(par, times, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8],
      l11= par[9],
      l12= par[10],
      l13= par[11],
      l14= par[12],
      l21= par[15],
      l22= par[16],
      l23= par[17],
      l24= par[18],
      l31= par[21],
      l32= par[22],
      l33= par[23],
      l34= par[24],
      l41= par[27],
      l42= par[28],
      l43= par[29],
      l44= par[30]);
  }
  
  state0 <- c(X1=0.4437582,X2=3.734670,X3=0.1264504,X4=0.37056179);
  y <- COMP.f3( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}


COMP.f3 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX1 <- r1*X1*(1-(X1/K1))+l11*X2^l12+l13*X3^l14
            dX2 <- r2*X2*(1-(X2/K2))+l21*X1^l22+l23*X3^l24
            dX3 <- r3*X3*(1-(X3/K3))+l31*X1^l32+l33*X2^l34
            dX4 <- r4*X4*(1-(X4/K4))+l41*X1^l42+l43*X2^l44
            list(c(dX1, dX2, dX3 ,dX4))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}



ind.get_mu <- function(par, times, x1,x2,x3,x4,options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8],
      l11= par[9],
      l12= par[10],
      l13= par[11],
      l14= par[12],
      l15= par[13],
      l16= par[14],
      l21= par[15],
      l22= par[16],
      l23= par[17],
      l24= par[18],
      l25= par[19],
      l26= par[20],
      l31= par[21],
      l32= par[22],
      l33= par[23],
      l34= par[24],
      l35= par[25],
      l36= par[26],
      l41= par[27],
      l42= par[28],
      l43= par[29],
      l44= par[30],
      l45= par[31],
      l46= par[32]);
  }
  
  state0 <- c(X1=x1,X2=x2,X3=x3,X4=x4);
  y <- COMP.f( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}




ind.get_mu1 <- function(par, times, x1,x2,x3,x4,options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8]);
  }
  
  state0 <- c(X1=x1,X2=x2,X3=x3,X4=x4);
  y <- COMP.f1( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}



ind.get_mu2 <- function(par, times,x1,x2,x3,x4)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8],
      l11= par[9],
      l12= par[10],
      l21= par[15],
      l22= par[16],
      l31= par[21],
      l32= par[22],
      l41= par[27],
      l42= par[28]);
  }
  
  state0 <- c(X1=x1,X2=x2,X3=x3,X4=x4);
  y <- COMP.f2( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}


ind.get_mu3 <- function(par, times, x1,x2,x3,x4,options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[3],
      r2 = par[4],
      K3 = par[5],
      r3 = par[6],
      K4 = par[7],
      r4 = par[8],
      l11= par[9],
      l12= par[10],
      l13= par[11],
      l14= par[12],
      l21= par[15],
      l22= par[16],
      l23= par[17],
      l24= par[18],
      l31= par[21],
      l32= par[22],
      l33= par[23],
      l34= par[24],
      l35= par[25],
      l36= par[26],
      l41= par[27],
      l42= par[28],
      l43= par[29],
      l44= par[30]);
  }
  
  state0 <- c(X1=x1,X2=x2,X3=x3,X4=x4);
  y <- COMP.f3( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:5] ) );
}
