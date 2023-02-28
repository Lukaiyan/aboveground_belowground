library(deSolve)
library(ggplot2)
library(mvtnorm)
library(tidyr)
library(xlsx)
library(reshape)
library(ggrepel)
library(patchwork)


load("./stem.Rdata")#dat1
load("./taproot.Rdata")#dat2
load("./lat_l.Rdata")#dat3
load("./lat_meanl.Rdata")#dat4


source("./com.ode.R")
source("./com.covar.R")
source("./com.curve.R")
source("./com.optim.R")



ret1<-com.est1(dat1,dat2,dat3,dat4,interval=c(7525,7525))



############################################cluster########################################
source("./fun.cluster.R")

sd_value<-SD(dat1,dat2,dat3,dat4,interval=c(1,8305))
times<-seq(1,57,4)
initial_par<-get_initial_par(effect=sd_value,k=130,r=5)
clus<-get_cluster(data=sd_value,k=130,r=5,init_sd_par=initial_par$init_sd_par,
            init_curve_par=initial_par$init_curve_par,init_pro=initial_par$init_pro)



source("network.all.fun.R")
all_para1 <- t(apply(sd_value1[,c(1:15)], 1,smooth.optim_ind,times=t,
                    para=rep(0.01,6)))
rownames(all_para1) <- c(1:dim(sd_value1)[1])

smooth_data1 <- t(apply(all_para1, 1, Legendre.model,t=seq(1,57,0.5)))
rownames(smooth_data1) <- c(1:dim(sd_value1)[1])

library(orthogonalsplinebasis)
Net<-optim_inter(all_cluster_value=smooth_data1,ind_Lpar11=all_para1,norder=5)
result1<-out_Netdata_cytoscape(optim_cluster=Net,which_cluster="cluster-1",clusterAA=smooth_data1)



######################################simulation###########################################
library(deSolve)
library(mvtnorm)
source("./com.sim.R")
parin <- c(ret[7539,c(1:80)],"lg1",7539)
allpar2<-sim.test(dat1,dat2,dat3,dat4,parin,h2=0.05,n=345,max_iter=10)
allpar4<-sim.test(dat1,dat2,dat3,dat4,parin,h2=0.1,n=345,max_iter=20)
allpar5<-sim.test(dat1,dat2,dat3,dat4,parin,h2=0.05,n=100,max_iter=10)
allpar6<-sim.test(dat1,dat2,dat3,dat4,parin,h2=0.1,n=100,max_iter=10)

