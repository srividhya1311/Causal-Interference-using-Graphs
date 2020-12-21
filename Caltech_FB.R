#
#
#
#******************************************************************************************************
#*************************Data for the Caltech Facebook network****************************************
#******************************************************************************************************
#
#
#
#.mat files from the Facebook100 (FB100) dataset should be placed in the data folder. 
#The Facebook100 (FB100) dataset is publicly available from the Internet Archive 
#at https://archive.org/details/oxford-2005-facebook-matrix and other public repositories.
cal_all=read.mat('Caltech36.mat')
str(cal_all)
adj=cal_all$A
g=graph_from_adjacency_matrix(adj)
params = purrr::cross(list(
  int = 0,
  direct = 1,
  spill = c(0.1, 0.25, 0.5, 0.75, 1),
  max_t = c(2),
  noise_sd = 1
))
#
#
#
#******************************************************************************************************
#*************************Sampling and estimation in Caltech FB network********************************
#******************************************************************************************************
#
#
#
#estimates in the complete network
mc<-5000
la_g<-replicate(mc,get_estimates(g,covariate_fns,params)[[1]])
dm_g<-replicate(mc,get_estimates(g,covariate_fns,params)[[2]])
ha_g<-replicate(mc,get_estimates(g,covariate_fns,params)[[3]])
#estimates in the sampled subgraphs
mc<-5
la_sg_final<-array(NA,dim=c(20,5,mc))
dm_sg_final<-array(NA,dim=c(20,5,mc))
ha_sg_final<-array(NA,dim=c(20,5,mc))
for (j in 1:mc){
  sample_graphs<-get_sg(g,FUN = random_edge)
  n_samples<-length(sample_graphs)
  la_sg<-matrix(NA,nrow=n_samples,ncol=5)
  dm_sg<-matrix(NA,nrow=n_samples,ncol=5)
  ha_sg<-matrix(NA,nrow=n_samples,ncol=5)
  for(i in 1:n_samples){
    sim=50
    sg<-sample_graphs[[i]]
    la_sg[i,]<-apply(replicate(sim,get_estimates(sg,covariate_fns,params)[[1]]),1,mean)
    dm_sg[i,]<-apply(replicate(sim,get_estimates(sg,covariate_fns,params)[[2]]),1,mean)
    ha_sg[i,]<-apply(replicate(sim,get_estimates(sg,covariate_fns,params)[[3]]),1,mean)
  }
  la_sg_final[,,j]<-la_sg
  dm_sg_final[,,j]<-dm_sg
  ha_sg_final[,,j]<-ha_sg
}
#column 1 with spillover = 0.1
final_la<-cbind(la_sg_final[,1,])
final_dm<-cbind(dm_sg_final[,1,])
final_ha<-cbind(ha_sg_final[,1,])
#plotting the results
p<-seq(0.25,0.99,length.out = 20)
par(mfrow=c(1,1),bg="grey100",xpd=F,cex=1.5)
plot(p,final_la[,1],type='l',col='white',
     ylim=c(0.31,0.36),lwd=2,lty=2,xlab = "Proportion of edges sampled from the complete graph",
     ylab = "GATE Estimates")
smoother<-loess.smooth(p,apply(final_la,1,mean),span=2/3)
lines(smoother$x,smoother$y,col="red")
abline(h=apply(la_g,1,mean)[1],col="red",lwd=2,lty=2)
abline(h=apply(dm_g,1,mean)[1],col="green",lwd=2,lty=2)
smoother<-loess.smooth(p,apply(final_dm,1,mean),span=1/3)
lines(smoother$x,smoother$y,col="green")
grid(20,7,col="grey75")