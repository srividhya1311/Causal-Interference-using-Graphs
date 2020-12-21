#
#
#
#******************************************************************************************************
#*************************Data for the Chinese Village network*****************************************
#******************************************************************************************************
#
#
#
#Data taken from cai-data should contain data from the paper'Social networks and the decision to insure'(2015). 
#The data is available for download from the publication website(https://www.aeaweb.org/articles?id=10.1257/app.20130442).
cai_all=read.dta('0422allinforawnet.dta')
cai_survey=read.dta('0422survey.dta')
cai_all[1:10,]
cai_edgelist = cai_all %>% 
  dplyr::filter(!is.na(network_id), network_id != 99, id %in% cai_survey$id, network_id %in% cai_survey$id) %>%
  dplyr::select(id, network_id)

ids=with(cai_edgelist,unique(c(id,network_id)))

g=cai_edgelist %>%
  as.matrix %>%
  apply(2,as.character) %>%
  graph_from_edgelist(directed=FALSE)

attributes = cai_survey %>% dplyr::filter(id %in% ids) %>% dplyr::select(id, takeup_survey, delay, intensive)

id_matcher = match(V(g)$name, attributes$id)
V(g)$y = attributes$takeup_survey[id_matcher]
V(g)$delay = attributes$delay[id_matcher]
V(g)$intensive = attributes$intensive[id_matcher]
V(g)$w = 1*(V(g)$intensive == 1)

y=V(g)$y
w=V(g)$w
#
#
#
#******************************************************************************************************
#*************************Sampling and estimation in Chinese Village***********************************
#******************************************************************************************************
#
#
#
#estimates in the complete network
x_obs=.build_obs_covariates(covariate_fns, g, w)
x_trt=.build_treated_covariates(covariate_fns, g)
x_ctrl=.build_control_covariates(covariate_fns, g)
data = list(y=y, w=w, x_obs=x_obs, x_trt=x_trt, x_ctrl=x_ctrl)
n = length(data$y)
n_folds = 5
fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
la_g<-linear_adjustment(data)
dm_g<-difference_in_means(data)
ha_g<-hajek(data,g,'frac1', threshold=0.75)
#estimates in the sampled subgraphs
mc<-40
results_la<-matrix(NA,nrow=50,ncol=mc)
results_dm<-matrix(NA,nrow=50,ncol=mc)
results_ha<-matrix(NA,nrow=50,ncol=mc)
for (j in 1:mc){
  df_sg<-data.frame()  
  sample_graphs<-get_sg(g,FUN = random_edge)
  n_samples<-length(sample_graphs)
  for(i in 1:n_samples){
    sg<-sample_graphs[[i]]
    
    id_matcher = match(V(sg)$name, attributes$id)
    V(sg)$y = attributes$takeup_survey[id_matcher]
    V(sg)$intensive = attributes$intensive[id_matcher]
    V(sg)$w = 1*(V(sg)$intensive == 1)
    
    y_sg=V(sg)$y
    w_sg=V(sg)$w
    x_obs_sg=.build_obs_covariates(covariate_fns, sg, w_sg)
    
    #new_treatment<-generate_covariate_data(sg,covariate_fns)
    #w_sg_new<-new_treatment$w
    #x_obs_sg_new<-new_treatment$x_obs
    x_trt_sg=.build_treated_covariates(covariate_fns, sg)
    x_ctrl_sg=.build_control_covariates(covariate_fns, sg)
    ####
    rows<-which(x_obs_sg$frac1=='NaN')
    if(length(rows)>0){
      x_obs_sg<-x_obs_sg[-(rows),]
      x_trt_sg<-x_trt_sg[-(rows),]
      x_ctrl_sg<-x_ctrl_sg[-(rows),]
      y_sg<-y_sg[-(rows)]
      w_sg<-w_sg[-(rows)] 
    }
    ###
    data_sg = list(y=y_sg, w=w_sg, x_obs=x_obs_sg, x_trt=x_trt_sg, x_ctrl=x_ctrl_sg)
    n = length(data_sg$y)
    df_sg[i,1]<-linear_adjustment(data_sg)
    df_sg[i,2]<-difference_in_means(data_sg)
    df_sg[i,3]<-hajek(data_sg,sg,'frac1', threshold=0.75)
  }
  colnames(df_sg)<-c("linear","diff_mean","hajek")
  results_la[,j]<-df_sg$linear
  results_dm[,j]<-df_sg$diff_mean
  results_ha[,j]<-df_sg$hajek
}
#plotting the results
p<-seq(0.25,0.99,length.out = 50)
par(mfrow=c(1,1),bg="grey100",box(which="plot"))
plot(presults_la[,1],type='l',col='white',
     ylim=c(0.025,0.175),lwd=2,lty=2,xlab = "Proportion of edges sampled from the complete graph",
     ylab = "GATE estimates")
grid(20,7,col="grey75")
smoother<-loess.smooth(1:n_samples,apply(results_la,1,mean),span=2/3)
lines(p,smoother$y,col="red")
abline(h=la_g,col="red",lwd=2,lty=2)
abline(h=dm_g,col="green",lwd=2,lty=2)
abline(h=ha_g,col="darkgreen",lwd=2,lty=2)
smoother<-loess.smooth(1:n_samples,apply(results_dm,1,mean),span=2/3)
lines(p,smoother$y,col="green")
smoother<-loess.smooth(1:n_samples,apply(results_ha,1,mean),span=2/3)
lines(p,smoother$y,col="darkgreen")
legend(2.5,0.225,
       legend=c("linear","diff_mean","hajek"),
       col = c("red","green","darkgreen"),pch=16)
#
#
#
#******************************************************************************************************
#*************************Comparison of distribution of covariates*************************************
#******************************************************************************************************
#
#
#
cv<-1   #column number of covariates
get_sample_graphs<-function(graph,hops=30,FUN){
  p<-seq(0.25,0.95,length.out = 9)
  sample_graphs<-list()
  for (i in 1:9){
    sg<-FUN(g,p[i])
    sample_graphs[[i]]<-sg
  }
  return(sample_graphs)
}
sample_graphs<-get_sample_graphs(g,FUN=random_edge)   #Example shown only for random edge method
n_samples<-length(sample_graphs)
par(mfrow=c(3,3))
for(i in 1:n_samples){
  sg<-sample_graphs[[i]]
  
  id_matcher = match(V(sg)$name, attributes$id)
  V(sg)$y = attributes$takeup_survey[id_matcher]
  V(sg)$intensive = attributes$intensive[id_matcher]
  V(sg)$w = 1*(V(sg)$intensive == 1)
  
  y_sg=V(sg)$y
  w_sg=V(sg)$w
  
  new_treatment<-generate_covariate_data(sg,covariate_fns)
  w_sg_new<-new_treatment$w
  
  x_obs_sg=.build_obs_covariates(covariate_fns, sg, w_sg)
  x_obs_maingraph=.build_obs_covariates(covariate_fns, g, w)
  x_obs_sg_new<-new_treatment$x_obs
  
  x1<-x_obs_maingraph[,cv]
  x2<-x_obs_sg[,cv]
  #x2<-x_obs_sg_new[,cv]
  p1 <- hist(x1)
  p2 <- hist(x2)                     
  x_max<-max(c(x1,x2),na.rm=TRUE)
  y_max<-max(c(p1$density,p2$density),na.rm=TRUE)  
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,x_max), ylim=c(0,y_max),freq = F,main = paste("Covariate",cv,"-(Main)vs(sampled_sg)")) 
  plot( p2, col=rgb(1,0,0,1/4), add=T,freq = F)
}