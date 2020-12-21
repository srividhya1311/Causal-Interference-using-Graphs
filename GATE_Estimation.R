#
#
#
#******************************************************************************************************
#*************************Functions to get covariates for regression adjustment************************
#******************************************************************************************************
#based on reference code at https://github.com/ajchin/regression-adjustments
#to reproduce the results in Chin (2019)
#
#
#
fraction_trt_nbrs = function(g, w) {as.vector(as_adj(g) %*% w) / degree(g)}
number_trt_nbrs = function(g, w) {as.vector(as_adj(g) %*% w)}
fraction_trt_nbrs2 = function(g, w) {
  adj = as_adj(g)
  adj2 = adj %*% adj
  as.vector(adj2 %*% w) / apply(adj2, 2, sum) 
}
fraction_trt_nbrs3 = function(g, w) {
  adj = as_adj(g)
  adj3 = adj %*% adj %*% adj
  as.vector(adj3 %*% w) / apply(adj3, 2, sum)
}
number_trt_nbrs2 = function(g, w) {
  adj = as_adj(g)
  adj2 = adj %*% adj
  as.vector(adj2 %*% w)
}
#Functions for data generation
.build_obs_covariates = function(covariate_fns, g, w) {
  x_obs = lapply(covariate_fns, function(f) {f(g, w)})
  do.call(data.frame, x_obs)
}
.build_counterfactual_covariates = function(covariate_fns, g, w){
  x = lapply(covariate_fns, function(f) {f(g, w)})
  do.call(data.frame, x)
}
.build_treated_covariates = function(covariate_fns, g) {
  .build_counterfactual_covariates(covariate_fns, g, rep(1, vcount(g)))#, '_T')
}
.build_control_covariates = function(covariate_fns, g) {
  .build_counterfactual_covariates(covariate_fns, g, rep(0, vcount(g)))#, '_C')
}
generate_covariate_data = function(g, covariate_fns) {
  w = rbinom(vcount(g), size=1, prob=0.5)
  x_observed = .build_obs_covariates(covariate_fns, g, w) 
  x_global_trt = .build_treated_covariates(covariate_fns, g)
  x_global_ctrl = .build_control_covariates(covariate_fns, g)
  list(w=w, x_obs=x_observed, x_trt=x_global_trt, x_ctrl=x_global_ctrl)
}
dynamic_lim_response<-function(w,g,param){
  adj_mat=as_adj(g)
  n_peers=degree(g)
  n=length(w)
  y=rep(0,n)
  
  max_t=param$max_t
  noise_sd=param$noise_sd
  int=param$int
  direct=param$direct
  spill=param$spill
  
  for(t in 1:max_t){
    avg_nbr_y=as.vector(adj%*%y/n_peers)
    y=int+direct*w+spill*avg_nbr_y+rnorm(n,sd=noise_sd)
    y=1*(y>0)
  }
  return(y)
}
get_response<-function(w,g,param){
  y<-matrix(NA,nrow=vcount(g),ncol=5)
  for (i in 1:5){
    y[,i]=dynamic_lim_response(w,g,param[[i]])
  }
  return(y)
}
covariate_fns=list(
  frac1 = fraction_trt_nbrs,
  frac2 = fraction_trt_nbrs2,
  num1 = number_trt_nbrs,
  num2 = number_trt_nbrs2
)
#
#
##******************************************************************************************************
#*************************Functions to get GATE Estimates**********************************************
#******************************************************************************************************
#
#
#
difference_in_means = function(data) {
  with(data, mean(y[w==1]) - mean(y[w==0]))
}
.indiv_hajek_weight = function(w, x, d, thresh, p_design=0.5) {
  if (w == 0) {
    if (x > 1 - thresh) return(0)
    p = (1 - p_design) * pbinom(floor(d * (1 - thresh)), size=d, prob=p_design)
    return(1 / p)
  }
  if (x < thresh) return(0)
  p = p_design * (1 - pbinom(floor(d * thresh), size=d, prob=p_design))
  return(1 / p)
}
hajek_weights = function(data, g, threshold_var_name, threshold, p_design=0.5) {
  w = data$w
  threshold_var = data$x_obs[, threshold_var_name]
  wts = sapply(1:nrow(data$x_obs), function(i) {
    .indiv_hajek_weight(w[i], threshold_var[i], degree(g)[i], threshold, p_design)
  })
  sum_trt_wt = sum(wts[w == 1])
  sum_ctrl_wt = sum(wts[w == 0])
  return(wts / sum_trt_wt * (w == 1) - wts / sum_ctrl_wt * (w == 0))
}
hajek = function(data, g, threshold_var_name, threshold, p_design=0.5) {
  wts = hajek_weights(data, g, threshold_var_name, threshold, p_design)
  return(sum(data$y * wts))
}
linear_adjustment = function(data, vars=NULL) {
  if (is.null(vars)) vars = names(data$x_obs)
  w = data$w
  y0 = data$y[w==0]
  y1 = data$y[w==1]
  x0 = data$x_obs %>% select(one_of(vars)) %>% filter(w==0)
  x1 = data$x_obs %>% select(one_of(vars)) %>% filter(w==1)
  beta0 = lm(y0 ~ ., data=x0) %>% coef
  beta1 = lm(y1 ~ ., data=x1) %>% coef
  omega0 = c(1, colMeans(data$x_ctrl %>% select(one_of(vars))))
  omega1 = c(1, colMeans(data$x_trt %>% select(one_of(vars))))
  sum(omega1 * beta1) - sum(omega0 * beta0)
}
get_estimates<-function(g,covariate_fns,params){
  cv_data<-generate_covariate_data(g,covariate_fns)
  x_obs=cv_data$x_obs
  response<-get_response(cv_data$w,g,params)
  rows<-which(x_obs$frac1=='NaN')
  rows_y<-which(is.na(response[,1]))
  if(length(rows)>0){
    x_obs<-cv_data$x_obs[-(rows),]
    x_trt<-cv_data$x_trt[-(rows),]
    x_ctrl<-cv_data$x_ctrl[-(rows),]
    y<-response[-(rows),]
    w<-cv_data$w[-(rows)] 
    rows_y<-which(is.na(y[,1]))
    x_obs<-x_obs[-(rows_y),]
    x_trt<-x_trt[-(rows_y),]
    x_ctrl<-x_ctrl[-(rows_y),]
    y<-y[-(rows_y),]
    w<-w[-(rows_y)] 
    data = list(y=y,w=w, x_obs=x_obs, x_trt=x_trt, x_ctrl=x_ctrl)
  } else if (length(rows_y)>0){
    x_obs<-cv_data$x_obs[-(rows_y),]
    x_trt<-cv_data$x_trt[-(rows_y),]
    x_ctrl<-cv_data$x_ctrl[-(rows_y),]
    y<-response[-(rows_y),]
    w<-cv_data$w[-(rows_y)] 
    data = list(y=y,w=w, x_obs=x_obs, x_trt=x_trt, x_ctrl=x_ctrl)    
  } else {
    data = list(y=response,w=cv_data$w, x_obs=cv_data$x_obs, x_trt=cv_data$x_trt, x_ctrl=cv_data$x_ctrl)    
  }
  la_est<-c()
  dm_est<-c()
  ha_est<-c()
  for(i in 1:5){
    y_tmp=data$y[,i]
    data_param=list(y=y_tmp,w=data$w,x_obs=data$x_obs,x_trt=data$x_trt,x_ctrl=data$x_ctrl)
    la_est[i]<-linear_adjustment(data_param)
    dm_est[i]<-difference_in_means(data_param)
    ha_est[i]<-hajek(data_param,g,'frac1', threshold=0.75)
  }
  return(list(la_est,dm_est,ha_est))
}
#
#
#