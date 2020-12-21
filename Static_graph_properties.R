#
#
#
#******************************************************************************************************
#************************* Functions for static graph properties***************************************
#******************************************************************************************************
#
#
#
#get degree distribution
degree_dist<-function(graph){
  G.degrees <- degree(graph)
  G.degree.freq <- as.data.frame(table(G.degrees))
  return(G.degree.freq[,2])
}
#clustering coefficient
clust_coeff<-function(graph){
  l<-transitivity(graph, type="local",isolates = c("zero"))
  cc<-c()
  n<-max(degree(graph))
  for (i in 1:n){
    if (length(which(degree(graph)==i))==0){
      cc[i]=0 
    } else {
      cc[i]<-mean(l[which(degree(graph)==i)]) 
    }
  }
  return(cc)
}
#number of reachable pairs of nodes at distance h = hops 
hop_plot<-function(graph,hops){
  reachable_nodes<-c()
  for (i in 1:hops){
    reachable_nodes[i]<- sum(ego_size(graph,order=i,nodes=V(graph)))
  }
  return(reachable_nodes)
}
#
#
#