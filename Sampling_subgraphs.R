#
#
#
#******************************************************************************************************
#*************************Functions to get sample graphs **********************************************
#******************************************************************************************************
#
#
#
get_sg<-function(graph,hops=30,FUN){
  p<-seq(0.25,0.99,length.out = 50)
  sample_graphs<-list()
  for (i in 1:50){
    sg<-FUN(graph,p[i])
    sample_graphs[[i]]<-sg
  }
  return(sample_graphs)
}
#
#
#