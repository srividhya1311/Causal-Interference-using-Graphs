#
#
#
#******************************************************************************************************
#************************* Functions for sampling from graphs******************************************
#******************************************************************************************************
#
#
#
random_node<-function(graph,sample_size){
  adj_mat<-as_adjacency_matrix(graph)
  random_nodes=sample(vcount(graph), vcount(graph)*sample_size)  
  sample_adj_mat<-adj_mat[random_nodes,random_nodes]
  sampled_graph<-graph_from_adjacency_matrix(sample_adj_mat,mode='undirected')
  return(sampled_graph)
}
random_degree_node<-function(graph,sample_size){
  adj_mat<-as_adjacency_matrix(graph)
  probs<-degree(graph)/max(degree(graph))
  random_nodes=sample(vcount(graph), vcount(graph)*sample_size,prob = probs)  
  sample_adj_mat<-adj_mat[random_nodes,random_nodes]
  sampled_graph<-graph_from_adjacency_matrix(sample_adj_mat,mode='undirected')
  return(sampled_graph)
}
random_rank_node<-function(graph,sample_size){
  adj_mat<-as_adjacency_matrix(graph)
  probs<-page_rank(graph)$vector/sum(page_rank(graph)$vector)
  random_nodes=sample(vcount(graph), vcount(graph)*sample_size,prob = probs)  
  sample_adj_mat<-adj_mat[random_nodes,random_nodes]
  sampled_graph<-graph_from_adjacency_matrix(sample_adj_mat,mode='undirected')
  return(sampled_graph)
}
random_edge<-function(graph,sample_size){
  sam_edges<-sample(E(graph),length(E(graph))*sample_size)
  c1<-ends(graph,sam_edges)[,1]
  c2<-ends(graph,sam_edges)[,2]
  edge_list<-cbind(c1,c2)
  sampled_graph<-graph_from_edgelist(edge_list,directed=FALSE)
  return(sampled_graph)
}
random_node_edge<-function(graph,sample_size){
  n<-vcount(graph)*sample_size
  random_nodes=sample(vcount(graph), n)
  sam_edges<-c()
  for(i in 1:n){
    edges<-incident_edges(graph,random_nodes[i])[[1]]
    sam_edges[i]<-sample(edges,1)  
  }
  c1<-ends(graph,sam_edges)[,1]
  c2<-ends(graph,sam_edges)[,2]
  edge_list<-cbind(c1,c2)
  sampled_graph<-graph_from_edgelist(edge_list,directed=FALSE)
  return(sampled_graph)
}
random_node_neighbor<-function(graph,sample_size){
  adj_mat<-as_adjacency_matrix(graph)
  n<-vcount(graph)*sample_size/mean(degree(graph))
  random_nodes<-sample(vcount(graph), n)  
  all_nodes<-adjacent_vertices(graph,random_nodes)
  vertices<-random_nodes
  for(i in 1:n){
    n_neigh<-length(all_nodes[[i]])
    list_of_nodes<-all_nodes[[i]]
    for (j in 1:n_neigh)
      vertices<-c(vertices,list_of_nodes[j])
  }
  sample_adj_mat<-adj_mat[vertices,vertices]
  sampled_graph<-graph_from_adjacency_matrix(sample_adj_mat,mode='undirected')
  return(sampled_graph)
}
random_node_jump<-function(graph,sample_size){
  n<-vcount(graph)*sample_size
  vertices<-c()
  initial_node<-sample(vcount(graph), 1)
  vertices[1]<-random_walk(graph,initial_node,1)
  for(i in 2:n){
    start_node<-sample(c(vertices[i-1],sample(vcount(graph),1)),1,prob = c(0.85,0.15))
    vertices[i]<-random_walk(graph,start_node,2)[2]
  }
  adj_mat<-as_adjacency_matrix(graph)
  sample_adj_mat<-adj_mat[vertices,vertices]
  sampled_graph<-graph_from_adjacency_matrix(sample_adj_mat,mode='undirected')
  return(sampled_graph)
}
random_node_walk<-function(graph,sample_size){
  n<-vcount(graph)*sample_size
  vertices<-c()
  initial_node<-sample(vcount(graph), 1)
  vertices[1]<-random_walk(graph,initial_node,1)
  for(i in 2:n){
    start_node<-sample(c(vertices[i-1],initial_node),1,prob = c(0.85,0.15))
    vertices[i]<-random_walk(graph,start_node,2)[2]
  }
  adj_mat<-as_adjacency_matrix(graph)
  sample_adj_mat<-adj_mat[vertices,vertices]
  sampled_graph<-graph_from_adjacency_matrix(sample_adj_mat,mode='undirected')
  return(sampled_graph)
}