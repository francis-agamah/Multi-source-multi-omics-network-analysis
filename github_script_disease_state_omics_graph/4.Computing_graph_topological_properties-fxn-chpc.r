#this function accept network file and disease state
#e.g topo_fxn("network-file","mild")

#X is the network file
topo_fxn<-function(x,state){
  #hG<-read.table(file.path(prepped_datas,"All_prot_trans_and_metabo_case_significant.txt"),header = T,sep = "\t")
  hG<-read.table(x,header = T,sep = "\t")
  hG<-graph.data.frame(d=hG,directed = F)
  node_degree<-degree(graph = hG,loops = F)
  node_degree1<- as.data.frame(as.table(node_degree)) #converting to a table and data frame
  colnames(node_degree1)=c("feature","node_degree")
  degree.distribution(graph = hG,cumulative = T);
  
  
  #Betweenness centrality of the edges
  estimate_betweenness(hG,directed = F,cutoff = -1);
  node_betwenness<-betweenness(hG,directed = F,normalized = F,nobigint = T,weights = NULL) #betweenness calculates vertex betweenness, 
  node_betwenness1<- as.data.frame(as.table(node_betwenness)) #converting to a table and data frame
  colnames(node_betwenness1)=c("feature","node_betweeness")
  edge_betweenness<-edge_betweenness(hG,e = E(hG), directed = F); #edge_betweenness calculates edge betweenness
  
  #closeness
  node_closeness<-closeness(hG,normalized = F) #Closeness centrality measures how many steps is required to access every other vertex from a given vertex.
  node_closeness1<- as.data.frame(as.table(node_closeness)) #converting to a table and data frame
  colnames(node_closeness1)=c("feature","node_closeness")
  #eigenvector
  node_eigen<-eigen_centrality(hG,directed = F,scale = T,options = arpack_defaults)$vector
  node_eigen1<- as.data.frame(as.table(node_eigen)) #converting to a table and data frame
  colnames(node_eigen1)=c("feature","node_eigen")
  #hist(centrality_metric_val$node_eigen, main="eigenvector distribution in healthy",xlab = "eigenvector", ylab = "Number of nodes",col="Pink")
  
  #centrality_metric<-c("node_degree","node_betwenness","node_closeness","node_eigen")
  centrality_metric_val<-data.frame(cbind(node_degree,node_betwenness,node_closeness,node_eigen))
  centrality_metric_val_table<-data.frame(node_degree1$feature, node_degree,node_betwenness1$node_betweeness,node_closeness1$node_closeness,node_eigen1$node_eigen, row.names = NULL)
  colnames(centrality_metric_val_table)<-c("feature","node_degree","node_betwenness", "node_closeness", "node_eigen1")
  #computing log of node and betweeness centrality metric
  
  #'''
  #Computiation of test statistic for node with no log transformation of node and betweeness centrality
  #'''
  #Bionetstat concept: computing test statistics for node comparison
  node_centrality_mean<-rowMeans(centrality_metric_val)
  #compute the distance between the centrality of nodes in graph and the average node centrality
  dist_between_centrality_of_nodes_for_node_degree<-abs(node_degree-node_centrality_mean)
  dist_between_centrality_of_nodes_for_node_betwenness<-abs(node_betwenness-node_centrality_mean)
  dist_between_centrality_of_nodes_for_node_closeness<-abs(node_closeness-node_centrality_mean)
  dist_between_centrality_of_nodes_for_node_eigen<-abs(node_eigen-node_centrality_mean)
  dist_matrix<-data.frame(cbind(dist_between_centrality_of_nodes_for_node_degree,dist_between_centrality_of_nodes_for_node_betwenness,
                                dist_between_centrality_of_nodes_for_node_closeness,dist_between_centrality_of_nodes_for_node_eigen  ))
  
  node_test_statistic<-rowMeans(dist_matrix) #this statistic measures the difference among centralities for each node
  node_test_statistic<-as.table(node_test_statistic)
  node_test_statistic<-data.frame(node_test_statistic)
  colnames(node_test_statistic)<-c("feature","statistic")
  node_test_statistic<-node_test_statistic[rev(order(test_fxn$node_test_statistic$statistic)), ]
  #node_test_statistic<-node_test_statistic[order(node_test_statistic$Freq),]
  #lines(density(centrality_metric_val$node_degree))
  #'''
  #Computiation of test statistic for node with log transformation of node and betweeness centrality
  #'''
  centrality_metric_val$log_node_degree<-log(centrality_metric_val$node_degree)
  centrality_metric_val$log_node_betwenness <-log(centrality_metric_val$node_betwenness)
  centrality_metric_val_log_transformed<-data.frame(cbind(centrality_metric_val$log_node_degree,centrality_metric_val$log_node_betwenness,node_closeness,node_eigen))
  colnames(centrality_metric_val_log_transformed)[1]<-"log_node_degree"
  colnames(centrality_metric_val_log_transformed)[2]<-"log_node_betweeness"
  centrality_metric_val_log_transformed_table<-data.frame(centrality_metric_val_table$feature, centrality_metric_val$node_degree,centrality_metric_val$node_betwenness, centrality_metric_val_log_transformed, row.names = NULL)
  
  node_centrality_mean_log_transformed<-rowMeans(centrality_metric_val_log_transformed)
  #compute the distance between the centrality of nodes in graph and the average node centrality
  dist_between_centrality_of_nodes_for_log_node_degree<-abs(centrality_metric_val$log_node_degree-node_centrality_mean_log_transformed)
  dist_between_centrality_of_nodes_for_log_node_betwenness<-abs(centrality_metric_val$log_node_betwenness-node_centrality_mean_log_transformed)
  dist_between_centrality_of_nodes_for_node_closeness<-abs(node_closeness-node_centrality_mean_log_transformed)
  dist_between_centrality_of_nodes_for_node_eigen<-abs(node_eigen-node_centrality_mean_log_transformed)
  
  dist_matrix_log_transformed<-data.frame(cbind(dist_between_centrality_of_nodes_for_log_node_degree,dist_between_centrality_of_nodes_for_log_node_betwenness,
                                                dist_between_centrality_of_nodes_for_node_closeness,dist_between_centrality_of_nodes_for_node_eigen  ))
  
  node_test_statistic_log_transformed<-rowMeans(dist_matrix_log_transformed) #this statistic measures the difference among centralities for each node
  node_test_statistic_log_transformed<-as.table(node_test_statistic_log_transformed)
  node_test_statistic_log_transformed<-data.frame(node_test_statistic_log_transformed)
  colnames(node_test_statistic_log_transformed)<-c("feature","statistic")
  node_test_statistic_log_transformed<-node_test_statistic_log_transformed[rev(order(node_test_statistic_log_transformed$statistic)), ]
  
  #'''
  #Combining the node centralities for both log transformed and non-log transformed
  #'''
  centrality_metric_val_combined<-as.data.frame(cbind(centrality_metric_val,centrality_metric_val$log_node_degree,centrality_metric_val$log_node_betwenness))
  #my_return_list<-list(centrality_metric_val, centrality_metric_val_combined, node_test_statistic, node_test_statistic_log_transformed)
  my_return_list<-list()
  my_return_list$centrality_metric_val<-centrality_metric_val_table #centrality_metric_val_table
  #my_return_list$centrality_metric_val_combined<-centrality_metric_val_combined
  my_return_list$node_test_statistic<-node_test_statistic
  my_return_list$node_test_statistic_log_transformed<-node_test_statistic_log_transformed
  write.table(node_test_statistic,paste(c("node_test_statistic",state,".txt"),collapse = "_"),sep="\t", col.names=T,row.names=F)
  write.table(node_test_statistic_log_transformed,paste(c("node_test_statistic_log_transformed",state,".txt"),sep="\t",collapse = "_"),col.names=T,row.names=F)
  write.table(centrality_metric_val_log_transformed_table,paste(c("centrality_metric_val",state,".txt"),collapse = "_"),sep="\t", col.names=T,row.names=F)
  
  return(my_return_list)
}

