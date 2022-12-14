#PDRWH
library(igraph)

get_mut_data=function(data_file){
  
  data=read.delim(data_file,header = T,as.is = T,check.names = F)
  especial_gene_index=which(is.na(data[,1])|data[,1]=="")
  if (length(especial_gene_index)>0){
    data=data[-especial_gene_index,]
  } else {data=data }
  
  rownames(data)=data[,1]
  data=data[,-1]
  mut_data=as.matrix(data[which(apply(data,1,sum)!=0),which(apply(data,2,sum)!=0)])
  
  return(mut_data)
}

get_vertex_W=function(sample_i,data_mut_com,com_mut_samples,Graph,sample_i_mut_genes,H){
  
  mut_data_i=data_mut_com[,sample_i]
  vertex_W=c()
  for(k in com_mut_samples){
    vertex_w=rep(0,length(sample_i_mut_genes))
    names(vertex_w)=sample_i_mut_genes
    
    mut_data_k=data_mut_com[,k] 
    mut_gene=names(mut_data_k[which(mut_data_k==1)])
    
    v_i=intersect(mut_gene,V(Graph)$name)
    v_diff=setdiff(mut_gene,V(Graph)$name)
    
    graph_i=induced_subgraph(Graph, v_i) 
    graph_i_degree=degree(graph_i)[v_i]
    graph_i_degree[which(graph_i_degree==0)]=0.01
    
    co_mut_genes=intersect(v_i,sample_i_mut_genes)
    
    vertex_w[co_mut_genes]=graph_i_degree[co_mut_genes]
    vertex_w[setdiff(sample_i_mut_genes,co_mut_genes)]=0.01
    vertex_W=cbind(vertex_W,vertex_w)
  }
  
  colnames(vertex_W)=com_mut_samples
  return(vertex_W)
}

get_P=function(H,vertex_W){
  
  Degree_v=apply(H,1,sum)
  D_v_inverse=diag(1/Degree_v)
  
  probability_hyperedge=D_v_inverse%*%H
  rownames(probability_hyperedge)=rownames(H)
  
  
  Degree_ve=apply(vertex_W,2,sum)
  D_ve_inverse=diag(1/Degree_ve)
  
  probability_vertex=D_ve_inverse%*%t(vertex_W)
  rownames(probability_vertex)=colnames(H)
  
  P = probability_hyperedge%*%probability_vertex
  return(as.matrix(P))
}

get_hyper_randomwalk=function(P,H,theta=0.85) {
  
  v0=rep(1/nrow(H),nrow(H))    #initial vector
  teleport=rep(1/nrow(H),nrow(H)) #
  
  Distance=c()
  vi=v0
  for(k in 1:20){
    
    vj=vi
    vi=theta*t(P)%*%vi+(1-theta)*teleport
    dis=sum(abs(vj-vi))
    Distance=append(Distance,dis)
    if (dis < 0.000001) {
      break
    }
  }
  plot(1:k,Distance)
  ##
  
  Importance=data.frame(vi[order(vi,decreasing = T),])
  colnames(Importance)=c("Importance")
  return(Importance)
}

hygraph_randomwalk=function(sample_i,data_mut_com,co_mut,Graph){
  
  sample_i_mut_genes=rownames(data_mut_com)[which(data_mut_com[,sample_i]==1)] 
  
  com_mut_sample_i=co_mut[sample_i,which(co_mut[sample_i,]!=0)] 
  com_mut_samples=names(com_mut_sample_i)
  
  H=data_mut_com[sample_i_mut_genes,com_mut_samples,drop=F] 
  
  Vertex_W=get_vertex_W(sample_i,data_mut_com,com_mut_samples,Graph,sample_i_mut_genes,H)
  
  P=get_P(H,Vertex_W)
  importance_score = get_hyper_randomwalk(P,H,theta=0.85)
  
  return(importance_score)
}

PDRWHscore=function(cancer_type,mut_data_file,network_file,outfile_dir){
  
  mut_data=get_mut_data(mut_data_file) 
  mut_data_morethan_2=colnames(mut_data)[which(apply(mut_data,2,sum)>=3)] #
  Network=read.table(network_file,header = T)
  Graph=graph.data.frame(Network) 
  
  com_samples=mut_data_morethan_2
  data_mut_com=mut_data[,com_samples] 
  #
  co_mut <- t(data_mut_com) %*% data_mut_com
  rownames(co_mut) <- com_samples
  colnames(co_mut) <- com_samples
  
  Importance_Score=list()
  for (i in 1:length(com_samples)){
    sample_i=com_samples[i]
    print(paste0(i," ",sample_i))
    importance_score=hygraph_randomwalk(sample_i,data_mut_com,co_mut,Graph)
    Importance_Score[[i]]=importance_score
  }
  names(Importance_Score)=com_samples
  
  save(Importance_Score,file=paste(outfile_dir, cancer_type,".Rdata",sep=""))
  #load(file = rdata_dir)
  #print("Congratrulations!")
}

resultTransform=function(cancer_type,input_dir,outfile_dir){
  load(paste0(input_dir,cancer_type,".Rdata"))
  
  all_genes <- c()
  for(i in (1:length(Importance_Score[]))) {
    genes <- rownames(Importance_Score[[i]])
    all_genes <- union(all_genes,genes)
  }
  
  res_matrix <- matrix(0,length(all_genes),length(Importance_Score[]))
  rownames(res_matrix) <- all_genes
  colnames(res_matrix) <- names(Importance_Score[])
  for(i in (1:length(Importance_Score[]))) {
    for (j in (1:length(rownames(Importance_Score[[i]])))) {
      res_matrix[rownames(Importance_Score[[i]])[j],i] <- Importance_Score[[i]][j,]
    }
  }  
  
  write.table(res_matrix,paste0(outfile_dir,"/PDRWH.txt"),quote = T,sep="\t")
}
