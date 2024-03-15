#
library('RColorBrewer')
library('ggplot2')
library("readxl")

mycolors<-c("darkred","#0072BD","#D95319","#EDB120","#4DBEEE")
#setwd("../")

get_data <- function(data_file) {
  
  data <- read.delim(data_file,header = T,as.is = T,check.names = F)
  especial_gene_index <- which(is.na(data[,1])|data[,1]=="")
  if (length(especial_gene_index) > 0){
    data[especial_gene_index,1] <- "NA"
  } else {data <- data }
  
  #data=data[-which(data[,1]==""),]
  rownames(data) <- data[,1]
  data <- data[,-1]
  data_mut <- as.matrix(data[which(apply(data,1,sum) != 0),which(apply(data,2,sum) != 0)])
  return(data_mut)
}

get_R_P_F <- function(score_data,genelist,data_mut,top_n,method,mut_fre,mut_degree,PR_type=1) {
  
  tmp <- colnames(data_mut)[2]
  if (nchar(tmp) > 12 & method == "PDRWH") {
    colnames(data_mut) <- substr(colnames(data_mut),1,12)
  }
  if (dim(score_data)[2] == 1) {
    samples <- colnames(data_mut)
  } else {
    samples <- colnames(score_data)
  }
  
  P_R_F1 <- list()
  for(k in 1:length(samples)){
    print(paste0(k,"  ",samples[k]))
    if(nchar(samples[k]) < 4){
      mut_data <- score_data[,samples[k]]
    } else if (PR_type == 1) {
      mut_data <- data_mut[,samples[k]]
      mut_gene <- names(mut_data[mut_data != 0])
    } else if (PR_type == 2) {
      mut_data <- score_data[,samples[k],drop=F]
      mut_gene <- rownames(mut_data)[which(mut_data != 0)]
    }
    
    score_label <- data.frame(matrix(0,nrow = length(mut_gene),ncol = 2))
    rownames(score_label) <- mut_gene
    com_genes <- intersect(genelist,mut_gene)
    if (dim(score_data)[2] == 1) {
      score_label[,1] <- score_data[mut_gene,1]
    } else {
      score_label[,1] <- score_data[mut_gene,samples[k]]
    }
    score_label[com_genes,2] <- 1
    colnames(score_label) <- c("score","label")
    score_label <- score_label[order(score_label[,1],decreasing = T),]
    
    n_com_genes <- length(com_genes)
    if (n_com_genes != 0) {
      Recall <- c()
      Precision <- c()
      F1 <- c()
      for(i in 1:nrow(score_label)){
        if (i > top_n) {
          break
        }
        com_n <- length(intersect(genelist,rownames(score_label)[1:i]))
        precision <- com_n / i
        recall <- com_n / n_com_genes
        if (precision * recall != 0){
          f1 <- 2 * precision * recall / (precision + recall)
        } else {f1 <- 0}
        
        Precision <- c(Precision,precision)
        Recall <- c(Recall,recall)
        F1 <- c(F1,f1)
      }
      P_R_F1_sample <- rbind(Precision,Recall,F1)
      
      if (nrow(score_label) < top_n) {
        P_R_F1[[k]] <- cbind(P_R_F1_sample,matrix(-1,3,top_n-ncol(P_R_F1_sample)))
      } else {
        P_R_F1[[k]] <- P_R_F1_sample
      }
      
    } else {
      P_R_F1[[k]] <- ""
    }
  }
  names(P_R_F1) <- samples
  
  null_names <- c()
  for(i in 1:length(P_R_F1)){
    if(is.character(P_R_F1[[i]])){
      null_names <- c(null_names,names(P_R_F1)[i])
    }
  }
  if (length(null_names)) {
    R_P_F_s <- P_R_F1[setdiff(names(P_R_F1),null_names)]
  } else {
    R_P_F_s <- P_R_F1
  }
  
  P_TOP <- c()
  R_TOP <- c()
  F_TOP <- c()
  for(i in 1:length(R_P_F_s)){
    P_TOP <- rbind(P_TOP,R_P_F_s[[i]][1,,drop=F])
    R_TOP <- rbind(R_TOP,R_P_F_s[[i]][2,,drop=F])
    F_TOP <- rbind(F_TOP,R_P_F_s[[i]][3,,drop=F])
  }
  # average_p <- apply(P_TOP,2,mean)
  # average_r <- apply(R_TOP,2,mean)
  # average_f <- apply(F_TOP,2,mean)
  P_TOP_idx <- (P_TOP > -1) + 0
  R_TOP_idx <- (R_TOP > -1) + 0
  F_TOP_idx <- (F_TOP > -1) + 0
  average_p <- apply(P_TOP * P_TOP_idx,2,sum) / apply(P_TOP_idx,2,sum)
  average_r <- apply(R_TOP * R_TOP_idx,2,sum) / apply(R_TOP_idx,2,sum)
  average_f <- apply(F_TOP * F_TOP_idx,2,sum) / apply(F_TOP_idx,2,sum)
  
  return(rbind(average_p,average_r,average_f))
}

get_top_n <- function(methods_nums,score_data,methods_names,genelist,data_mut,genelist_name,top_N) {
  
  plot_data <- list()
  for(i in (1:methods_nums)) {
    print(paste0(i,"  ",methods_names[i]))
    method_i_plot_data <- get_R_P_F(score_data[[i]],genelist,data_mut,top_N,methods_names[i])
    plot_data[[i]] <- method_i_plot_data
  }
  
  pdf(paste("../out/",genelist_name,"_topN.pdf",sep=""),width = 9.6,height = 3.2)
  par(mfrow = c(1,3))
  m <- seq(from = 1, to = top_N, by = 1)
  
  plot(plot_data[[1]][1,],type = "o",ylim = c(0,1),xlim = c(1,top_N),pch = 16,col = mycolors[1],
       cex = 0.6,frame.plot = F,xlab = "Top N genes",ylab = "Average precision",
       xaxt="n")
  for ( i in (1:methods_nums)) {
    points(plot_data[[i]][1,],type = "o",pch = 15 + i,col = mycolors[i],lty = i,cex = 0.6)
  }
  legend(6,1,methods_names,col = mycolors[1:methods_nums],text.col = mycolors[1:methods_nums],
         lty = c(1:methods_nums),pch = c(16:(15+methods_nums)),bty = "n",cex = 0.8)
  axis(1,m)
  
  plot(plot_data[[1]][2,],type = "o",ylim = c(0,1),xlim = c(1,top_N),pch = 16,col = mycolors[1],
       cex = 0.6,frame.plot = F,xlab = "Top N genes",ylab = "Average recall",main = genelist_name,
       xaxt="n")
  for ( i in (1:methods_nums)) {
    points(plot_data[[i]][2,],type = "o",pch = 15 + i,col = mycolors[i],lty = i,cex = 0.6)
  }
  axis(1,m)
  
  plot(plot_data[[1]][3,],type = "o",ylim = c(0,0.5),xlim = c(1,top_N),pch = 16,col = mycolors[1],
       cex = 0.6,frame.plot = F,xlab = "Top N genes",ylab = "Average F1",
       xaxt="n")
  for ( i in (1:methods_nums)) {
    points(plot_data[[i]][3,],type = "o",pch = 15 + i,col = mycolors[i],lty = i,cex = 0.6)
  }
  axis(1,m)
  
  dev.off()
}

process_prodigy <- function(results,data_mut) {
  
  genes <- union(unique(unlist(results)),rownames(data_mut))
  prodigy <- matrix(0,length(genes),length(results))
  rownames(prodigy) <- genes
  colnames(prodigy) <- names(results)
  for( i in 1:length(results)){
    if(length(results[[i]]) != 0){
      for (j in 1:length(results[[i]])){
        prodigy[results[[i]][j],i] = length(results[[i]]) + 1 - j
      }
    }else{}
  }
  return(prodigy)
}


###plot
genelist <- read.csv("../data/general_driver_list.csv",header = T)[,1]
cancer_list <- c("BRCA","KIRC","LIHC","GBM","STAD")

for (k in 1:length(cancer_list)) {
  cancer <- cancer_list[k]
  print("----------------------------------------------------")
  print(paste0("Precision, Recall, F1-score in ",cancer,"..."))
  print("----------------------------------------------------")
  prdwh <- read.table(paste0("../out/",cancer,"/PDRWH.txt"),header = T,check.names = F)
  data_mut <- get_data(paste0("../data/",cancer,"/",cancer,"_mc3_gene_level.txt"))
  colnames(data_mut) <- substr(colnames(data_mut),1,12)
  score_data <- list(prdwh)

  get_top_n(methods_nums = 1,
            score_data = score_data,
            methods_names = "PDRWH",
            genelist,data_mut,paste0(cancer,""),
            top_N = 8)
  if (k == length(cancer_list)) {
    print("----------------------------------------------------")
    print(paste0("Precision, Recall, F1-score in ",cancer,": Achieved!"))
    print("----------------------------------------------------")
  }
}


