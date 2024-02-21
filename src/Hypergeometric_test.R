#Hypergeometric test

get_phyper_res <- function(reference_gene_list,individualized_result,snv,k) {
  
  phyper_res <- c()
  individualized_result <- as.matrix(individualized_result)
  for (i in (1:ncol(individualized_result))){
    sample <- colnames(individualized_result)[i]
    mutant_genes <- intersect(rownames(individualized_result),rownames(snv)[which(snv[,sample] != 0)])
    
    temp <- individualized_result[mutant_genes,sample]
    if (k < length(mutant_genes)) {
      select_drivers <- names(temp[order(temp,decreasing = T)])[1:k]      
    } else {
      select_drivers <- names(temp[order(temp,decreasing = T)])[1:length(mutant_genes)]      
    }
    overlapping_genes_between_genelist_and_sample <- intersect(rownames(individualized_result),reference_gene_list)
    overlapping_genes_between_genelist_and_select_drivers <- intersect(select_drivers,reference_gene_list)
    
    N <- nrow(individualized_result)
    M <- length(overlapping_genes_between_genelist_and_sample)
    n <- length(select_drivers)
    m <- length(overlapping_genes_between_genelist_and_select_drivers)
    
    phyper_res <- c(phyper_res,phyper(m - 1,M,(N - M),n,lower.tail = F))
    mut_num <- length(rownames(snv)[which(snv[,sample] != 0)])
    mut_spe <- length(intersect(rownames(snv)[which(snv[,sample] != 0)],reference_gene_list))
  }  
  ratio_of_significant_samples <- 1 - length(which(phyper_res > 0.05)) / (dim(individualized_result)[2])
  
  return(ratio_of_significant_samples)
}

random_get_phyper <- function(reference_gene_list,individualized_result,snv,k) {
  
  phyper_res <- c()
  individualized_result <- as.matrix(individualized_result)
  for (i in (1:ncol(individualized_result))){
    sample <- colnames(individualized_result)[i]
    mutant_genes <- intersect(rownames(individualized_result),rownames(snv)[which(snv[,sample] != 0)])
    
    temp <- individualized_result[mutant_genes,sample]
    if (k < length(mutant_genes)) {
      select_drivers <- names(temp[order(temp,decreasing = T)])[sample(1:length(mutant_genes),k)]      
    } else {
      select_drivers <- names(temp[order(temp,decreasing = T)])[1:length(mutant_genes)]      
    }
    overlapping_genes_between_genelist_and_sample <- intersect(rownames(individualized_result),reference_gene_list)
    overlapping_genes_between_genelist_and_select_drivers <- intersect(select_drivers,reference_gene_list)
    
    N <- nrow(individualized_result)
    M <- length(overlapping_genes_between_genelist_and_sample)
    n <- length(select_drivers)
    m <- length(overlapping_genes_between_genelist_and_select_drivers)
    
    phyper_res <- c(phyper_res,phyper(m - 1,M,(N - M),n,lower.tail = F))
  }  
  ratio_of_significant_samples <- 1 - length(which(phyper_res > 0.05)) / (dim(individualized_result)[2])
  
  return(ratio_of_significant_samples)
}

library('stats')
library('readxl')
library('ggplot2')
library('reshape2')
library('dplyr')

#setwd("../")

cancer_type <- c("BRCA","KIRC","LIHC","GBM","STAD")
TOP <- c(8,10,12,8,16)
cancer_specific_index <- c("G3:G41","Z3:Z24","L3:L33","AE3:AE28","M3:M37")

#for cancer-specific
specific_matrix <- matrix(0,5,2)
colnames(specific_matrix) <- c("PDRWH","random")

#
for (i in 1:length(cancer_type)) {
  str0 <- cancer_type[i]
  print(paste0(cancer_type[i],"..."))
  cancer_specific <- c(t(read_excel("../data/Cancer_specific.xls",col_names = F,range = cancer_specific_index[i])))
   
  pdrwh <- read.table(paste0("../out/",str0,"/PDRWH.txt"),header = T,check.names = F)
  snv <- read.delim2(paste0("../data/",str0,"/",str0,"_mc3_gene_level.txt"),check.names = F)
  {
    snv[is.na(snv)] <- "NA"
    rownames(snv) <- snv[,1]
    snv <- snv[,-1]
    # colnames(snv) <- substr(colnames(snv),1,12)
  }
  #Hypergeometric test
  top_N <- TOP[i]
  #cancer_specific
  specific_matrix[i,"PDRWH"] <- get_phyper_res(reference_gene_list = cancer_specific,
                                   individualized_result = pdrwh,snv = snv,k = top_N)
  specific_matrix[i,"random"] <- random_get_phyper(reference_gene_list = cancer_specific,
                                                  individualized_result = pdrwh,snv = snv,k = top_N)
}

#
comparison <- data.frame(cancer_type,specific_matrix)
names(comparison) <- c("Type","PDRWH","random")
comparison <- melt(comparison,variable.name = "Method",value.name = "Percent")

ggplot(comparison, aes(x = Type, y = Percent, fill = Method)) + 
  geom_bar(stat ="identity",width = 0.6,position = "dodge") + 
  scale_fill_manual(values = c("red3","gray90")) + 
  labs(x = "",y = "", title = "tumor-specific") + 
  ylim(0,1) + 
  guides(fill = guide_legend(reverse = F)) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 20,face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),        
        legend.text = element_text(size = 16, face = "bold"), 
        legend.position = 'top',  
        legend.key.size = unit(0.4,'cm'))

rownames(specific_matrix) <- cancer_type
write.table(specific_matrix,"../out/Hypergeometric test.txt",sep = "\t",quote = F)


