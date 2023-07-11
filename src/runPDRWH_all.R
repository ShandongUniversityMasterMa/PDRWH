#run PDRWH in five different cohorts

getwd()
source("PDRWH_software.R")

cancer <- c("BRCA","KIRC","LIHC","GBM","STAD")

for (i in 1:length(cancer)) {
  cancer_type <- cancer[i]
  print("----------------------------------------------------")
  print(paste0("Run PDRWH for cancer type:",cancer_type,"..."))
  cat("\n\n")
  
  snv <- paste0("../data/",cancer_type,"/",cancer_type,"_mc3_gene_level.txt")
  outfile_dir <- paste0("../out/",cancer_type,"/")
  PDRWHscore(cancer_type,snv,"../data/STRINGv10.txt",outfile_dir)
  resultTransform(cancer_type,outfile_dir,outfile_dir)
  # PDRWH_cohort_score(outfile_dir)
  
  print("----------------------------------------------------")
  print(paste0("PDRWH for ",cancer_type,":  Achieved!"))
  print("----------------------------------------------------")
  cat("\n\n")
  if (i==length(cancer)) {
    print("----------------------------------------------------")
    print(paste0("Congratrulations! All achieved!"))
    print("----------------------------------------------------")
  }
}



