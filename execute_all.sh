#!/bin/env Rscript

# This bash script can be used to apply PDRWH method on five different cohorts.
# The required R library are igraph. 
# The personalized prediction results must be saved in a folder 
# with the name of the cancer type or the method in the following order:
#########################################################
#
#		out/BRCA/BRCA.Rdata
#
#               out/BRCA/PDRWH.txt
#
#########################################################

#########################################################
#
#	Apply PDRWH on BRCA, KIRC, LIHC, GBM, STAD
#
#########################################################

cd ./src

Rscript runPDRWH_all.R

