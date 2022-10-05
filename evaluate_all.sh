#!/bin/env Rscript

# This bash script can be used to evaluate the results of PDRWH method in five different cohorts.
# The personalized prediction results must be saved in the following order:
#########################################################
#
#		out/
#
#########################################################

#########################################################
#
#	Evaluation on BRCA, KIRC, LIHC, GBM, STAD
#
#########################################################

cd ./src/

Rscript P_R_F1.R

Rscript Hypergeometric_test.R

