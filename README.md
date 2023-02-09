# PDRWH: prioritizing Personalized cancer Driver genes via Random Walk on Hypergraph model

### This is the original repository for the PDRWH paper. PDRWH requires R >= 4.2.1.

**Installing the dependencies**

```
R

install.packages(c("igraph","stats","readxl","ggplot2","RColorBrewer","reshape2","dplyr"))

```

## **Input**

There are two files for PDRWH: the Protein-Protein Interaction (PPI) network edges file and the somatic mutation data of tumor patients. 
Both files are located in data folder

1. The PPI network file:

The file is located at 'data/STRINGv10.txt'.

```
v1(node_gene_id) v2(node_gene_id) weight
ARF5	CFTR	0.15015015015015
ARF5	CYP51A1	0.215215215215215
ARF5	RALA	0.223223223223223

```

2. The somatic mutation data:

The file is located at 'data/**/**_mc3_gene_level.txt'.
'**' represents a certain cancer type.

For example, 'data/BRCA/BRCA_mc3_gene_level.txt':
```
sample	TCGA-3C-AAAU-01	TCGA-3C-AALI-01	TCGA-3C-AALJ-01
UBE2Q2	0	0	0
CHMP1B	0	0	0
PSMA2P1	0	0	0
SHQ1P1	0	0	0
CPHL1P	0	0	0

```

### Known general driver genes and tumor-specific driver genes

The two files are located in data folder. 

```
'data/general_driver_list.csv'
'data/Cancer_specific.xls'

```

These files are used to find genes predicted by PDRWH or other methods that are also present in the known driver list.


## **Run**

There are two bash scripts, one for running PDRWH algorithm, the other for evaluation.

1. Personalized prediction:

```
cd ./PDRWH-1.0/

sh execute_all.sh

```

2. Evaluation

```
cd ./PDRWH-1.0/

sh evaluation_all.sh

```


## **Outputs**

The 'execute_all.sh' file will output personalized prediction results of patients in the 'out/**/' folder.
'**' represents a certain cancer type. For example, BRCA:

```
'out/BRCA/BRCA.Rdata'
'out/BRCA/PDRWH.txt'

```

The 'evaluation_all.sh' file will output the average precision, recall and F1 score for the top N predicted driver genes by a known general reference driver gene list. 
The PDF files can be found in 'out/**_topN.pdf'. '**' represents a certain cancer type. For example, BRCA:

```
'out/BRCA_topN.pdf'

```

The 'evaluation_all.sh' file also show the result of hypergeometric test, which indicates that the predicted driver genesfor certainpatient is significantly enriched in known tumor-specific driver genes.

```
'out/Hypergeometric test.txt'

```

The picture will be showed in 'RStudio Plots' window.


