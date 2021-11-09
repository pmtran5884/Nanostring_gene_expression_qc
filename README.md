# Nanostring_gene_expression_qc
R code to perform QC and normalization of multiple batches of nCounter gene expression data

Nanostring RCC data must first be downloaded from the digital analyzer and then read into nSolver software or R package. Then output the results as a csv file to begin the steps below.

Step 0: Visualize raw data to detect batch effects and outliers
Step 1: Remove samples and genes with low counts
Step 2: background thresholding
Step 3: Normalize samples using Positive Control and Housekeeping Genes (remove low quality positive control and housekeeping genes from normalization estimate; use geometric mean for normalization)
Step 4: Visualize
Step 5: Batch Effect Removal (uses ComBat algorithm as implemented in "sva" package)
