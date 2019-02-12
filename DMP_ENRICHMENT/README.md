# DMP enrichment for scDEGs
The directories contain:
       
       1) Input data from differential expression and methylation for NeuN and OLIG2.
       2) Codes to set up the enrichment analysis.
       
### Usage
# 1) Run DEG analysis: 
**R CMD BATCH --vanilla NeuN_Methylation_Enrichment.R**

**R CMD BATCH --vanilla OLIG2_Methylation_Enrichment.R**

### Details
The script involves 2 steps
- **Data formatting**
- **Enrichment test** based on *Fisher's exact test*
- **Output** as heatmap containing Odd Ratios and FDR from the *Fisher's exact test*
