# Schizophrenia DGE
Scripts for RNA seq analysis.

### Usage
# 1) Run DGE analysis: 
**rmarkdown::render("NeuN_Dge_Analysis.R")**

**rmarkdown::render("OLIG2_Dge_Analysis.R")**

### Details
The script involves 3 steps
- **Data transformation** using *log2(CPM)*
- **QC** of tranformed counts
- **Modeling** of transformed counts based on *limma*

# 2) After DGE analysis use:
**R CMD BATCH --vanilla CrossValidation_NeuN.R**

**R CMD BATCH --vanilla CrossValidation_OLIG2.R**

### Details
The script involves 2 steps
- **Leave-One-Out CV** for *limma* based on 200 bootstraps
- **Permutation CV** for *limma* based on 200 permutations
