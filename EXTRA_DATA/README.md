# BrainSeq/CMC DEG
Scripts for RNA seq analysis.

### Usage
# 1) Run DEG analysis: 
**rmarkdown::render("BrainSeq_DEG_Analysis.R")**

**rmarkdown::render("CMC_DEG_Analysis.R")**

### Details
The script involves 3 steps
- **Data transformation** using *log2(CPM)*
- **QC** of tranformed counts
- **Modeling** of transformed counts based on *limma*

# 2) After DEG analysis enter the COMPARATIVE directory:
**R CMD BATCH --vanilla Comparative_Analysis.R**

### Details
The script will:
- **Compare statistics** for *NeuN* fold changes vs *CMC/BrainSeq* fold changes 
- **Compare statistics** for *OLIG2* fold changes vs *CMC/BrainSeq* fold changes 
- **Output Databases and Scatterplot**
