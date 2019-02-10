# Schizophrenia DGE ShinyApp
Scripts for visualization of Schizophrenia CellType DGE.

### Usage
# Run ShinyApp loader: 
**R CMD BATCH --vanilla ToRun.R**
   
    1) Load the expression
    2) Load the demographic
    3) Have fun!

### Details
The app will ask:
- **Load Expression data** coming from *DGE analysis* and based on *Expression cleaned for covariates and surrogate variables*
- **Load demographic data** which contains *predictors + covariates + surrogate variables*
- Searching for the **Gene** of interest

The app will show:
- **Violin plot and Boxplot** with statistics. 
