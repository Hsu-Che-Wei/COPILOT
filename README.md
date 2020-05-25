# COPILOT
Single cell RNA-seq preprocessing tool for gene-by-cell matrices of UMI counts. It is recommended to use the raw spliced and unpliced counts matrices produced by scKB pipeline as the input of copilot.

## Dependencies

R (>= 3.6.1)

R packages: Seurat (>= 3.1.4), Matrix, rjson, gridExtra, ggplot2, DoubletFinder, DropletUtils, R2HTML

Python modules: umap-learn, leidenalg (optional)

## Installation (in R/RStudio)
   ```
devtools::install_github('Hsu-Che-Wei/COPILOT')
   ```

## Tutorial

Please check out the jupyter notebook named "COPILOT_tutorial.ipynb" or "COPILOT_tutorial.html" for convenience.

## COPILOT summary file

Please check an example of COPILOT summary file in the folder named "COPILOT_summary_file".
