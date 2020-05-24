# COPILOT
Single cell RNA-seq preprocessing tool for gene-by-cell matrices of UMI counts. It is recommended to use the raw spliced and unpliced counts matrices produced by scKB pipeline as the input of copilot.

## Dependencies

R (>= 3.6.1)

R packages: Seurat, Matrix, rjson, gridExtra, ggplot2, DoubletFinder, DropletUtils, R2HTML

umap-learn #installed by "pip install umap-learn"

leidenalg #optional, required if Leiden clustering is chosen

