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

## Usage

copilot( sample.name, spliced.mtx = NULL, unspliced.mtx = NULL, total.mtx = NULL, filtered.mtx.output.dir = NULL, species.name = "Not Provided", transcriptome.name = "Not Provided", sample.stats = NULL, mt.pattern = NULL, mt.threshold = 5, cp.pattern = NULL, top.percent = 1, filtering.ratio = 1, estimate.doublet.rate = TRUE, doublet.rate = NULL, remove.doublet = TRUE, do.seurat = TRUE, do.annotation = FALSE, unwanted.genes = NULL, HVG = FALSE, HVGN = 200, dir_to_bulk = NULL, clustering_alg = 3, res = 0.5, min.UMI.low.quality = 100, min.UMI.high.quality = 300, legend.position = c(0.8, 0.8) )

## Arguments

sample.name: User defined sample name (character), which should be the same as the name of directory that contains spliced and spliced matrices if you are following scKB pipeline to produce raw counts matrices.

spliced.mtx: Gene by cell matrix of spliced counts, which should have column and row names, Default is NULL.

unspliced.mtx: Gene by cell matrix of unspliced counts, which should have column and row names. Default is NULL.

total.mtx: Gene by cell matrix of total counts, which should have column and row names. Default is NULL.

filtered.mtx.output.dir: Output directory for quality filtered matrices. Default is NULL.

species.name: Species name (character). Default is "Not Provided".

transcriptome.name: Name of transcriptome annotation file. (e.g. TAIR10 for Arabidopsis). Default is "Not Provided".

sample.stats: Meta data of the sample in data.frame format. Default is NULL.

mt.pattern: Pattern of mitochondrial gene names/ids (character; e.g. "ATMG") or list of mitochondrial genes (character vector). Default is NULL, however this argument is required to run copilot.

mt.threshold: Threshold of mitochondrial expression percentage. Cell would be treated as dying cell if it has mitochodrial expression percentage higherthan this threshold (numeric). Default is 5.

cp.pattern: Pattern of chloroplast gene names/ids (character; e.g. "ATCG") or list of chloroplast genes (character vector). Default is NULL.

top.percent: Percentage of cells that contain high numer of UMIs filtered (numeric). Default is 1.

filtering.ratio: Metric that controls the stringency of cell filtering (lenient: 1; strict:0; moderate: 0 < filtering.ratio < 1; numeric). Default is 1.

estimate.doublet.rate: Whether or not to estimate doublet rate according to 10X Genomics' estimation (boolean). Default is TRUE.

doublet.rate: User specified doublet rate (numeric). Default is NULL.

remove.doublet: Whether or not to remove doublets after quality filtering of gene and cell (boolean). Default is TRUE.

do.seurat: Whether or not to perform normalization, PCA, UMAP and clustering using Seurat and output a Seurat object (boolean). Default is TRUE.

do.annotation: Whether or not to do annotation (boolean). COPILOT only supports annotation on root of Arabidopsis thaliana and Oryza sativa. Default is FALSE.

unwanted.genes: Gene IDs/names of unwanted genes (character vector, e.g. cell cycle related genes, organelle genes ... etc). Default is NULL.

HVG: Whether or not to select highly variable genes (boolean). Default is FALSE.

HVGN: Number of highly variable genes selected (numeric). Defalut is 200.

dir_to_bulk: Directory to reference expression profile for annotation. Default is NULL.

clustering_alg: Algorithm for clustering (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm, which requires the leidenalg python; numeric). Default is 3.

res: Resolution used for clustering (numeric). Default is 0.5.

min.UMI.low.quality: Minimum UMIs for a barcode to be considered as cell (numeric). Default is 100.

min.UMI.high.quality: Minimum UMIs for a cell to be considered as high quality cell (numeric). Default is 300.

legend.position: x y position of the legend on UMI histogram plot (numeric vector of length 2). Default is c(0.8,0.8).

## Caution

COPILOT relies on detecting the distribution of cells enriched in mitochondrial expression (putative dying cells) to identify the initial filtering threshold for low quality cells. Thus, even though user is confident that the library does not contain dying cells, it is still required for user to provide list of mitochondrial gene IDs or a character pattern that is shared among those gene IDs (See argument "mt.pattern"). Alternatively, user can provide genes that represent the signal of noise or low quality cells.

If the algorithm detects no putative dying cells, the initial filtering threshold will then be set as what user feed to argument "min.UMI.high.quality". It is recommended to play around with this argument until the number of high quality cell remains after filtering makes sense.
