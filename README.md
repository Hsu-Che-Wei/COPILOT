# COPILOT (Cell preprOcessing PIpeline kaLlistO busTools) 

 🚨 **Update!! (Sept 2025):** 
 
The easiest way to run COPILOT (unfiltered matrix -> filtered matrix) now is to follow the notebook in the folder [COPILOT_2025_update]("./COPILOT_2025_update").

---

Single cell RNA-seq preprocessing tool for gene-by-cell matrices of UMI counts. The utility and summary file offered by this tool is directly comparable to CellRanger 3.1, the tool developed by 10X Genomics. It is recommended to use the raw spliced and unpliced counts matrices produced by scKB pipeline as the input of COPILOT.

## Citation

Shahan, R., Hsu, C.-W., Nolan, T.M., Cole, B.J., Taylor, I.W., Greenstreet, L., Zhang, S., Afanassiev, A., Vlot, A.H.C., Schiebinger, G., Benfey, P.N., Ohler, U., 2022. A single-cell Arabidopsis root atlas reveals developmental trajectories in wild-type and cell identity mutants. Developmental Cell S1534580722000338. https://doi.org/10.1016/j.devcel.2022.01.008

Hsu CW, Shahan R, Nolan TM, Benfey PN, Ohler U. Protocol for fast scRNA-seq raw data processing using scKB and non-arbitrary quality control with COPILOT. STAR Protoc. 2022 Dec 16;3(4):101729. doi: 10.1016/j.xpro.2022.101729. Epub 2022 Sep 30. PMID: 36181683; PMCID: PMC9530667.


## Advantages over CellRanger 3.1

1. Since scKB and COPILOT incorporates kallisto, which is a read aligner using peudoalignment, they are about ~ 30 times faster than CellRanger 3.1 in processing time. 

2. If user wants to calculate RNA velocity or infer trajectory based on splicing dynamics in downstream analysis, then COPILOT along with scKB are the right choice, since together they produce spliced and unspliced quality-filtered count matrices. 

3. If user wants to use specific sets of genes that represent the signal of noise or low quality cells for cell-filtering, COPILOT offers such function.

## Dependencies

R (>= 3.6.1)

R packages: Seurat (>= 3.1.4), Matrix, scales, rjson, gridExtra, ggplot2, DoubletFinder, DropletUtils, R2HTML

Python modules: umap-learn (required), leidenalg (optional)

## Installation (in R/RStudio)
   ```
devtools::install_github('Hsu-Che-Wei/COPILOT')
   ```
   
## Docker Image
   ```
#Pull image from docker hub   
docker pull cheweihsu/copilot:latest

#Run image
docker run -it -d cheweihsu/copilot:latest

#Check running image ID
docker ps -a

#Enter the container
docker exec -it {image ID} bash

#In the container, activate conda environment "copilot"
conda activate copilot

#Then you can run copilot in R!
   ```

## Tutorial

[0-COPILOT_tutorial](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/0-COPILOT_tutorial.ipynb)

[0-COPILOT_tutorial_toy_data](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/0-COPILOT_tutorial_toy_data.ipynb)

## COPILOT summary file

![COPILOT_Summary](/images/COPILOT_Summary.png)

Please check out examples of COPILOT summary file [here](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/COPILOT_summary_file/COPILOT_summary_files.pdf), or check out the folder named "COPILOT_summary_file".


## Usage

copilot( sample.name, spliced.mtx = NULL, unspliced.mtx = NULL, total.mtx = NULL, filtered.mtx.output.dir = NULL, species.name = "Not Provided", transcriptome.name = "Not Provided", sample.stats = NULL, mt.pattern, mt.threshold = 5, cp.pattern = NULL, top.percent = 1, filtering.ratio = 1, estimate.doublet.rate = TRUE, doublet.rate = NULL, remove.doublet = TRUE, do.seurat = TRUE, do.annotation = FALSE, unwanted.genes = NULL, HVG = FALSE, HVGN = 200, dir_to_bulk = NULL, dir_to_color_scheme = NULL, clustering_alg = 3, res = 0.5, min.UMI.low.quality = 100, min.UMI.high.quality = 300, legend.position = c(0.8, 0.8) )

## Arguments

> + ```sample.name```: User defined sample name (character), which should be the same as the name of directory that contains spliced and unspliced matrices if you are following scKB pipeline to produce raw counts matrices.
> + ```spliced.mtx```: Gene by cell matrix of spliced counts, which should have column and row names, Default is NULL.
> + ```unspliced.mtx```: Gene by cell matrix of unspliced counts, which should have column and row names. Default is NULL.
> + ```total.mtx```: Gene by cell matrix of total counts, which should have column and row names. Default is NULL.
> + ```filtered.mtx.output.dir```: Output directory for quality filtered matrices. Default is NULL.
> + ```species.name```: Species name (character). Default is "Not Provided".
> + ```transcriptome.name```: Name of transcriptome annotation file. (e.g. TAIR10 for _Arabidopsis_). Default is "Not Provided".
> + ```sample.stats```: Meta data of the sample in data.frame format. Default is NULL.
> + ```mt.pattern```: Pattern of mitochondrial gene names/ids (character; e.g. "ATMG") or list of mitochondrial genes (character vector). This argument is required to run copilot.
> + ```mt.threshold```: Threshold of mitochondrial expression percentage. Cell would be treated as dying cell if it has mitochodrial expression percentage higher than this threshold (numeric). Default is 5.
> + ```cp.pattern```: Pattern of chloroplast gene names/ids (character; e.g. "ATCG") or list of chloroplast genes (character vector). Default is NULL.
> + ```top.percent```: Percentage of cells that contain high numer of UMIs filtered (numeric). Default is 1.
> + ```filtering.ratio```: Metric that controls the stringency of cell filtering (lenient: 1; strict:0; moderate: 0 < filtering.ratio < 1; numeric). Default is 1.
> + ```estimate.doublet.rate```: Whether or not to estimate doublet rate according to 10X Genomics' estimation (boolean). Default is TRUE.
> + ```doublet.rate```: User specified doublet rate (numeric). Default is NULL.
> + ```remove.doublet```: Whether or not to remove doublets after quality filtering of gene and cell (boolean). Default is TRUE.
> + ```do.seurat```: Whether or not to perform normalization, PCA, UMAP and clustering using Seurat and output a Seurat object (boolean). Default is TRUE.
> + ```do.annotation```: Whether or not to do annotation (boolean). COPILOT only supports annotation on root of _Arabidopsis thaliana_. Default is FALSE.
> + ```unwanted.genes```: Gene IDs/names of unwanted genes (character vector, e.g. cell cycle related genes, organelle genes ... etc). Default is NULL.
> + ```HVG```: Whether or not to select highly variable genes (boolean). Default is FALSE.
> + ```HVGN```: Number of highly variable genes selected (numeric). Defalut is 200.
> + ```dir_to_bulk```: Directory to reference expression profile for annotation. Default is NULL.
> + ```dir_to_color_scheme```: Directory to color scheme file for annotation. Default is NULL.
> + ```clustering_alg```: Algorithm for clustering (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm, which requires the leidenalg python; numeric). Default is 3.
> + ```res```: Resolution used for clustering (numeric). Default is 0.5.
> + ```min.UMI.low.quality```: Minimum UMIs for a barcode to be considered as cell (numeric). Default is 100.
> + ```min.UMI.high.quality```: Minimum UMIs for a cell to be considered as high quality cell (numeric). Default is 300.
> + ```legend.position```: x y position of the legend on UMI histogram plot (numeric vector of length 2). Default is c(0.8,0.8).

## Caution

It is recommended to use the raw spliced and unpliced counts matrices produced by scKB pipeline as the input of copilot, and run copilot in the working directory that contains folder with raw matrices inside. Note that the folder name should be the same as the argument "sample.name". Alternatively, user can directly feed the matrices to copilot in R by filling in the arguments "spliced.mtx", "unspliced.mtx" and "total.mtx". Note that if "total.mtx" is provided, then please let the "spliced.mtx", "unspliced.mtx" remain as NULL, and vise versa.

COPILOT relies on detecting the distribution of cells enriched in mitochondrial expression (putative dying cells) to identify the initial filtering threshold for low quality cells. Thus, even though user is confident that the library does not contain dying cells, it is still required for user to provide list of mitochondrial gene IDs or a character pattern that is shared among those gene IDs (See argument "mt.pattern"). Alternatively, user can provide genes that represent the signal of noise or low quality cells.

If the algorithm detects no putative dying cells, the initial filtering threshold will then be set as what user feed to argument "min.UMI.high.quality". It is recommended to play around with this argument until the number of high quality cell remains after filtering makes sense.

For users that have scRNA-seq data of _Arabidopsis thaliana_ root, and wish to perform correlation-based annotation, please make sure that all necessary files are downloaded from folder named "supp_data".  

## Application 

One can find here codes demonstrating how some typical downstream analysis of scRNA-seq were performed on COPILOT-preprocessed data. The codes also demonstrate how analysis in [Shahan & Hsu et al. 2020](https://doi.org/10.1101/2020.06.29.178863) were done.

[1-Basic_Seurat_Analysis_&_Correlation_Based_Annotation](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/1-Correlation_Based_Annotation.ipynb)

[2-Integration](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/2-Integration.ipynb)

[3-1-novoSpaRc](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/3-1-novoSpaRc.ipynb)

[3-2-novoSpaRc_to_Atlas](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/3-2-novoSpaRc_to_Atlas.ipynb)

[4-1-Prepare_SEMITONES_Input](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/4-1-Prepare_SEMITONES_Input.ipynb)

[4-2-Prepare_SEMITONES_Input_2](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/4-2-Prepare_SEMITONES_Input_2.ipynb)

[4-3-Marker_Annotation_By_SEMITONES](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/4-3-Marker_Annotation_By_SEMITONES.ipynb)

[5-ICI_Computation](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/5-ICI_Computation.ipynb)

[6-CytoTRACE](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/6-CytoTRACE.ipynb)

[7-Combine_and_Finalize_Annotation](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/7-Combine_and_Finalize_Annotation.ipynb)

[8-Extract_Tissue_Lineage](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/8-Extract_Tissue_Lineage.ipynb)

[9-scVelo_Latent_Time](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/9-scVelo_Latent_Time.ipynb)

[10-Consensus_Time](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/10-Consensus_Time.ipynb)

[11_Label_Transfer](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/11_Label_Transfer.ipynb)

[12-1-Atlas_FindMarkers](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/12-1-Atlas_FindMarkers.ipynb)

[12-2-Endo_Cortex_Heatmap](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/12-2-Endo_Cortex_Heatmap.ipynb)

[12-3-WT_SCR_SHR_Heatmap](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/12-3-WT_SCR_SHR_Heatmap.ipynb)

[13-Differential_Abundance_Analysis](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/13-Differential_Abundance_Analysis.ipynb)

[14-Ploidy_Annotation](https://github.com/Hsu-Che-Wei/COPILOT/tree/master/jupyter_notebook/14-Ploidy_Annotation.ipynb)







