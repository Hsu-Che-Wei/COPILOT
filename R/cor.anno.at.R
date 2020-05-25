#' Correlation-based annotation for Arabidopsis Root
#' @param seu Seurat object
#' @param dir_to_bulk Directory to reference expression profile
#' @param unwanted.genes Unwanted genes to exclude from analysis (character vector)
#' @param clustering_alg Algorithm for clustering (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm, which requires the leidenalg python; numeric)
#' @param res Resolution used for Leiden clustering (numeric)
#' @param mt.pattern Pattern of mitochondrial gene names/ids (character) ex. "ATMG" or list of mitochondrial genes (character vector)
#' @param cp.pattern Pattern of chloroplast gene names/ids (character) ex. "ATCG" or list of chloroplast genes (character vector)
#' @import Seurat

cor.anno.at <- function (seu, dir_to_bulk, unwanted.genes, clustering_alg, res, mt.pattern, cp.pattern) {

  library(Seurat)

  load(file=dir_to_bulk)

  #### Bulk RNA-seq ####
  rc <- as.matrix(seu@assays$SCT@data)
  time <- Reduce(merge.rownames, list(time,rc))
  celltype <- Reduce(merge.rownames, list(celltype,rc))

  time_label=c("Elongation", "Maturation", "Meristem")
  celltype_label=c("phloem & companion cells", "developing cortex", "hair cells", "matured cortex",
                   "matured endodermis", "non-hair cells", "columella", "phloem pole pericycle",
                   "matured xylem pole", "protophloem & metaphloem","developing xylem", "endodermis & QC cells", "LRC & non-hair cells","QC cells")

  time_stat <- suppressWarnings(sapply(4:ncol(time), function(i) sapply(1:3, function(j) cor.test(time[,i],time[,j],method = "spearman")[c(3,4)])))
  time_cor <- time_stat[seq(2,nrow(time_stat),2),]
  time_pvalue <- time_stat[seq(1,nrow(time_stat)-1,2),]
  time_max <- sapply(1:(ncol(time)-3), function(i) max(as.numeric(time_cor[,i])))
  time_ident <- sapply(1:(ncol(time)-3), function(i) time_label[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  time_maxp <- sapply(1:(ncol(time)-3), function(i) as.numeric(time_pvalue[,i])[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  names(time_max) <- time_ident

  celltype_stat <- suppressWarnings(sapply(15:ncol(celltype), function(i) sapply(1:14, function(j) cor.test(celltype[,i],celltype[,j],method = "spearman")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-14), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-14), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-14), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$celltype.ID <- celltype_ident
  seu@meta.data$timezone.ID <- time_ident
  seu@meta.data$celltype.cor <- celltype_max
  seu@meta.data$timezone.cor <- time_max
  seu@meta.data$celltype.pvalue <- celltype_maxp
  seu@meta.data$timezone.pvalue <- time_maxp

  time_stat <- suppressWarnings(sapply(4:ncol(time), function(i) sapply(1:3, function(j) cor.test(time[,i],time[,j],method = "pearson")[c(3,4)])))
  time_cor <- time_stat[seq(2,nrow(time_stat),2),]
  time_pvalue <- time_stat[seq(1,nrow(time_stat)-1,2),]
  time_max <- sapply(1:(ncol(time)-3), function(i) max(as.numeric(time_cor[,i])))
  time_ident <- sapply(1:(ncol(time)-3), function(i) time_label[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  time_maxp <- sapply(1:(ncol(time)-3), function(i) as.numeric(time_pvalue[,i])[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  names(time_max) <- time_ident

  celltype_stat <- suppressWarnings(sapply(15:ncol(celltype), function(i) sapply(1:14, function(j) cor.test(celltype[,i],celltype[,j],method = "pearson")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-14), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-14), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-14), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$celltype.ID.P <- celltype_ident
  seu@meta.data$timezone.ID.P <- time_ident
  seu@meta.data$celltype.cor.P <- celltype_max
  seu@meta.data$timezone.cor.P <- time_max
  seu@meta.data$celltype.pvalue.P <- celltype_maxp
  seu@meta.data$timezone.pvalue.P <- time_maxp

  #### Micro Array ####
  time <- Reduce(merge.rownames, list(Long,rc))
  celltype <- Reduce(merge.rownames, list(Rad,rc))

  time_label=c("Columella", "Meri-1", "Meri-2", "Meri-3", "Meri-4", "Meri-5", "Meri-6", "Elong-7", "Elong-8", "Mat-9", "Mat-10", "Mat-11", "Mat-12")
  celltype_label=c("QC", "Hair Cell", "Cortex", "Non-Hair Cell", "Xylem Pole Pericycle", "LRC",
                   "Columella", "Phloem Pole Pericycle", "Mat.Xylem", "Meri.Xylem", "Phloem [S32]", "Endodermis", "Phloem [SUC2]")

  time_stat <- suppressWarnings(sapply(14:ncol(time), function(i) sapply(1:13, function(j) cor.test(time[,i],time[,j],method = "spearman")[c(3,4)])))
  time_cor <- time_stat[seq(2,nrow(time_stat),2),]
  time_pvalue <- time_stat[seq(1,nrow(time_stat)-1,2),]
  time_max <- sapply(1:(ncol(time)-13), function(i) max(as.numeric(time_cor[,i])))
  time_ident <- sapply(1:(ncol(time)-13), function(i) time_label[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  time_maxp <- sapply(1:(ncol(time)-13), function(i) as.numeric(time_pvalue[,i])[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  names(time_max) <- time_ident

  celltype_stat <- suppressWarnings(sapply(14:ncol(celltype), function(i) sapply(1:13, function(j) cor.test(celltype[,i],celltype[,j],method = "spearman")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-13), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-13), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-13), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$Rad.ID <- celltype_ident
  seu@meta.data$Long.ID <- time_ident
  seu@meta.data$Rad.cor <- celltype_max
  seu@meta.data$Long.cor <- time_max
  seu@meta.data$Rad.pvalue <- celltype_maxp
  seu@meta.data$Long.pvalue <- time_maxp

  time_stat <- suppressWarnings(sapply(14:ncol(time), function(i) sapply(1:13, function(j) cor.test(time[,i],time[,j],method = "pearson")[c(3,4)])))
  time_cor <- time_stat[seq(2,nrow(time_stat),2),]
  time_pvalue <- time_stat[seq(1,nrow(time_stat)-1,2),]
  time_max <- sapply(1:(ncol(time)-13), function(i) max(as.numeric(time_cor[,i])))
  time_ident <- sapply(1:(ncol(time)-13), function(i) time_label[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  time_maxp <- sapply(1:(ncol(time)-13), function(i) as.numeric(time_pvalue[,i])[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  names(time_max) <- time_ident

  celltype_stat <- suppressWarnings(sapply(14:ncol(celltype), function(i) sapply(1:13, function(j) cor.test(celltype[,i],celltype[,j],method = "pearson")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-13), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-13), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-13), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$Rad.ID.P <- celltype_ident
  seu@meta.data$Long.ID.P <- time_ident
  seu@meta.data$Rad.cor.P <- celltype_max
  seu@meta.data$Long.cor.P <- time_max
  seu@meta.data$Rad.pvalue.P <- celltype_maxp
  seu@meta.data$Long.pvalue.P <- time_maxp

  # Run PCA, UMAP, Clustering
  use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)),grep(paste(cp.pattern, collapse = "|"),rownames(seu)),sort(match(unwanted.genes, rownames(seu))))]
  seu <- RunPCA(seu, verbose = FALSE, approx = FALSE, npcs = 50, features=use.genes)
  suppressMessages(suppressWarnings(
    seu <- RunUMAP(seu, reduction = "pca", dims = 1:50, umap.method = "umap-learn", metric = "correlation")
  ))
  suppressMessages(suppressWarnings(
    seu <- FindNeighbors(seu, reduction = "pca",dims = 1:50)
  ))
  suppressMessages(suppressWarnings(
    seu <- FindClusters(seu, resolution = res, algorithm = clustering_alg)
  ))

  return(seu)
}
