#' Correlation-based annotation for Oryza Root
#' @param seu Seurat object
#' @param dir_to_bulk Directory to reference expression profile
#' @param unwanted.genes Unwanted genes to exclude from analysis (character vector)
#' @param clustering_alg Algorithm for clustering (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm, which requires the leidenalg python; numeric)
#' @param res Resolution used for Leiden clustering (numeric)
#' @param mt.pattern Pattern of mitochondrial gene names/ids (character) ex. "ATMG" or list of mitochondrial genes (character vector)
#' @param cp.pattern Pattern of chloroplast gene names/ids (character) ex. "ATCG" or list of chloroplast genes (character vector)
#' @import Seurat

cor.anno.os <- function (seu, dir_to_bulk, unwanted.genes, clustering_alg, res, mt.pattern, cp.pattern) {

  library(Seurat)

  load(file=dir_to_bulk)

  #### Bulk RNA-seq ####
  rc <- as.matrix(seu@assays$SCT@data)
  use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)),grep(paste(cp.pattern, collapse = "|"),rownames(seu)),sort(match(unwanted.genes, rownames(seu))))]
  rc <- rc[use.genes,]
  time <- Reduce(merge.rownames, list(time_MSU,rc))
  celltype <- Reduce(merge.rownames, list(celltype_MSU,rc))
  celltype <- celltype[,-c(3,6:12)]

  time_label=c("Meristem", "Elongation", "Maturation1", "Maturation2")
  celltype_label=c("Meri.Cortex.Plate", "QC.Plate",
                   "Meri.Endodermis.Plate", "Stele.LRP.Plate")


  time_stat <- suppressWarnings(sapply(5:ncol(time), function(i) sapply(1:4, function(j) cor.test(time[,i],time[,j],method = "spearman")[c(3,4)])))
  time_cor <- time_stat[seq(2,nrow(time_stat),2),]
  time_pvalue <- time_stat[seq(1,nrow(time_stat)-1,2),]
  time_max <- sapply(1:(ncol(time)-4), function(i) max(as.numeric(time_cor[,i])))
  time_ident <- sapply(1:(ncol(time)-4), function(i) time_label[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  time_maxp <- sapply(1:(ncol(time)-4), function(i) as.numeric(time_pvalue[,i])[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  names(time_max) <- time_ident

  celltype_stat <- suppressWarnings(sapply(5:ncol(celltype), function(i) sapply(1:4, function(j) cor.test(celltype[,i],celltype[,j],method = "spearman")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-4), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-4), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-4), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$celltype.plate.ID <- celltype_ident
  seu@meta.data$timezone.ID <- time_ident
  seu@meta.data$celltype.plate.cor <- celltype_max
  seu@meta.data$timezone.cor <- time_max
  seu@meta.data$celltype.plate.pvalue <- celltype_maxp
  seu@meta.data$timezone.pvalue <- time_maxp

  celltype <- Reduce(merge.rownames, list(celltype_MSU,rc))
  celltype <- celltype[,-c(1:5,11)]

  celltype_label=c( "Root.Hair.Field", "Pericycle.Field", "Endo.Exodermis.Field",
                    "Phloem.Pole.Field","QC.Field", "Stele.LRP.Field")

  celltype_stat <- suppressWarnings(sapply(7:ncol(celltype), function(i) sapply(1:6, function(j) cor.test(celltype[,i],celltype[,j],method = "spearman")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-6), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-6), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-6), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$celltype.field.ID <- celltype_ident
  seu@meta.data$celltype.field.cor <- celltype_max
  seu@meta.data$celltype.field.pvalue <- celltype_maxp

  #### Bulk RNA-seq Pearson ####

  celltype <- Reduce(merge.rownames, list(celltype_MSU,rc))
  celltype <- celltype[,-c(3,6:12)]

  time_label=c("Meristem", "Elongation", "Maturation1", "Maturation2")
  celltype_label=c("Meri.Cortex.Plate", "QC.Plate",
                   "Meri.Endodermis.Plate", "Stele.LRP.Plate")


  time_stat <- suppressWarnings(sapply(5:ncol(time), function(i) sapply(1:4, function(j) cor.test(time[,i],time[,j],method = "pearson")[c(3,4)])))
  time_cor <- time_stat[seq(2,nrow(time_stat),2),]
  time_pvalue <- time_stat[seq(1,nrow(time_stat)-1,2),]
  time_max <- sapply(1:(ncol(time)-4), function(i) max(as.numeric(time_cor[,i])))
  time_ident <- sapply(1:(ncol(time)-4), function(i) time_label[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  time_maxp <- sapply(1:(ncol(time)-4), function(i) as.numeric(time_pvalue[,i])[which(as.numeric(time_cor[,i])==max(as.numeric(time_cor[,i])))])
  names(time_max) <- time_ident

  celltype_stat <- suppressWarnings(sapply(5:ncol(celltype), function(i) sapply(1:4, function(j) cor.test(celltype[,i],celltype[,j],method = "pearson")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-4), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-4), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-4), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$celltype.plate.ID.P <- celltype_ident
  seu@meta.data$timezone.ID.P <- time_ident
  seu@meta.data$celltype.plate.cor.P <- celltype_max
  seu@meta.data$timezone.cor.P <- time_max
  seu@meta.data$celltype.plate.pvalue.P <- celltype_maxp
  seu@meta.data$timezone.pvalue.P <- time_maxp

  celltype <- Reduce(merge.rownames, list(celltype_MSU,rc))
  celltype <- celltype[,-c(1:5,11)]

  celltype_label=c( "Root.Hair.Field", "Pericycle.Field", "Endo.Exodermis.Field",
                    "Phloem.Pole.Field","QC.Field", "Stele.LRP.Field")

  celltype_stat <- suppressWarnings(sapply(7:ncol(celltype), function(i) sapply(1:6, function(j) cor.test(celltype[,i],celltype[,j],method = "pearson")[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-6), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-6), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-6), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident

  seu@meta.data$celltype.field.ID.P <- celltype_ident
  seu@meta.data$celltype.field.cor.P <- celltype_max
  seu@meta.data$celltype.field.pvalue.P <- celltype_maxp

  # Run PCA, UMAP, Clustering
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
