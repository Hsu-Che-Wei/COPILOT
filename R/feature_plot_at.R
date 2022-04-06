#' Feature plot for Arabidopsis Root
#' @param seu Seurat object
#' @param res Resolution used for clustering (numeric)
#' @param doublet.rate Doublet rate estimated (numeric)
#' @param dir_to_color_scheme Directory to color scheme
#' @param save Output directory for the feature plot
#' @import Seurat
#' @import ggplot2
#' @import scales
#' @import egg

feature_plot_at <- function(seu,res,doublet.rate,dir_to_color_scheme,save){
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(egg)
  load(dir_to_color_scheme)
  color <- celltypepalette[sort(match(unique(seu$celltype.ID.P),celltypeorder))];
  p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype.ID.P", cols = color)+ggtitle("RNA Seq Annotation")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p2 <- DimPlot(seu, reduction = "umap", group.by = "timezone.ID.P", order = c("Maturation","Elongation","Meristem"),cols = c("#DCEDC8", "#42B3D5", "#1A237E"))+ggtitle(" ")+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  color.rad <- radpalette[sort(match(unique(seu$Rad.ID.P),radorder))];
  p3 <- DimPlot(seu, reduction = "umap", group.by = "Rad.ID.P", cols = color.rad)+ggtitle("Microarray Annotation")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  seu@meta.data$Long.ID.P <- factor(seu@meta.data$Long.ID.P, levels = longorder[sort(match(unique(seu$Long.ID.P),longorder))])
  color.long <- longpalette[sort(match(unique(seu$Long.ID.P),longorder))];
  p4 <- DimPlot(seu, reduction = "umap", group.by = "Long.ID.P",cols = color.long)+ggtitle(" ")+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  seu@meta.data$log_nCount_RNA <- log10(seu@meta.data$nCount_RNA);
  p5 <- FeaturePlot(seu, reduction = "umap", features = c("log_nCount_RNA"))+ggtitle("log10 UMI Counts")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  seu@meta.data$log_nFeature_RNA <- log10(seu@meta.data$nFeature_RNA);
  p6 <- FeaturePlot(seu, reduction = "umap", features = c("log_nFeature_RNA"))+ggtitle("log10 Number of Genes")+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p7 <- FeaturePlot(seu, reduction = "umap", features = c("percent.mt"))+ggtitle("Percent Mitochondrial")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p8 <- DimPlot(seu, reduction = "umap", group.by = colnames(seu@meta.data)[grep("DF.classifications",colnames(seu@meta.data))], order = c("Doublet","Singlet"),cols = c("#cccccc", "#ff4040"))+ggtitle(paste0("Doublet Rate ",doublet.rate," %"))+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p9 <- DimPlot(seu, reduction = "umap", label = TRUE)+ggtitle(paste0("Clustering (resolution ",res,")"))+theme(legend.position="none",plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  gl <- lapply(list(p1, p2, p3, p4, p5, p6, p7, p8, p9), ggplotGrob)
  gwidth <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
  gl <- lapply(gl, "[[<-", "widths", value = gwidth)

  png(save, width = 22, height = 40, units = 'in', res=300)
  print(
    #gridExtra::grid.arrange(grobs=gl, ncol=2)
    ggarrange(gl, ncol=2)
  )
  dev.off()
}
