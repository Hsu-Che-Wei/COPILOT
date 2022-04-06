#' Feature plot
#' @param seu Seurat object
#' @param res Resolution used for clustering (numeric)
#' @param doublet.rate Doublet rate estimated (numeric)
#' @param save Output directory for the feature plot
#' @import Seurat
#' @import ggplot2
#' @import scales
#' @import egg

feature_plot <- function(seu,res,doublet.rate,save){
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(egg)
  seu@meta.data$log_nCount_RNA <- log10(seu@meta.data$nCount_RNA);
  p1 <- FeaturePlot(seu, reduction = "umap", features = c("log_nCount_RNA"))+ggtitle("log10 UMI Counts")+theme(plot.margin = unit(c(2.4,0.8,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  seu@meta.data$log_nFeature_RNA <- log10(seu@meta.data$nFeature_RNA);
  p2 <- FeaturePlot(seu, reduction = "umap", features = c("log_nFeature_RNA"))+ggtitle("log10 Number of Genes")+theme(plot.margin = unit(c(2.4,0,0.8,0.8), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p3 <- FeaturePlot(seu, reduction = "umap", features = c("percent.mt"))+ggtitle("Percent Mitochondrial")+theme(plot.margin = unit(c(2.4,0.8,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p4 <- DimPlot(seu, reduction = "umap", group.by = colnames(seu@meta.data)[grep("DF.classifications",colnames(seu@meta.data))], order = c("Doublet","Singlet"),cols = c("#cccccc", "#ff4040"))+ggtitle(paste0("Doublet Rate ",doublet.rate," %"))+theme(plot.margin = unit(c(2.4,0,0.8,0.8), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p5 <- DimPlot(seu, reduction = "umap", label = TRUE)+ggtitle(paste0("Clustering (resolution ",res,")"))+theme(legend.position="none",plot.margin = unit(c(2.4,0.8,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  gl <- lapply(list(p1, p2, p3, p4, p5), ggplotGrob)
  gwidth <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
  gl <- lapply(gl, "[[<-", "widths", value = gwidth)

  png(save, width = 22, height = 24, units = 'in', res=300)
  print(
    #gridExtra::grid.arrange(grobs=gl, ncol=2)
    ggarrange(gl, ncol=2)
  )
  dev.off()
}
