#' Feature plot for Oryza Root
#' @param seu Seurat object
#' @param res Resolution used for clustering (numeric)
#' @param doublet.rate Doublet rate estimated (numeric)
#' @param dir_to_color_scheme Directory to color scheme
#' @param save Output directory for the feature plot
#' @import Seurat
#' @import ggplot2
#' @import gridExtra

feature_plot_os <- function(seu,res,doublet.rate,dir_to_color_scheme,save){
  load(dir_to_color_scheme)
  color <- celltypepalette[sort(match(unique(seu$celltype.plate.ID.P),celltypeorder))];
  p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype.plate.ID.P", cols = color)+ggtitle("TRAP Seq Annotation")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  color <- celltypepalette[sort(match(unique(seu$celltype.field.ID.P),celltypeorder))];
  p2 <- DimPlot(seu, reduction = "umap", group.by = "celltype.field.ID.P", cols = color)+ggtitle(" ")+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p3 <- DimPlot(seu, reduction = "umap", group.by = "timezone.ID.P", order = c("Maturation2","Maturation1","Elongation","Meristem"), cols = c("#ebf8e3", "#51c8bd", "#009ac8", "#005fa8"))+ggtitle("RNA Seq Annotation")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  seu@meta.data$log_nCount_RNA <- log10(seu@meta.data$nCount_RNA);
  p4 <- FeaturePlot(seu, reduction = "umap", features = c("log_nCount_RNA"))+ggtitle("log10 UMI Counts")+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  seu@meta.data$log_nFeature_RNA <- log10(seu@meta.data$nFeature_RNA);
  p5 <- FeaturePlot(seu, reduction = "umap", features = c("log_nFeature_RNA"))+ggtitle("log10 Number of Genes")+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p6 <- FeaturePlot(seu, reduction = "umap", features = c("percent.mt"))+ggtitle("Percent Mitochondrial")+theme(plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p7 <- DimPlot(seu, reduction = "umap", group.by = colnames(seu@meta.data)[grep("DF.classifications",colnames(seu@meta.data))], order = c("Doublet","Singlet"),cols = c("#cccccc", "#ff4040"))+ggtitle(paste0("Doublet Rate ",doublet.rate," %"))+theme(plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  p8 <- DimPlot(seu, reduction = "umap", label = TRUE)+ggtitle(paste0("Clustering (resolution ",res,")"))+theme(legend.position="none",plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
  gl <- lapply(list(p1, p2, p3, p4, p5, p6, p7, p8), ggplotGrob)
  gwidth <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
  gl <- lapply(gl, "[[<-", "widths", value = gwidth)

  png(save, width = 22, height = 32, units = 'in', res=300)
  print(
    gridExtra::grid.arrange(grobs=gl, ncol=2)
  )
  dev.off()
}
