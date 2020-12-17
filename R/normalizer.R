

normalizer <- function(sample.name, seu, HVG = FALSE, HVGN = 200){
  library(Seurat)


  # Whether to apply highly variable gene selection
  if (HVG) {
    vfn <- HVGN
  } else {
    vfn <- nrow(seu)
  }
  message("SCTransforming...")

  # Normalization using SCTransform
  suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, assay = "RNA", new.assay.name = "SCT", verbose = T) )
  #suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, assay = "spliced_RNA", new.assay.name = "spliced_SCT", verbose = FALSE))
  #suppressWarnings(seu <- SCTransform(seu, variable.features.n = vfn, assay = "unspliced_RNA", new.assay.name = "unspliced_SCT", verbose = FALSE))


  #Save seuart object
  saveRDS(seu, file = paste0(paste0("./",sample.name,"/"),sample.name,"_sct.rds"))
  return
}
