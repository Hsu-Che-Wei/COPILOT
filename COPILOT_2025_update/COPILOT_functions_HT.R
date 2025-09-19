ToMonocle3 <- function(seurat_object,
                       scale_all = FALSE,
                       assay = "SCT",
                       reduction_for_projection = "PCA",
                       UMAP_cluster_slot = NULL){
  
  if(scale_all){
    message("Getting residuals for all Seurat genes in chosen assay slot and placing in scale.data")
    seurat_genes <- rownames(seurat_object[[assay]])
    remaining_genes <- setdiff(seurat_genes, rownames(seurat_object[[assay]]@scale.data))
    if(assay == "SCT"){
      seurat_object <- Seurat::GetResidual(seurat_object, features = remaining_genes, assay = assay, umi.assay = "RNA")
    } else {
      seurat_object <- Seurat::ScaleData(seurat_object, features = rownames(seurat_object[[assay]]))
    }
  }
  
  #We prep the seurat object by creating gene loadings for ALL genes in the Seurat scale.data slot. This is done to allow downstream monocle3 functions on gene_modules to work appropriately.
  message("Projecting gene loadings for all Seurat genes in scale.data slot")
  seurat_object <- Seurat::ProjectDim(seurat_object, reduction = reduction_for_projection, assay = assay)
  
  ##################
  
  message("Initializing CDS object")
  
  #Extract Seurat's log-transformed values
  expression_matrix <- Seurat::GetAssayData(seurat_object, assay = assay, slot = "data")
  #Extract Seurat meta_data
  meta_data <- seurat_object@meta.data
  #Extract gene names from Seurat object SCT slot to make CDS
  seurat_genes <- data.frame(gene_short_name = rownames(seurat_object[[assay]]),
                             row.names = rownames(seurat_object[[assay]]))
  new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)

  
    
  ##################
  
  message("Making an SCE object from the Seurat object to facilitate transfer of information from SCE to CDS")
  sce <- as.SingleCellExperiment(seurat_object, assay = assay)
  message("Loading in all Seurat reductions (PCA, HARMONY, UMAP, etc.) into CDS")
  SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)
  message("Loading in specified Seurat assay into CDS")
  SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)
  message("Loading in Seurat gene names into CDS")
  SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
  SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)
  message("Loading in Seurat gene loadings into CDS")
  new_cds@preprocess_aux$gene_loadings <- seurat_object@reductions[[reduction_for_projection]]@feature.loadings.projected
  
  ##################
  
  message("Get user specified selected clusters (or active idents) from Seurat and load into CDS")
  if(is.null(UMAP_cluster_slot)){
    list_cluster <- Idents(seurat_object)
  } else {
    Idents(seurat_object) <- UMAP_cluster_slot
    list_cluster <- Idents(seurat_object)
  }
  new_cds@clusters[["UMAP"]]$clusters <- list_cluster
  #The next two commands are run in order to allow "order_cells" to be run in monocle3
  rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
  colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL
  
  ##################
  
  message("Setting all cells as belonging to one partition (multiple partitions not supported yet)")
  recreate_partition <- c(rep(1, length(new_cds@colData@rownames)))
  names(recreate_partition) <- new_cds@colData@rownames
  recreate_partition <- as.factor(recreate_partition)
  new_cds@clusters[["UMAP"]]$partitions <- recreate_partition
  
  ##################
  message("Done")
  new_cds
}

IntToMonocle3 <- function(seurat_object,
                       scale_all = FALSE,
                       assay = "integrated",
                       reduction_for_projection = "PCA",
                       UMAP_cluster_slot = NULL){
  
  if(scale_all){
    message("Getting residuals for all Seurat genes in chosen assay slot and placing in scale.data")
    seurat_genes <- rownames(seurat_object[[assay]])
    remaining_genes <- setdiff(seurat_genes, rownames(seurat_object[[assay]]@scale.data))
    if(assay == "SCT"){
      seurat_object <- Seurat::GetResidual(seurat_object, features = remaining_genes, assay = assay, umi.assay = "RNA")
    } else {
      seurat_object <- Seurat::ScaleData(seurat_object, features = rownames(seurat_object[[assay]]))
    }
  }
  
  #We prep the seurat object by creating gene loadings for ALL genes in the Seurat scale.data slot. This is done to allow downstream monocle3 functions on gene_modules to work appropriately.
  message("Projecting gene loadings for all Seurat genes in scale.data slot")
  seurat_object <- Seurat::ProjectDim(seurat_object, reduction = reduction_for_projection, assay = assay)
  
  ##################
  
  message("Initializing CDS object")
  
  #Extract Seurat's log-transformed values
  expression_matrix <- seurat_object@assays$integrated@data
  expression_matrix[which(expression_matrix < 0)]=0
  #Extract Seurat meta_data
  meta_data <- seurat_object@meta.data
  #Extract gene names from Seurat object integrated slot to make CDS
  seurat_genes <- data.frame(gene_short_name = rownames(seurat_object[[assay]]),
                             row.names = rownames(seurat_object[[assay]]))
  new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)

  
    
  ##################
  
  message("Making an SCE object from the Seurat object to facilitate transfer of information from SCE to CDS")
  sce <- as.SingleCellExperiment(seurat_object, assay = assay)
  message("Loading in all Seurat reductions (PCA, HARMONY, UMAP, etc.) into CDS")
  SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)
  message("Loading in specified Seurat assay into CDS")
  SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)
  message("Loading in Seurat gene names into CDS")
  SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
  SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)
  message("Loading in Seurat gene loadings into CDS")
  new_cds@preprocess_aux$gene_loadings <- seurat_object@reductions[[reduction_for_projection]]@feature.loadings.projected
  
  ##################
  
  message("Get user specified selected clusters (or active idents) from Seurat and load into CDS")
  if(is.null(UMAP_cluster_slot)){
    list_cluster <- Idents(seurat_object)
  } else {
    Idents(seurat_object) <- UMAP_cluster_slot
    list_cluster <- Idents(seurat_object)
  }
  new_cds@clusters[["UMAP"]]$clusters <- list_cluster
  #The next two commands are run in order to allow "order_cells" to be run in monocle3
  rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
  colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL
  
  ##################
  
  message("Setting all cells as belonging to one partition (multiple partitions not supported yet)")
  recreate_partition <- c(rep(1, length(new_cds@colData@rownames)))
  names(recreate_partition) <- new_cds@colData@rownames
  recreate_partition <- as.factor(recreate_partition)
  new_cds@clusters[["UMAP"]]$partitions <- recreate_partition
  
  ##################
  message("Done")
  new_cds
}

plain <- function(x,...) {
    prettyNum(x,big.mark = ",", scientific = FALSE)
}

bc_rank_plot <- function(lnc_gg,save){
    png(save, width = 7, height = 5, units = 'in', res=300)
    print({
     ggplot(lnc_gg, aes(x=rank, y=lnc, col = select, alpha=select)) +   
     geom_point(size=3) + 
     scale_x_log10(labels = plain, breaks = trans_breaks("log10", function(x) round(10^x, 0))) + 
     scale_y_log10(labels = plain, breaks = trans_breaks('log10', function(x) floor(10^x)))+
     scale_color_manual(values = c("#003B6D","#EBEDF3","#d96459","#f2e394"), name = NULL, guide = guide_legend(reverse = TRUE)) +
     scale_alpha_manual(values = c(0.8,0.8,0.8,0.8)) +
     labs(x = 'Barcodes', y = 'UMI Counts') +
     guides(alpha = FALSE, colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
     theme_bw() + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
            legend.background = element_rect(fill = 'transparent'),
            legend.position = "none",legend.key=element_blank(),
            plot.margin=unit(c(0,1,0,0),"cm"))
    })
    dev.off()
}

UMI_hist_plot <- function(lnc_gg,save,legend.position){
    png(save, width = 7, height = 5, units = 'in', res=300)
    print({
    ggplot(lnc_gg, aes(lnc, fill = select)) +
    geom_histogram(bins = 100)+
    scale_fill_manual(values = c("#003B6D","#EBEDF3","#d96459","#f2e394"))+
    xlab("UMI Counts")+
    ylab("Freqency")+
    scale_y_continuous(labels = plain)+
    scale_x_log10(labels = plain, breaks = trans_breaks("log10", function(x) round(10^x, 0))) + 
    theme(panel.background=element_blank(),plot.background=element_blank(),legend.title = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.position=legend.position,
            plot.margin=unit(c(0,1,0,0),"cm"))
    })
    dev.off()
}

Gene_hist_plot <- function(lng_gg,save){
    png(save, width = 7, height = 5, units = 'in', res=300)
    print({
    ggplot(lng_gg, aes(lng, fill = select)) +
    geom_histogram(bins = 100)+
    scale_fill_manual(values = c("#003B6D","#EBEDF3","#d96459","#f2e394"))+
    xlab("Number of Genes")+
    ylab("Freqency")+
    scale_y_continuous(labels = plain)+
    scale_x_log10(labels = plain, breaks = trans_breaks("log10", function(x) round(10^x, 0))) + 
    theme(panel.background=element_blank(),plot.background=element_blank(),legend.title = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),legend.position="none",
            plot.margin=unit(c(0,1,0,0),"cm"))
    })
    dev.off()
}

print_HTML_no_analysis <- function(parameters, cell_stats, seq_stats, sample_stats, dir, sample.name){
    system(paste0('base64 ', dir, '/UMI_hist_plot.png > ', dir, '/UMI_hist.txt'))
    system(paste0('base64 ', dir, '/Gene_hist_plot.png > ', dir, '/Gene_hist.txt'))
    system(paste0('base64 ', dir, '/bc_rank_plot.png > ', dir, '/barcode_rank.txt'))
    b64_uh <- readChar(paste0(dir, '/UMI_hist.txt'), file.info(paste0(dir, '/UMI_hist.txt'))$size)
    b64_gh <- readChar(paste0(dir, '/Gene_hist.txt'), file.info(paste0(dir, '/Gene_hist.txt'))$size)
    b64_bc <- readChar(paste0(dir, '/barcode_rank.txt'), file.info(paste0(dir, '/barcode_rank.txt'))$size)

    target <- HTMLInitFile(dir, filename=paste0(sample.name, '_summary'))
    HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Calibri">', file=target)
    HTML("<div id='Page1'>", file=target)
    HTML("<div class='topnav'>", file=target)
    HTML(paste0('<a class="notactive" ', 'href="#" ', 'onclick="return show(',"'Page2','Page1'",')"', ";>Analysis</a>"), file=target)
    HTML(paste0('<a class="active" ', 'href="#">Summary</a>'), file=target)
    HTML("</div>", file=target)
    HTML("<div class='title'>", file=target)
    HTML.title(paste0(sample.name,' Summary'), HR=1, file = target)
    HTML.title('Processed by COPILOT', HR=3, file = target)
    HTML("</div>", file = target)
    
    HTML("<div id='wrapper'>", file=target)
    
    HTML("<div class='boxed' id='left' align='center'>", file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Parameters', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', parameters$para, '</td> <td align="right" style="color:#72BE4B;">', parameters$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Cell Stats', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', cell_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', cell_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Sequencing Stats', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', seq_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', seq_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Sample Stats', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', sample_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', sample_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML("</div>", file = target)
    
    HTML("<div class='boxed' id='right' align='center'>", file=target)
    HTML("<div class='subtitle'>", file=target)
    HTML.title('UMI Counts Histogram', HR=3, file=target)
    HTML("</div>", file = target)
    HTML(paste0("<img src='data:image/png;base64,", b64_uh, "' width=90%>"), file=target)
    HTML("<div class='subtitle'>", file=target)
    HTML.title('Number of Genes Histogram', HR=3, file=target)
    HTML("</div>", file = target)
    HTML(paste0("<img src='data:image/png;base64,", b64_gh, "' width=90%>"), file=target)
    HTML("<div class='subtitle'>", file=target)
    HTML.title('Barcode Rank Plot', HR=3, file=target)
    HTML("</div>", file = target)
    HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
    HTML("</div>", file = target)

    HTML("</div>", file = target)
    HTML("</div>", file = target)
    
    HTML('<div id="Page2" style="display:none">', file=target)
    HTML("<div class='topnav'>", file=target)
    HTML(paste0('<a class="active" ', 'href="#">Analysis</a>'), file=target)
    HTML(paste0('<a class="notactive" ', 'href="#" ', 'onclick="return show(',"'Page1','Page2'",')"',";>Summary</a>"), file=target)
    HTML("</div>", file=target)
    HTML("<div class='title'>", file=target)
    HTML.title(paste0(sample.name,' Analysis'), HR=1, file = target)
    HTML.title('Processed by COPILOT', HR=3, file = target)
    HTML("</div>", file = target)
    
    HTML("<div id='wrapper'>", file=target)
    HTML("<div class='boxed' id='center' align='center'>", file=target)
    HTML("</div>", file = target)
    HTML("</div>", file = target)
    HTML("</div>", file = target) 
    HTML('<style type="text/css">
        .title {
            background-color: #003B6D;
            color: white;
            padding: 0px;
            padding-left: 20px;
            padding-bottom: 0px;
            position: fixed;
            top: 0;
            left: 0;
            z-index: 100;
            width: 100%;
            text-align: left;
            LINE-HEIGHT:20px;
            margin: 0px;
        }
        caption {
            text-align: left;
            font-weight: bolder;
        }
        .subtitle{
            #background-color: #ff9900;
            padding-top: 36px;
            margin: 10px;
            text-align: left;
            font-weight: bolder;
        }
        .boxed {
            #background-color: #ff9900;
            #border: 1px solid #868D96;
            padding: 10px;
            margin: 30px;
        }
        h1 {
            font-family: Calibri;
            font-size: 40px;
        }
        h2 {
            font-family: Calibri;
            font-size: 32px;
        }
        h3 {
            font-family: Calibri;
            font-size: 26px;
        }
        #wrapper {
            display: flex;
            align-items: flex-start; 
        }
        #left {
            width: 50%;
        }
        #right {
            width: 50%;
        }
        #center {
            width: 100%;
        }
        table {
            #background-color: #ff4040;
            padding: 0px;
            margin: 0px;
            font-family: "Calibri";
            font-size: 22px;
            #border: 1px solid #CDCDCD;
        }
        td {
            border-bottom: 1px solid #CDCDCD;
        }
        #mathplayer{
            height: 80px;
        }
        /* Add a black background color to the top navigation */
        .topnav {
          background-color: transparent;
          overflow: hidden;
          top: 56px;
          right: 0;
          width: 100%;
          position: fixed;
          padding: 0px;
          margin: 0px;
          z-index: 999;
          text-decoration: none;
        }

        /* Style the links inside the navigation bar */
        .topnav a {
          #border-bottom: 1px solid white;
          float: right;
          color: white;
          text-align: center;
          text-decoration: none;
          font-size: 26px;
          font-family: Calibri;
          padding: 20px;
          padding-top: 10px;
          padding-bottom: 6px;
          margin: 0px; 
        }

        /* Change the color of links on hover */
        .topnav a:hover {
          #border-bottom: 1px solid black;
          background-color: white;
          color: black;
          float: right;
          color: "#EBEDF3";
          text-align: center;
          text-decoration: none;
          font-size: 26px;
          font-family: Calibri;
          padding: 20px;
          padding-top: 10px;
          padding-bottom: 6px;
          margin: 0px; 
        }
        
        /* Add a color to the active/current link */
        .topnav a.active {
          #border-bottom: 1px solid white;
          background-color: white;
          color: black;
          float: right;
          color: "#EBEDF3";
          text-align: center;
          text-decoration: none;
          font-size: 26px;
          font-family: Calibri;
          padding: 20px;
          padding-top: 10px;
          padding-bottom: 6px;
          margin: 0px; 
        }
        </style> </head>', file=target)
    HTMLEndFile()
    system(paste0("sed -i '/<hr/d' ", dir, "/", sample.name, "_summary.html"))
    system(paste0("sed -i '/<font/d' ", dir, "/", sample.name, "_summary.html"))
    system(paste0("sed -i '/R2HTML/d' ", dir, "/", sample.name, "_summary.html"))
    system(paste0("sed -i 's/nequations=0;/function show(shown, hidden) {document.getElementById(shown).style.display=",'"block"',';document.getElementById(hidden).style.display="none"',";return false;}/g' ", dir, "/", sample.name, "_summary.html"))
}

print_HTML <- function(parameters, cell_stats, seq_stats, sample_stats, dir, sample.name){
    system(paste0('base64 ', dir, '/UMI_hist_plot.png > ', dir, '/UMI_hist.txt'))
    system(paste0('base64 ', dir, '/Gene_hist_plot.png > ', dir, '/Gene_hist.txt'))
    system(paste0('base64 ', dir, '/bc_rank_plot.png > ', dir, '/barcode_rank.txt'))
    system(paste0('base64 ', dir, '/feature_plot.png > ', dir, '/feature.txt'))
    b64_uh <- readChar(paste0(dir, '/UMI_hist.txt'), file.info(paste0(dir, '/UMI_hist.txt'))$size)
    b64_gh <- readChar(paste0(dir, '/Gene_hist.txt'), file.info(paste0(dir, '/Gene_hist.txt'))$size)
    b64_bc <- readChar(paste0(dir, '/barcode_rank.txt'), file.info(paste0(dir, '/barcode_rank.txt'))$size)
    b64_ft <- readChar(paste0(dir, '/feature.txt'), file.info(paste0(dir, '/feature.txt'))$size)

    target <- HTMLInitFile(dir, filename=paste0(sample.name, '_summary'))
    HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Calibri">', file=target)
    HTML("<div id='Page1'>", file=target)
    HTML("<div class='topnav'>", file=target)
    HTML(paste0('<a class="notactive" ', 'href="#" ', 'onclick="return show(',"'Page2','Page1'",')"', ";>Analysis</a>"), file=target)
    HTML(paste0('<a class="active" ', 'href="#">Summary</a>'), file=target)
    HTML("</div>", file=target)
    HTML("<div class='title'>", file=target)
    HTML.title(paste0(sample.name,' Summary'), HR=1, file = target)
    HTML.title('Processed by COPILOT', HR=3, file = target)
    HTML("</div>", file = target)
    
    HTML("<div id='wrapper'>", file=target)
    
    HTML("<div class='boxed' id='left' align='center'>", file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Parameters', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', parameters$para, '</td> <td align="right" style="color:#72BE4B;">', parameters$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Cell Stats', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', cell_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', cell_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Sequencing Stats', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', seq_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', seq_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML('<table style="width:100%">', file=target)
    HTML('<caption>',file=target)
    HTML.title('Sample Stats', HR=3, file=target)
    HTML('</caption>',file=target)
    HTML(paste('<tr> <td style="color:#444444;">', sample_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', sample_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML("</div>", file = target)
    
    HTML("<div class='boxed' id='right' align='center'>", file=target)
    HTML("<div class='subtitle'>", file=target)
    HTML.title('UMI Counts Histogram', HR=3, file=target)
    HTML("</div>", file = target)
    HTML(paste0("<img src='data:image/png;base64,", b64_uh, "' width=90%>"), file=target)
    HTML("<div class='subtitle'>", file=target)
    HTML.title('Number of Genes Histogram', HR=3, file=target)
    HTML("</div>", file = target)
    HTML(paste0("<img src='data:image/png;base64,", b64_gh, "' width=90%>"), file=target)
    HTML("<div class='subtitle'>", file=target)
    HTML.title('Barcode Rank Plot', HR=3, file=target)
    HTML("</div>", file = target)
    HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
    HTML("</div>", file = target)

    HTML("</div>", file = target)
    HTML("</div>", file = target)
    
    HTML('<div id="Page2" style="display:none">', file=target)
    HTML("<div class='topnav'>", file=target)
    HTML(paste0('<a class="active" ', 'href="#">Analysis</a>'), file=target)
    HTML(paste0('<a class="notactive" ', 'href="#" ', 'onclick="return show(',"'Page1','Page2'",')"',";>Summary</a>"), file=target)
    HTML("</div>", file=target)
    HTML("<div class='title'>", file=target)
    HTML.title(paste0(sample.name,' Analysis'), HR=1, file = target)
    HTML.title('Processed by COPILOT', HR=3, file = target)
    HTML("</div>", file = target)
    
    HTML("<div id='wrapper'>", file=target)
    HTML("<div class='boxed' id='center' align='center'>", file=target)
    HTML(paste0("<img src='data:image/png;base64,", b64_ft, "' width=100%>"), file=target)
    HTML("</div>", file = target)
    HTML("</div>", file = target)
    HTML("</div>", file = target) 
    HTML('<style type="text/css">
        .title {
            background-color: #003B6D;
            color: white;
            padding: 0px;
            padding-left: 20px;
            padding-bottom: 0px;
            position: fixed;
            top: 0;
            left: 0;
            z-index: 100;
            width: 100%;
            text-align: left;
            LINE-HEIGHT:20px;
            margin: 0px;
        }
        caption {
            text-align: left;
            font-weight: bolder;
        }
        .subtitle{
            #background-color: #ff9900;
            padding-top: 36px;
            margin: 10px;
            text-align: left;
            font-weight: bolder;
        }
        .boxed {
            #background-color: #ff9900;
            #border: 1px solid #868D96;
            padding: 10px;
            margin: 30px;
        }
        h1 {
            font-family: Calibri;
            font-size: 40px;
        }
        h2 {
            font-family: Calibri;
            font-size: 32px;
        }
        h3 {
            font-family: Calibri;
            font-size: 26px;
        }
        #wrapper {
            display: flex;
            align-items: flex-start; 
        }
        #left {
            width: 50%;
        }
        #right {
            width: 50%;
        }
        #center {
            width: 100%;
        }
        table {
            #background-color: #ff4040;
            padding: 0px;
            margin: 0px;
            font-family: "Calibri";
            font-size: 22px;
            #border: 1px solid #CDCDCD;
        }
        td {
            border-bottom: 1px solid #CDCDCD;
        }
        #mathplayer{
            height: 80px;
        }
        /* Add a black background color to the top navigation */
        .topnav {
          background-color: transparent;
          overflow: hidden;
          top: 56px;
          right: 0;
          width: 100%;
          position: fixed;
          padding: 0px;
          margin: 0px;
          z-index: 999;
          text-decoration: none;
        }

        /* Style the links inside the navigation bar */
        .topnav a {
          #border-bottom: 1px solid white;
          float: right;
          color: white;
          text-align: center;
          text-decoration: none;
          font-size: 26px;
          font-family: Calibri;
          padding: 20px;
          padding-top: 10px;
          padding-bottom: 6px;
          margin: 0px; 
        }

        /* Change the color of links on hover */
        .topnav a:hover {
          #border-bottom: 1px solid black;
          background-color: white;
          color: black;
          float: right;
          color: "#EBEDF3";
          text-align: center;
          text-decoration: none;
          font-size: 26px;
          font-family: Calibri;
          padding: 20px;
          padding-top: 10px;
          padding-bottom: 6px;
          margin: 0px; 
        }
        
        /* Add a color to the active/current link */
        .topnav a.active {
          #border-bottom: 1px solid white;
          background-color: white;
          color: black;
          float: right;
          color: "#EBEDF3";
          text-align: center;
          text-decoration: none;
          font-size: 26px;
          font-family: Calibri;
          padding: 20px;
          padding-top: 10px;
          padding-bottom: 6px;
          margin: 0px; 
        }
        </style> </head>', file=target)
    HTMLEndFile()
    system(paste0("sed -i '/<hr/d' ", dir, "/", sample.name, "_summary.html"))
    system(paste0("sed -i '/<font/d' ", dir, "/", sample.name, "_summary.html"))
    system(paste0("sed -i '/R2HTML/d' ", dir, "/", sample.name, "_summary.html"))
    system(paste0("sed -i 's/nequations=0;/function show(shown, hidden) {document.getElementById(shown).style.display=",'"block"',';document.getElementById(hidden).style.display="none"',";return false;}/g' ", dir, "/", sample.name, "_summary.html"))
}

cor.anno.at <- function (seu, dir_to_bulk, unwanted.genes, res, res.pseudo.bulk, mt.pattern, cp.pattern) {

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
    seu <- FindClusters(seu, resolution = res.pseudo.bulk, algorithm = 3) 
        ))
  
	# Pseudo-bulk annotation
	# Pool the expression value in each cluster
	load(file=dir_to_bulk)                         
	
	afm <- as.matrix(seu@assays$SCT@data)
	pooled <- matrix(nrow=nrow(afm), ncol = 0)
	for (i in 1:(length(unique(seu@meta.data$seurat_clusters)))) {
	m <- afm[,which(seu@meta.data$seurat_clusters==i)]
	pooled <- cbind(pooled, rowSums(m)/ncol(m))
	}
	 
	time <- Reduce(merge.rownames, list(time,pooled))
	celltype <- Reduce(merge.rownames, list(celltype,pooled)) 
	
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
	
	seucluster.celltype.ID <- celltype_ident
	seucluster.timezone.ID <- time_ident
	seucluster.celltype.cor <- celltype_max
	seucluster.timezone.cor <- time_max
	seucluster.celltype.pvalue <- celltype_maxp
    seucluster.timezone.pvalue <- time_maxp
	
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
	
	seucluster.celltype.ID.P <- celltype_ident
	seucluster.timezone.ID.P <- time_ident
	seucluster.celltype.cor.P <- celltype_max
	seucluster.timezone.cor.P <- time_max
	seucluster.celltype.pvalue.P <- celltype_maxp
    seucluster.timezone.pvalue.P <- time_maxp   
	            
	time <- Reduce(merge.rownames, list(Long,pooled))
	celltype <- Reduce(merge.rownames, list(Rad,pooled)) 
	                         
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
	
	seucluster.Rad.ID <- celltype_ident
	seucluster.Long.ID <- time_ident
	seucluster.Rad.cor <- celltype_max
	seucluster.Long.cor <- time_max
	seucluster.Rad.pvalue <- celltype_maxp
    seucluster.Long.pvalue <- time_maxp
	                         
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
	
	seucluster.Rad.ID.P <- celltype_ident
	seucluster.Long.ID.P <- time_ident
	seucluster.Rad.cor.P <- celltype_max
	seucluster.Long.cor.P <- time_max
	seucluster.Rad.pvalue.P <- celltype_maxp
    seucluster.Long.pvalue.P <- time_maxp
	  
	#store values
	rc.ori.ident <- seu@active.ident
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.ID <- as.character(seu@active.ident)
	                         
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.ID.P <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.ID <- as.character(seu@active.ident)
	                         
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.ID.P <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Rad.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Rad.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Rad.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Rad.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Rad.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Rad.ID <- as.character(seu@active.ident)
	                         
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Rad.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Rad.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Rad.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Rad.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Rad.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Rad.ID.P <- as.character(seu@active.ident)
	                         
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Long.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Long.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Long.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Long.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Long.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Long.ID <- as.character(seu@active.ident)
	                         
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Long.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Long.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Long.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Long.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.Long.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.Long.ID.P <- as.character(seu@active.ident)
	
	suppressMessages(suppressWarnings(
	seu <- FindClusters(seu, resolution = res, algorithm = 4) 
	    ))
    
    return(seu)                             
}

cor.anno.os <- function (seu, dir_to_bulk, unwanted.genes, res, res.pseudo.bulk, mt.pattern, cp.pattern) {

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
    seu <- FindClusters(seu, resolution = res.pseudo.bulk, algorithm = 3) 
        ))
        
	# Pseudo-bulk annotation
	# Pool the expression value in each cluster
	load(file=dir_to_bulk)                         
	
	afm <- as.matrix(seu@assays$SCT@data)
	pooled <- matrix(nrow=nrow(afm), ncol = 0)
	for (i in 1:(length(unique(seu@meta.data$seurat_clusters)))) {
	m <- afm[,which(seu@meta.data$seurat_clusters==i)]
	pooled <- cbind(pooled, rowSums(m)/ncol(m))
	}
	
	pooled <- pooled[use.genes,]
	 
	time <- Reduce(merge.rownames, list(time_MSU,pooled))
	celltype <- Reduce(merge.rownames, list(celltype_MSU,pooled)) 
	celltype <- celltype[,-c(3,6:12)]   
	
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
	
	seucluster.celltype.plate.ID.P <- celltype_ident
	seucluster.timezone.ID.P <- time_ident
	seucluster.celltype.plate.cor.P <- celltype_max
	seucluster.timezone.cor.P <- time_max
	seucluster.celltype.plate.pvalue.P <- celltype_maxp
	seucluster.timezone.pvalue.P <- time_maxp
	
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
	
	seucluster.celltype.plate.ID <- celltype_ident
	seucluster.timezone.ID <- time_ident
	seucluster.celltype.plate.cor <- celltype_max
	seucluster.timezone.cor <- time_max
	seucluster.celltype.plate.pvalue <- celltype_maxp
	seucluster.timezone.pvalue <- time_maxp
	
	celltype <- Reduce(merge.rownames, list(celltype_MSU,pooled))
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
	
	seucluster.celltype.field.ID.P <- celltype_ident
	seucluster.celltype.field.cor.P <- celltype_max
	seucluster.celltype.field.pvalue.P <- celltype_maxp  
	                  
    celltype_stat <- suppressWarnings(sapply(7:ncol(celltype), function(i) sapply(1:6, function(j) cor.test(celltype[,i],celltype[,j],method = "spearman")[c(3,4)])))
    celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
    celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]                        
    celltype_max <- sapply(1:(ncol(celltype)-6), function(i) max(as.numeric(celltype_cor[,i])))
    celltype_ident <- sapply(1:(ncol(celltype)-6), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
    celltype_maxp <- sapply(1:(ncol(celltype)-6), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])   
    names(celltype_max) <- celltype_ident 
	
	seucluster.celltype.field.ID <- celltype_ident
	seucluster.celltype.field.cor <- celltype_max
	seucluster.celltype.field.pvalue <- celltype_maxp                   
	
	
	#store values
	rc.ori.ident <- seu@active.ident
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.plate.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.plate.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.plate.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.plate.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.plate.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.plate.ID.P <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.field.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.field.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.field.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.field.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.field.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.field.ID.P <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.cor.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.cor.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.pvalue.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.pvalue.P <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.ID.P
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.ID.P <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.plate.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.plate.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.plate.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.plate.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.plate.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.plate.ID <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.field.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.field.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.field.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.field.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.celltype.field.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.celltype.field.ID <- as.character(seu@active.ident)
	
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.cor
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.cor <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.pvalue
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.pvalue <- as.numeric(levels(seu@active.ident)[seu@active.ident])
	seu@active.ident <- rc.ori.ident
	new.cluster.ids <- seucluster.timezone.ID
	names(new.cluster.ids) <- levels(seu)
	seu <- RenameIdents(object = seu, new.cluster.ids)
	seu@meta.data$seucluster.timezone.ID <- as.character(seu@active.ident)
	
	suppressMessages(suppressWarnings(
	seu <- FindClusters(seu, resolution = res, algorithm = 4) 
	    ))

    return(seu)                             
}

feature_plot_at <- function(seu,res,doublet.rate,dir_to_color_scheme,save){
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
    p9 <- DimPlot(seu, reduction = "umap", label = TRUE)+ggtitle(paste0("Leiden Clustering (resolution ",res,")"))+theme(legend.position="none",plot.margin = unit(c(2.4,0,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    gl <- lapply(list(p1, p2, p3, p4, p5, p6, p7, p8, p9), ggplotGrob)
    gwidth <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
    gl <- lapply(gl, "[[<-", "widths", value = gwidth)
    
    png(save, width = 22, height = 40, units = 'in', res=300)
    print(
    gridExtra::grid.arrange(grobs=gl, ncol=2)
    )
    dev.off()
}

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
    p8 <- DimPlot(seu, reduction = "umap", label = TRUE)+ggtitle(paste0("Leiden Clustering (resolution ",res,")"))+theme(legend.position="none",plot.margin = unit(c(2.4,0,0.8,0), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    gl <- lapply(list(p1, p2, p3, p4, p5, p6, p7, p8), ggplotGrob)
    gwidth <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
    gl <- lapply(gl, "[[<-", "widths", value = gwidth)
    
    png(save, width = 22, height = 32, units = 'in', res=300)
    print(
    gridExtra::grid.arrange(grobs=gl, ncol=2)
    )
    dev.off()
}

feature_plot <- function(seu,res,doublet.rate,save){
    seu@meta.data$log_nCount_RNA <- log10(seu@meta.data$nCount_RNA);
    p1 <- FeaturePlot(seu, reduction = "umap", features = c("log_nCount_RNA"))+ggtitle("log10 UMI Counts")+theme(plot.margin = unit(c(2.4,0.8,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    seu@meta.data$log_nFeature_RNA <- log10(seu@meta.data$nFeature_RNA);  
    p2 <- FeaturePlot(seu, reduction = "umap", features = c("log_nFeature_RNA"))+ggtitle("log10 Number of Genes")+theme(plot.margin = unit(c(2.4,0,0.8,0.8), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    p3 <- FeaturePlot(seu, reduction = "umap", features = c("percent.mt"))+ggtitle("Percent Mitochondrial")+theme(plot.margin = unit(c(2.4,0.8,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    p4 <- DimPlot(seu, reduction = "umap", group.by = colnames(seu@meta.data)[grep("DF.classifications",colnames(seu@meta.data))], order = c("Doublet","Singlet"),cols = c("#cccccc", "#ff4040"))+ggtitle(paste0("Doublet Rate ",doublet.rate," %"))+theme(plot.margin = unit(c(2.4,0,0.8,0.8), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    p5 <- DimPlot(seu, reduction = "umap", label = TRUE)+ggtitle(paste0("Leiden Clustering (resolution ",res,")"))+theme(legend.position="none",plot.margin = unit(c(2.4,0.8,0.8,4), "cm"),plot.title = element_text(hjust = 0,vjust=2,size= 22, face="plain"));
    gl <- lapply(list(p1, p2, p3, p4, p5), ggplotGrob)
    gwidth <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
    gl <- lapply(gl, "[[<-", "widths", value = gwidth)
    
    png(save, width = 22, height = 24, units = 'in', res=300)
    print(
    gridExtra::grid.arrange(grobs=gl, ncol=2)
    )
    dev.off()
}

find_modes <- function(x) {
  modes <- NULL
  for ( i in 2:(length(x)-1) ){
    if ( (x[i] > x[i-1]) & (x[i] > x[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  return(modes)
}

est.doub.rate <- function (ncell) {
    round(ncell*4e-04 + 1.685e-15,2)
}

copilot_prepros <- function(sample.name, spliced.mtx = NULL, unspliced.mtx = NULL, dir_to_h5 = NULL, filtered.mtx.output.dir = NULL, species.name = "Not Provided", transcriptome.name = "Not Provided", sample.stats = NULL, mt.pattern = NULL, mt.threshold = 5, remove.mt.enrich.cell = TRUE, keep.top.x.cells = NULL, cp.pattern = NULL, top.percent = 1, filtering.ratio = 1, estimate.doublet.rate = TRUE, doublet.rate = 5, remove.doublet = TRUE, filter.by.knee = FALSE,
                            do.seurat = TRUE, do.annotation = TRUE, unwanted.genes = NULL, HVG = FALSE, HVGN = 2000, dir_to_bulk, res = 0.5, res.pseudo.bulk = 50, do.dynamic.filtering = TRUE, min.UMI.low.quality = 100, min.UMI.high.quality = 300, legend.position = c(0.8,0.8)){
  
  if(is.null(spliced.mtx) && is.null(unspliced.mtx) && is.null(dir_to_h5)){
    # load raw mtx
    spliced <- readMM(paste0("./",sample.name,"/spliced.mtx"))
    rownames(spliced) <- as.character(read.table(paste0("./",sample.name,"/spliced.barcodes.txt"), header=F)$V1)
    colnames(spliced) <- as.character(read.table(paste0("./",sample.name,"/spliced.genes.txt"), header=F)$V1)
    spliced <- t(spliced)
    unspliced <- readMM(paste0("./",sample.name,"/unspliced.mtx"))
    rownames(unspliced) <- as.character(read.table(paste0("./",sample.name,"/unspliced.barcodes.txt"), header=F)$V1)
    colnames(unspliced) <- as.character(read.table(paste0("./",sample.name,"/unspliced.genes.txt"), header=F)$V1)
    unspliced <- t(unspliced)
    # Create matrix spliced, unspliced and combined 
    bcs_use <- intersect(colnames(spliced),colnames(unspliced))
    tot_genes <- Matrix::rowSums(spliced)
    genes_use <- rownames(spliced)[tot_genes > 0]
    sf <- spliced[genes_use, bcs_use]
    uf <- unspliced[genes_use, bcs_use]
    afr <- sf+uf
    tot_gene <- Matrix::colSums(afr)
    bcs_use <- colnames(afr)[tot_gene > min.UMI.low.quality]
    sf <- spliced[genes_use, bcs_use]
    uf <- unspliced[genes_use, bcs_use]
    afr <- sf+uf
  } else if (is.null(dir_to_h5)) {
    spliced <- spliced.mtx
    unspliced <- unspliced.mtx
    # Create matrix spliced, unspliced and combined 
    bcs_use <- intersect(colnames(spliced),colnames(unspliced))
    tot_genes <- Matrix::rowSums(spliced)
    genes_use <- rownames(spliced)[tot_genes > 0]
    sf <- spliced[genes_use, bcs_use]
    uf <- unspliced[genes_use, bcs_use]
    afr <- sf+uf
    tot_gene <- Matrix::colSums(afr)
    bcs_use <- colnames(afr)[tot_gene > min.UMI.low.quality]
    sf <- spliced[genes_use, bcs_use]
    uf <- unspliced[genes_use, bcs_use]
    afr <- sf+uf
  } else {
    # Open the H5 file
    h5_file_f <- H5Fopen(dir_to_h5)
    # Get the names of the datasets in the H5 file
    dataset_names <- h5ls(h5_file_f)$name
    # Read a dataset from the H5 file into a matrix
    matrix_data <- h5read(h5_file_f, dataset_names[1])
    # Construct sparse matrix
    afr <- sparseMatrix(i = matrix_data$indices + 1, p = matrix_data$indptr, x = matrix_data$data, dim = matrix_data$shape, dimnames = list(matrix_data$features$id, matrix_data$barcodes))
    # Close the H5 file
    H5Fclose(h5_file_f)
    tot_genes <- Matrix::rowSums(afr)
    genes_use <- rownames(afr)[tot_genes > 0]
    tot_gene <- Matrix::colSums(afr)
    bcs_use <- colnames(afr)[tot_gene > min.UMI.low.quality]
    afr <- afr[which(rownames(afr) %in% genes_use), which(colnames(afr) %in% bcs_use)]
  }
  
  
  nc <- colSums(afr)
  ng <- colSums(afr>0)
  lnc <- log10(nc)
  pmt <- (colSums(afr[grep(paste(mt.pattern, collapse = "|"),rownames(afr)),])/nc)*100 # percent mt
  lnc.mt <- lnc[pmt >= mt.threshold]
  if (length(lnc.mt)==0){
    mmt <- log10(min.UMI.high.quality)# If no high mt cell, then set the first filtering threshold at min UMI counts threshold         
  } else if (length(lnc.mt)==1) {
    mmt <- lnc.mt # log10 mean UMI counts for cells above threshold
  } else {
    dd <- density(lnc.mt)
    mmt <- max(dd$x[find_modes(dd$y)][dd$x[find_modes(dd$y)] < max(boxplot.stats(lnc.mt)$stats)]) # find the maximum mode
  }
  if (do.dynamic.filtering){
    lcu <- as.numeric(log10(cumsum(rev(sort(nc))))) # log culmulative UMI counts
    snc <- as.numeric(rev(sort(nc))) # sorted UMI counts
    lin <- log10(seq(1,ncol(afr),1))
    dYu <- diff(lcu)/diff(lin)  # the derivative  
    dX <- rowMeans(embed(lin,2)) # centers the X values for plotting
    local_min <- min(dYu[1:(length(which(lnc>mmt))-1)]) #local min derivative
    cnu <- which(dYu<=local_min)[1]
    if (snc[cnu] < min.UMI.high.quality){
      cnu <- length(which(nc > min.UMI.high.quality))# If no local minimum, then set the first filtering threshold at min UMI counts threshold  
    }
    
    ct <- sort(nc,decreasing = TRUE)[cnu] # UMI threshold at elbow point 
    
    # Filtering on UMI counts (first round)
    filter.UMI.thres.n=c(min.UMI.low.quality,ct)
    filter.UMI.thres=c(ct,1000000)
    nf <- afr[,(colSums(afr) > filter.UMI.thres.n[1])&(colSums(afr) <= filter.UMI.thres.n[2])] # lower quality cells
    af <- afr[,(colSums(afr) > filter.UMI.thres[1])&(colSums(afr) < filter.UMI.thres[2])] # higher quality cells
    nnf <- log1p((nf/colSums(nf))*10000) #log normalization
    naf <- log1p((af/colSums(af))*10000) #log normalization
    np <- rowSums(nnf)/ncol(nnf) # low quality cell profile
    sp <- rowSums(naf)/ncol(naf) # high quality cell profile
    
    cornp <- apply(naf,2,cor,y=np)
    corsp <- apply(naf,2,cor,y=sp)
    
    sidx <- c()
    nidx <- c()
    for (i in 1:ncol(af)){
      if (corsp[i] > cornp[i]){
        sidx <- c(sidx, i)
      } else {
        nidx <- c(nidx, i)
      }
    }
    
    #Decide filtering threshold: lenient: 1, strict:0, moderate: 0 < filtering.ratio < 1
    nidx_thres <- length(nidx)*filtering.ratio #moderate
    print(paste0("threshold cell number: ",nidx_thres))
    n_iteration <- 0
    print(paste0("removed cells: ",length(nidx)))
    while(length(nidx)!=0 & length(nidx) >= nidx_thres){
      n_iteration <- n_iteration+1 
      #nf <- cbind(nf,afr[,colnames(af)[nidx]])
      nf <- cbind(nf,afr[,which(colnames(afr) %in% colnames(af)[nidx])])
      #af <- afr[,colnames(af)[sidx]]
      af <- afr[,which(colnames(afr) %in% colnames(af)[sidx])]
      nnf <- log1p((nf/colSums(nf))*10000) #log normalization
      naf <- log1p((af/colSums(af))*10000) #log normalization
      np <- rowSums(nnf)/ncol(nnf) # low quality cell profile
      sp <- rowSums(naf)/ncol(naf) # high quality cell profile
      cornp <- apply(naf,2,cor,y=np)
      corsp <- apply(naf,2,cor,y=sp)   
      sidx <- c()
      nidx <- c()
      for (i in 1:ncol(af)){            
        if (corsp[i] > cornp[i]){                
          sidx <- c(sidx, i)
        } else {                
          nidx <- c(nidx, i)
        }
      }
      print(paste0("iteration: ",n_iteration))
      print(paste0("removed cells: ",length(nidx)))
    }
    if (is.null(dir_to_h5)){
      sf <- sf[,colnames(af)]
      uf <- uf[,colnames(af)] 
    }
    
    message("Iteration finished")
  } else {
    # Filtering on UMI counts (first round)
    if (is.null(keep.top.x.cells) && filter.by.knee==FALSE){
      ct <- sort(nc,decreasing = TRUE)[length(which(nc > min.UMI.high.quality))] # UMI threshold at elbow point
    } else if (filter.by.knee==FALSE) {
      ct <- sort(nc,decreasing = TRUE)[keep.top.x.cells] # UMI threshold at elbow point
    } else {
      br.out <- barcodeRanks(afr)
      ct <- sort(nc,decreasing = TRUE)[length(which(nc > metadata(br.out)$knee))]
    }
    filter.UMI.thres.n=c(min.UMI.low.quality,ct)
    filter.UMI.thres=c(ct,1000000)
    nf <- afr[,(colSums(afr) > filter.UMI.thres.n[1])&(colSums(afr) <= filter.UMI.thres.n[2])] # lower quality cells
    af <- afr[,(colSums(afr) > filter.UMI.thres[1])&(colSums(afr) < filter.UMI.thres[2])] # higher quality cells
    if (is.null(dir_to_h5)){
      sf <- sf[,colnames(af)]
      uf <- uf[,colnames(af)]
    }
    n_iteration <- 0
  }
  
  #prepare data.frame for ggplot2
  nc <- colSums(af)
  ng <- colSums(af>0)
  lnc <- log10(nc)
  lng <-log10(ng)
  pmts <- (colSums(af[grep(paste(mt.pattern, collapse = "|"),rownames(af)),])/nc)*100 # percent mt for high quality cell
  ncn <- colSums(nf)
  ngn <- colSums(nf>0)
  lncn <- log10(ncn)
  lngn <- log10(ngn)
  pmtn <- (colSums(nf[grep(paste(mt.pattern, collapse = "|"),rownames(nf)),])/ncn)*100 # percent mt for low quality cell
  
  if (remove.mt.enrich.cell){
    select <- rep(paste("high quality"),length(lnc))
    label_mt <- paste0("percent mt >= ",mt.threshold)
    select[which(pmts >= mt.threshold)] <- label_mt # index of percent mt > mt.threshold
    quantile_idx <- which(lnc < quantile(lnc[which(select=="high quality")],probs = (1-top.percent/100))) # quantile top 1 % tail to determine outliers
    ssidx <- intersect(quantile_idx, which(select=="high quality")) # index of final high quality cell
    select[ssidx] <- paste("high quality (",length(ssidx),")",sep="")
    select[which(select=="high quality")] <- paste("top ",top.percent,"% (",length(which(select=="high quality")),")",sep="")
    nonselect <- rep("low quality",length(lncn))
    nonselect[which(pmtn >= mt.threshold)] <- label_mt # index of percent mt > mt.threshold
    select <- c(select,nonselect)
    select[which(select==label_mt)] <- paste("percent mt >= ",mt.threshold," (",length(which(select==label_mt)),")",sep="")
    select[which(select=="low quality")] <- paste0("low quality (",length(which(select=="low quality")),")")
    lncc <- c(lnc,lncn)
    lngc <- c(lng,lngn)
    lnc_gg <- data.frame(select=select,lnc=lncc)
    lnc_gg$rank <- match(lnc_gg$lnc,sort(lnc_gg$lnc,decreasing = TRUE))
    lnc_gg$lnc <- 10^lnc_gg$lnc
    lnc_gg$select <- factor(lnc_gg$select, levels = c(paste0("high quality (",length(grep("high",select)),")"), 
                                                      paste0("low quality (",length(grep("low",select)),")"),
                                                      paste0("percent mt >= ",mt.threshold," (",length(grep("percent",select)),")"),
                                                      paste0("top ",top.percent,"% (",length(grep("top",select)),")")))
    lng_gg <- data.frame(select=select,lng=lngc)
    lng_gg$rank <- match(lng_gg$lng,sort(lng_gg$lng,decreasing = TRUE))
    lng_gg$lng <- 10^lng_gg$lng
    lng_gg$select <- factor(lng_gg$select, levels = c(paste0("high quality (",length(grep("high",select)),")"), 
                                                      paste0("low quality (",length(grep("low",select)),")"),
                                                      paste0("percent mt >= ",mt.threshold," (",length(grep("percent",select)),")"),
                                                      paste0("top ",top.percent,"% (",length(grep("top",select)),")")))
    
    af <- af[,ssidx]
    if (is.null(dir_to_h5)){
      sf <- sf[,ssidx]
      uf <- uf[,ssidx]
    }
  } else {
    select <- rep(paste("high quality"),length(lnc))
    #label_mt <- paste0("percent mt >= ",mt.threshold)
    #select[which(pmts >= mt.threshold)] <- label_mt # index of percent mt > mt.threshold
    quantile_idx <- which(lnc < quantile(lnc[which(select=="high quality")],probs = (1-top.percent/100))) # quantile top 1 % tail to determine outliers
    ssidx <- intersect(quantile_idx, which(select=="high quality")) # index of final high quality cell
    select[ssidx] <- paste("high quality (",length(ssidx),")",sep="")
    select[which(select=="high quality")] <- paste("top ",top.percent,"% (",length(which(select=="high quality")),")",sep="")
    nonselect <- rep("low quality",length(lncn))
    #nonselect[which(pmtn >= mt.threshold)] <- label_mt # index of percent mt > mt.threshold
    select <- c(select,nonselect)
    #select[which(select==label_mt)] <- paste("percent mt >= ",mt.threshold," (",length(which(select==label_mt)),")",sep="")
    select[which(select=="low quality")] <- paste0("low quality (",length(which(select=="low quality")),")")
    lncc <- c(lnc,lncn)
    lngc <- c(lng,lngn)
    lnc_gg <- data.frame(select=select,lnc=lncc)
    lnc_gg$rank <- match(lnc_gg$lnc,sort(lnc_gg$lnc,decreasing = TRUE))
    lnc_gg$lnc <- 10^lnc_gg$lnc
    lnc_gg$select <- factor(lnc_gg$select, levels = c(paste0("high quality (",length(grep("high",select)),")"),
                                                      paste0("low quality (",length(grep("low",select)),")"),
                                                      #paste0("percent mt >= ",mt.threshold," (",length(grep("percent",select)),")"),
                                                      paste0("top ",top.percent,"% (",length(grep("top",select)),")")))
    lng_gg <- data.frame(select=select,lng=lngc)
    lng_gg$rank <- match(lng_gg$lng,sort(lng_gg$lng,decreasing = TRUE))
    lng_gg$lng <- 10^lng_gg$lng
    lng_gg$select <- factor(lng_gg$select, levels = c(paste0("high quality (",length(grep("high",select)),")"),
                                                      paste0("low quality (",length(grep("low",select)),")"),
                                                      #paste0("percent mt >= ",mt.threshold," (",length(grep("percent",select)),")"),
                                                      paste0("top ",top.percent,"% (",length(grep("top",select)),")")))
    
    af <- af[,ssidx]
    if (is.null(dir_to_h5)){
      sf <- sf[,ssidx]
      uf <- uf[,ssidx]
    }
  }
  
  if(is.null(dir_to_h5)){
    #kb stats
    kb_stats <- c(fromJSON(file = paste0("./",sample.name,"/inspect.json")), fromJSON(file = paste0("./",sample.name,"/run_info.json"))) # load run info 
    tech <- strsplit(kb_stats$call, '\\s')[[1]][8]
    seq_stats <- data.frame(stat = c('Number of Reads Processed', 'Reads Pseudoaligned', 'Reads on Whitelist', 'Total UMI Counts','Sequencing Technology', 'Species', 'Transcriptome'), 
                            value = prettyNum(c(kb_stats$n_processed, paste0(kb_stats$p_pseudoaligned, ' %'), paste0(round(kb_stats$percentageReadsOnWhitelist,2),' %'), sum(afr), tech, species.name, transcriptome.name), big.mark = ','))
  } else {
    seq_stats <- data.frame(stat = c('Number of Reads Processed', 'Reads Pseudoaligned', 'Reads on Whitelist', 'Total UMI Counts','Sequencing Technology', 'Species', 'Transcriptome'), 
                            value = prettyNum(c("Not Provided", "Not Provided", "Not Provided", sum(afr), "Not Provided", species.name, transcriptome.name), big.mark = ','))
  }
  
  
  #cell stats
  if (estimate.doublet.rate){
    p_cnts_in_cells <- round((sum(af)/sum(afr))*100, 2) # calculate cell stats and save to df
    med_cnts_cell <- median(colSums(af))
    med_genes_cell <- median(apply(af, 2, function(x) sum(x >= 1)))
    tot_genes_detected <- sum(rowSums(af)>=1)
    if (length(table(lnc_gg$select)[grep("mt",names(table(lnc_gg$select)))]) == 0){
      p_cnts_dying <- 0
    } else {
      p_cnts_dying <- as.numeric(round(table(lnc_gg$select)[grep("mt",names(table(lnc_gg$select)))]/nrow(lnc_gg)*100,2)) 
    }
    p_cnts_high <-  as.numeric(round(table(lnc_gg$select)[grep("high",names(table(lnc_gg$select)))]/nrow(lnc_gg)*100,2))
    doublet.rate <- est.doub.rate(ncol(af))
    cell_stats <- data.frame(stat = c('Estimated Number of High Quality Cell', 'High Quality Cell', 'Total UMI Counts in High Quality Cell','UMI Counts in High Quality Cell',
                                      'Median UMI Counts per High Quality Cell', 'Median Genes per High Quality Cell', 'Total Genes Detected in High Quality Cell'
                                      ,'Cell above Mitochondrial Expression Threshold','Estimated Doublet Rate in High Quality Cell'), 
                             value = prettyNum(c(ncol(af), paste0(p_cnts_high,' %'), sum(af), paste0(p_cnts_in_cells,' %'), med_cnts_cell,
                                                 med_genes_cell, tot_genes_detected, paste0(p_cnts_dying,' %'),paste0(doublet.rate,' %')), big.mark = ','))
    
  } else {
    p_cnts_in_cells <- round((sum(af)/sum(afr))*100, 2) # calculate cell stats and save to df
    med_cnts_cell <- median(colSums(af))
    med_genes_cell <- median(apply(af, 2, function(x) sum(x >= 1)))
    tot_genes_detected <- sum(rowSums(af)>=1)
    if (length(table(lnc_gg$select)[grep("mt",names(table(lnc_gg$select)))]) == 0){
      p_cnts_dying <- 0
    } else {
      p_cnts_dying <- as.numeric(round(table(lnc_gg$select)[grep("mt",names(table(lnc_gg$select)))]/nrow(lnc_gg)*100,2)) 
    }
    p_cnts_high <-  as.numeric(round(table(lnc_gg$select)[grep("high",names(table(lnc_gg$select)))]/nrow(lnc_gg)*100,2))
    cell_stats <- data.frame(stat = c('Estimated Number of High Quality Cell', 'High Quality Cell', 'Total UMI Counts in High Quality Cell','UMI Counts in High Quality Cell',
                                      'Median UMI Counts per High Quality Cell', 'Median Genes per High Quality Cell', 'Total Genes Detected in High Quality Cell'
                                      ,'Cell above Mitochondrial Expression Threshold'), 
                             value = prettyNum(c(ncol(af), paste0(p_cnts_high,' %'), sum(af), paste0(p_cnts_in_cells,' %'), med_cnts_cell,
                                                 med_genes_cell, tot_genes_detected, paste0(p_cnts_dying,' %')), big.mark = ','))
  }
  
  #parameters
  if (remove.doublet){
    if (estimate.doublet.rate){
      parameters <- data.frame(para = c('Iteration of Filtering', 'Mitochondrial Expression Threshold', 'Top High Quality Cell Filtered', 'Doublet Removed'), 
                               value = c(n_iteration,paste0(mt.threshold,' %'),paste0(top.percent,' %'), 'Yes'))
    } else {
      parameters <- data.frame(para = c('Iteration of Filtering', 'Mitochondrial Expression Threshold', 'Top High Quality Cell Filtered', 'Expected Doublet Rate in High Quality Cell', 'Doublet Removed'), 
                               value = c(n_iteration,paste0(mt.threshold,' %'),paste0(top.percent,' %'),paste0(doublet.rate,' %'), 'Yes'))
    }
  } else {
    if (estimate.doublet.rate){
      parameters <- data.frame(para = c('Iteration of Filtering', 'Mitochondrial Expression Threshold', 'Top High Quality Cell Filtered', 'Doublet Removed'), 
                               value = c(n_iteration,paste0(mt.threshold,' %'),paste0(top.percent,' %'), 'No'))
    } else {
      parameters <- data.frame(para = c('Iteration of Filtering', 'Mitochondrial Expression Threshold', 'Top High Quality Cell Filtered', 'Expected Doublet Rate in High Quality Cell', 'Doublet Removed'), 
                               value = c(n_iteration,paste0(mt.threshold,' %'),paste0(top.percent,' %'),paste0(doublet.rate,' %'), 'No'))
    }
  }
  
  
  
  #Plots                               
  bc_rank_plot(lnc_gg=lnc_gg, save = paste0("./",sample.name,"/bc_rank_plot.png"))
  UMI_hist_plot(lnc_gg=lnc_gg, save = paste0("./",sample.name,"/UMI_hist_plot.png"), legend.position=legend.position)
  Gene_hist_plot(lng_gg=lng_gg, save = paste0("./",sample.name,"/Gene_hist_plot.png"))
  
  if(med_cnts_cell < 1000) {
    do.annotation <- FALSE
  }
  
  #do.analysis or not
  if (!do.seurat){
    #Output html 
    print_HTML_no_analysis(parameters = parameters, cell_stats = cell_stats, seq_stats = seq_stats, sample_stats = sample.stats, dir = paste0("./",sample.name), sample.name = sample.name)
    #Save filtered raw counts
    if (is.null(filtered.mtx.output.dir)){
      if (is.null(dir_to_h5)){
        write10xCounts(paste0("./",sample.name,'/spliced_counts_filtered'), sf, overwrite=T)
        write10xCounts(paste0("./",sample.name,'/unspliced_counts_filtered'), uf, overwrite=T)
      } else {
        write10xCounts(paste0("./",sample.name,'/counts_filtered'), af, overwrite=T)
      }
    } else {
      system(paste0("mkdir ",filtered.mtx.output.dir))
      if (is.null(dir_to_h5)){
        write10xCounts(paste0(filtered.mtx.output.dir,'/spliced_counts_filtered'), sf, overwrite=T)
        write10xCounts(paste0(filtered.mtx.output.dir,'/unspliced_counts_filtered'), uf, overwrite=T)
      } else {
        write10xCounts(paste0(filtered.mtx.output.dir,'/counts_filtered'), af, overwrite=T)
      }
    }
    
  } else {
    if (is.null(dir_to_h5)){
      #Create Seurat Objects
      seu <- suppressWarnings(CreateSeuratObject(counts = af, assay = "RNA", project = sample.name))
      seu[["spliced_RNA"]] <- CreateAssayObject(sf)
      seu[["unspliced_RNA"]] <- CreateAssayObject(uf)
      DefaultAssay(seu) <- "RNA" 
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = paste(mt.pattern, collapse = "|"))
      seu[["percent.cp"]] <- PercentageFeatureSet(seu, pattern = paste(cp.pattern, collapse = "|"))
    } else {
      #Create Seurat Objects
      seu <- suppressWarnings(CreateSeuratObject(counts = af, assay = "RNA", project = sample.name))
      DefaultAssay(seu) <- "RNA" 
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = paste(mt.pattern, collapse = "|"))
      seu[["percent.cp"]] <- PercentageFeatureSet(seu, pattern = paste(cp.pattern, collapse = "|"))
    }
    
    # Whether to apply highly variable gene selection
    if (HVG) {
      vfn <- HVGN
    } else {
      vfn <- nrow(seu) 
    }
    
    # Normalization using SCTransform
    suppressWarnings(
      seu <- SCTransform(seu, variable.features.n = vfn, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
    )
    #suppressWarnings(
    #seu <- SCTransform(seu, variable.features.n = vfn, assay = "spliced_RNA", new.assay.name = "spliced_SCT", verbose = FALSE)
    #)
    #suppressWarnings(
    #seu <- SCTransform(seu, variable.features.n = vfn, assay = "unspliced_RNA", new.assay.name = "unspliced_SCT", verbose = FALSE)
    #)
    
    # Doublet estimation
    DefaultAssay(seu) <- "SCT"
    
    if(!is.null(mt.pattern) && !is.null(cp.pattern) && !is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)),grep(paste(cp.pattern, collapse = "|"),rownames(seu)),sort(match(unwanted.genes, rownames(seu))))]        
    } else if (!is.null(mt.pattern) && !is.null(cp.pattern) && is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)),grep(paste(cp.pattern, collapse = "|"),rownames(seu)))]
    } else if (!is.null(mt.pattern) && is.null(cp.pattern) && !is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)),sort(match(unwanted.genes, rownames(seu))))]
    } else if (!is.null(mt.pattern) && is.null(cp.pattern) && is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)))]
    } else if (is.null(mt.pattern) && is.null(cp.pattern) && !is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(sort(match(unwanted.genes, rownames(seu))))]
    } else if (is.null(mt.pattern) && !is.null(cp.pattern) && is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(grep(paste(cp.pattern, collapse = "|"),rownames(seu)))]
    } else if (is.null(mt.pattern) && !is.null(cp.pattern) && !is.null(unwanted.genes)){
      use.genes <- rownames(seu)[-c(grep(paste(mt.pattern, collapse = "|"),rownames(seu)),grep(paste(cp.pattern, collapse = "|"),rownames(seu)))]
    } else {
      use.genes <- rownames(seu)
    }
    
    seu <- RunPCA(seu, verbose = FALSE, approx = FALSE, npcs = 10, features=use.genes)
    nExp_poi <- round((doublet.rate/100)*ncol(seu))
    
    suppressMessages(suppressWarnings(
      seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.15, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    ))
    
    # Doublet removal
    if (remove.doublet){
      seu <- subset(seu, cells = colnames(seu)[which(seu@meta.data[,grep("DF.classifications",colnames(seu@meta.data))]=="Singlet")])
      # Normalization using SCTransform
      suppressWarnings(
        seu <- SCTransform(seu, variable.features.n = vfn, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
      )
      #suppressWarnings(
      #  seu <- SCTransform(seu, variable.features.n = vfn, assay = "spliced_RNA", new.assay.name = "spliced_SCT", verbose = FALSE)
      #)
      #suppressWarnings(
      #  seu <- SCTransform(seu, variable.features.n = vfn, assay = "unspliced_RNA", new.assay.name = "unspliced_SCT", verbose = FALSE)
      #)
    }
    DefaultAssay(seu) <- "SCT"
    
    #store misc
    seu@misc$percent_mt_raw <- list(pmt)
    seu@misc$log_nCount_raw <- list(lnc_gg)
    seu@misc$log_nFeature_raw <- list(lng_gg)
    seu@misc$parameters <- list(parameters)
    seu@misc$seq_stats <- list(seq_stats)
    seu@misc$cell_stats <- list(cell_stats)
    if (!is.null(sample.stats)){
      seu@misc$sample_stats <- list(sample.stats)
    }
    
    if (do.annotation){            
      #Annotation whole transcriptome correlation
      if (species.name=="Arabidopsis thaliana"){
        seu <- cor.anno.at(seu = seu, dir_to_bulk = dir_to_bulk, unwanted.genes = unwanted.genes, res = res, res.pseudo.bulk = res.pseudo.bulk, mt.pattern = mt.pattern, cp.pattern = cp.pattern)
        #Feature Plot
        feature_plot_at(seu = seu, res = res, doublet.rate = doublet.rate, dir_to_color_scheme = "./color_scheme_at.RData",save = paste0("./",sample.name,"/feature_plot.png"))
        #Output HTML
        print_HTML(parameters = parameters, cell_stats = cell_stats, seq_stats = seq_stats, sample_stats = sample.stats, dir = paste0("./",sample.name), sample.name = sample.name)
      } else if (species.name=="Oryza sativa"){
        seu <- cor.anno.os(seu = seu, dir_to_bulk = dir_to_bulk, unwanted.genes = unwanted.genes, res = res, res.pseudo.bulk = res.pseudo.bulk, mt.pattern = mt.pattern, cp.pattern = cp.pattern)
        #Feature Plot
        feature_plot_os(seu = seu, res = res, doublet.rate = doublet.rate, dir_to_color_scheme = "./color_scheme_os.RData",save = paste0("./",sample.name,"/feature_plot.png"))
        #Output HTML
        print_HTML(parameters = parameters, cell_stats = cell_stats, seq_stats = seq_stats, sample_stats = sample.stats, dir = paste0("./",sample.name), sample.name = sample.name)
      }
    } else {
      # Run PCA, UMAP, Clustering                         
      seu <- RunPCA(seu, verbose = FALSE, approx = FALSE, npcs = 50, features=use.genes)
      suppressMessages(suppressWarnings(                         
        seu <- RunUMAP(seu, reduction = "pca", dims = 1:50, umap.method = "umap-learn", metric = "correlation")
      ))
      suppressMessages(suppressWarnings(
        seu <- FindNeighbors(seu, reduction = "pca",dims = 1:50)
      ))
      suppressMessages(suppressWarnings(
        seu <- FindClusters(seu, resolution = res, algorithm = 4) 
      ))
      #Feature Plot
      feature_plot(seu = seu, res = res, doublet.rate = doublet.rate, save = paste0("./",sample.name,"/feature_plot.png"))
      #Output HTML
      print_HTML(parameters = parameters, cell_stats = cell_stats, seq_stats = seq_stats, sample_stats = sample.stats, dir = paste0("./",sample.name), sample.name = sample.name)
    }
    
    #Save filtered raw counts
    if (is.null(filtered.mtx.output.dir)){
      if (is.null(dir_to_h5)){
        write10xCounts(paste0("./",sample.name,'/spliced_counts_filtered'), sf, overwrite=T)
        write10xCounts(paste0("./",sample.name,'/unspliced_counts_filtered'), uf, overwrite=T)
      } else {
        write10xCounts(paste0("./",sample.name,'/counts_filtered'), af, overwrite=T)
      }
    } else {
      system(paste0("mkdir ",filtered.mtx.output.dir))
      if (is.null(dir_to_h5)){
        write10xCounts(paste0(filtered.mtx.output.dir,'/spliced_counts_filtered'), sf, overwrite=T)
        write10xCounts(paste0(filtered.mtx.output.dir,'/unspliced_counts_filtered'), uf, overwrite=T)
      } else {
        write10xCounts(paste0(filtered.mtx.output.dir,'/counts_filtered'), af, overwrite=T)
      }
    }
    #Save seuart object
    
    saveRDS(seu, file = paste0(paste0("./",sample.name,"/"),sample.name,"_COPILOT.rds"))
  }   
} 
