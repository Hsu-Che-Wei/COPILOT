download_GEO_data <- function(sample_metadata) {
  GEO_accessions = sample_metadata %>%
    filter(Database == "GEO") %>%
    pull(Accession) %>%
    unique()
  
  GEO_download = purrr::map(
    .x = GEO_accessions,
    .f = GEOquery::getGEOSuppFiles,
    baseDir = "supp_data/published_datasets",
    fetch_files = T,
    makeDirectory = F,
    filter_regex = "*.CEL")
  
  GEO_files = GEO_download %>%
    purrr::reduce(rbind) %>%
    rownames() %>%
    stringr::str_match(pattern = 'supp_published_datasets.(.*)') %>% .[,2] %>%
    tibble::enframe(name = NULL, value = "Filename") %>%
    dplyr::filter(Filename %in% sample_metadata$Filename) %>%
    dplyr::pull(Filename)
  
  return(GEO_files)
}

download_other_affy <- function(files, url) {
  download.file(url, method = 'curl', destfile = "supp_data/published_datasets/pubaffy.zip")
  extracted <- unzip("supp_data/published_datasets/pubaffy.zip", exdir = "supp_data/published_datasets", files = files)
  extracted_base <- stringr::str_match(extracted, pattern = "supp_data/published_datasets/(.*)")[,2]
  
  return(extracted_base)
}

download_AREX_data <- function(sample_metadata, url) {
  AREX_files = sample_metadata %>%
    filter(Database == "AREX") %>%
    pull(Filename) %>%
    download_other_affy(url = url)
}

read_and_normalize_affy <- function(sample_metadata_sub = NULL, probe_data) {
  if(is.null(sample_metadata_sub)) return(NULL)
  filenames = sample_metadata_sub %>%
    pull(Filename) %>%
    paste0("supp_data/published_datasets/", .)
  
  abobj_raw = affy::read.affybatch(filenames)
  
  abobj_norm = gcrma::gcrma(
    object = abobj_raw,
    affinity.source = "local",
    normalize = T)
 
  array_exprs = affy::exprs(abobj_norm) %>% 
    tibble::as_tibble(rownames = "array_element_name") %>% 
    dplyr::left_join(probe_data) %>%
    dplyr::filter(!stringr::str_detect(array_element_name, "_s_|_x_")) %>% 
    dplyr::filter(!stringr::str_detect(Alt_Locus, "///")) %>% 
    dplyr::group_by(Alt_Locus) %>% 
    dplyr::mutate(n_probes = n()) %>% 
    dplyr::filter(n_probes == 1) %>% 
    dplyr::select(-Locus, -array_element_name, -n_probes) %>% 
    tibble::column_to_rownames("Alt_Locus")

  return(array_exprs)
}

read_and_normalize_rnaseq <- function(sample_metadata_sub = NULL, rnaseq_norm_method = "DESeq2") {
  if(is.null(sample_metadata_sub)) return(NULL)
  
  rnaseq_data = readr::read_csv(paste0("supp_data/published_datasets/RNASeq_counts.csv")) %>%
    dplyr::select(Locus, one_of(sample_metadata_sub$Filename)) %>%
    tibble::column_to_rownames("Locus")
  
  colData <- tibble::tibble(Filename = colnames(rnaseq_data)) %>%
    left_join(sample_metadata_sub) %>%
    mutate(Filename = factor(Filename, levels = colnames(rnaseq_data))) %>%
    arrange(Filename)
  
  if(rnaseq_norm_method == "DESeq2") {
    dds <- DESeq2::DESeqDataSetFromMatrix(rnaseq_data, colData = colData, design = ~ Marker) %>%
      DESeq2::estimateSizeFactors() %>%
      DESeq2::estimateDispersions() %>%
      DESeq2::vst()

    rnaseq_exprs <- SummarizedExperiment::assay(dds)
  }

  if(rnaseq_norm_method == "edgeR") {
    erobj_raw <- edgeR::DGEList(rnaseq_data)
    erobj_norm = erobj_raw %>% 
      edgeR::estimateDisp() %>% 
      edgeR::calcNormFactors(method = "upperquartile")

    rnaseq_exprs = edgeR::cpm(erobj_norm)
  }
  
  return(rnaseq_exprs)
}

merge_expression_data <- function(rnaseq_exprs, array_exprs, universe) {
  if(is.null(rnaseq_exprs)) return(array_exprs[universe,])
  if(is.null(array_exprs)) return(rnaseq_exprs[universe,])
  
  rnaseq_exprs_fsqn = FSQN::quantileNormalizeByFeature(
    matrix_to_normalize = rnaseq_exprs[universe,] %>% t(),
    target_distribution_matrix = array_exprs[universe,] %>% t()) %>%
    t()
  
  expression_data_mat = cbind(
    rnaseq_exprs_fsqn[universe,],
    array_exprs[universe,])
  
  return(expression_data_mat)
  
}

build_expression_table <- function(sample_metadata, exclude_loci, probe_data = NULL, rnaseq_norm_method = "DESeq2") {
  affy_metadata <- sample_metadata %>% dplyr::filter(Type == "ATH1")
  rnaseq_metadata <- sample_metadata %>% dplyr::filter(Type == "RNASeq")
  array_exprs = NULL
  rnaseq_exprs = NULL
  
  if(dim(affy_metadata)[1] > 0) {
    array_exprs = read_and_normalize_affy(sample_metadata_sub = affy_metadata, probe_data = probe_data)
    if(rnaseq_norm_method == "edgeR") array_exprs = array_exprs^2
  }
  
  if(dim(rnaseq_metadata)[1] > 0) {
    rnaseq_exprs = read_and_normalize_rnaseq(sample_metadata_sub = rnaseq_metadata, rnaseq_norm_method = rnaseq_norm_method)
  }
  if(is.null(array_exprs) | is.null(rnaseq_exprs)) {
    universe <- union(row.names(array_exprs), row.names(rnaseq_exprs))
  } else {
    universe <- intersect(row.names(array_exprs), row.names(rnaseq_exprs))
  }
  
  universe <- tibble::tibble(Locus = universe) %>%
    dplyr::filter(stringr::str_detect(Locus, "A[Tt][1-5][Gg].....")) %>%
    dplyr::filter(!(Locus %in% exclude_loci)) %>%
    pull(Locus)
  
  expression_data_mat = merge_expression_data(array_exprs, rnaseq_exprs, universe)
  expression_data = expression_data_mat %>%
    tibble::as_tibble(rownames = "Locus") %>%
    tidyr::gather("Filename", "Expression", -Locus) %>%
    dplyr::left_join(sample_metadata) %>%
    dplyr::mutate(Sample_Name = Filename) #%>%
#    dplyr::group_by(Locus) %>%
#    dplyr::mutate(var = var(Expression)) %>%
#    dplyr::filter(var > 0.01) %>%
#    dplyr::select(-var)

  return(expression_data)
}

compute_ici_wrapper <- function(sobj_filename, spec_table, slot = "scale.data",
                                       assay = "SCT", information_level = 20, 
                                       write_file = NULL) {
  ## Read in dataset
  sobj <- readr::read_rds(sobj_filename)
  
  expression_data <- Seurat::GetAssayData(sobj, assay = assay, slot = slot) %>% 
    tibble::as_tibble(rownames = "Locus")
  future::plan("sequential")
  ici <- ICITools::compute_ici_scores(expression_data = expression_data,
                                      spec_table = spec_table,
                                      sig = T,
                                      information_level = information_level,
                                      min_spec_score = 0.15)
  
  ici <- dplyr::mutate(ici, analysis = unique(spec_table$analysis)) %>%
    dplyr::mutate(subset = unique(spec_table$subset)) %>% 
    dplyr::mutate(Data_Type = slot) %>% 
    dplyr::group_by(Cell) %>% 
    dplyr::mutate(ici_score_scale = scale(ici_score))
  
  if(!is.null(write_file)) readr::write_csv(ICI_scores, path = write_file)
  return(ici)
}

