# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate markers of cells
#' using edgeR and pseudobulk profiles of all samples

library(SingleCellExperiment)
library(scater)
library(edgeR)
library(tidyverse)

# Parse args
args <- commandArgs(trailingOnly = F)
pb_data_file <- args[6]
coldata_file <- args[7]
out_file <- args[8]

# Importing pb data
pb_data <- read_csv(pb_data_file, show_col_types = FALSE) %>%
  as.data.frame()
rownames(pb_data) <- gsub("_c", "c", pb_data[,1])
pb_data <- pb_data[, -1]
pb_data <- as.matrix(pb_data) %>% t()

# Importing coldata of the matrices
coldat <- read_csv(coldata_file,
                   show_col_types = FALSE)[,-1] %>%
  dplyr::mutate(colname = gsub("_c", "c", colname)) %>%
  dplyr::rename(cell_types = "leiden") %>%
  dplyr::mutate(cell_types = gsub("_c", "c", cell_types)) %>%
  column_to_rownames("colname") %>%
  dplyr::rename(cell_counts = "counts")


pb_data <- pb_data[,rownames(coldat)]

# Defining cts
cts <- coldat$cell_types %>% 
  unique() %>%
  set_names()

# Pipeline for differential expression

de_res <- map(cts, function(ct) {
  print(ct)
  ct_meta_data <- coldat %>%
    mutate(test_column = ifelse(cell_types == ct, ct, "rest"))
  
  dat <- DGEList(pb_data, samples = DataFrame(ct_meta_data))
  
  keep <- filterByExpr(dat, group = ct_meta_data$test_column)
  
  dat <- dat[keep,]
  
  dat <- calcNormFactors(dat)
  
  design <- model.matrix(~factor(test_column,
                                 levels = c("rest",ct)), dat$samples)
  
  colnames(design) <- c("int", ct)
  
  dat <- estimateDisp(dat, design)
  
  fit <- glmQLFit(dat, design, robust=TRUE)
  
  res <- glmQLFTest(fit, coef=ncol(design))
  
  de_res <- topTags(res, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  return(de_res)
  
})

de_res <- de_res %>% 
  enframe() %>%
  unnest()

de_res %>%
  dplyr::filter(logFC > 0) %>%
  arrange(name, FDR, - logFC) %>%
  write_csv(file = out_file)
