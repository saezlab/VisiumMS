library(zellkonverter)
library(liana)
library(tidyverse)


# Read adata
adata <- readH5AD('data/prc/sc/raw.h5ad')
inters <- read.csv('data/prc/pca_misty/inters/ChronicActive|Control_pos_intra.csv')

# Filter by lesion types
msk <-  adata@colData$lesion_type %in% c('Chronic Active', 'Control')
adata <- adata[, msk]

# Extract unique cell types
cell_types <- unique(unname(unlist(inters)))
cell_types <- gsub('\\.', '_', cell_types)

# Filter by cell types
msk <-  adata@colData$leiden %in% c('Oligos', 'Astros', 'Microglia', 'Macrophages_f')#cell_types
adata <- adata[, msk]

# Rename assays
names(adata@assays) <- c("counts")

# Normalize
adata <- scuttle::logNormCounts(adata)

# Assign label
SingleCellExperiment::colLabels(adata) <- adata@colData$leiden


run_liana <- function(adata, lesion_type){

    print(lesion_type)
    # Filter by lesion_type
    msk <- adata@colData$lesion_type == lesion_type
    tmp <- adata[, msk]

    # Run liana
    res <- liana_aggregate(liana_wrap(tmp)) %>%
        filter(aggregate_rank < 0.05) %>%
        select(source, target, ligand.complex, receptor.complex,
               aggregate_rank, natmi.edge_specificity, sca.LRscore)

    # Write
    lesion_type <- stringr::str_replace_all(lesion_type, " ", "")
    res_path <- 'data/prc/liana/'
    dir.create(res_path, showWarnings = FALSE, recursive = TRUE)
    write.csv(res, sprintf('%s/%s.csv', res_path, lesion_type), row.names = FALSE)
}


run_liana(adata, 'Chronic Active')
run_liana(adata, 'Control')
