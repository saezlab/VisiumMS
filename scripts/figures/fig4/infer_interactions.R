library(zellkonverter)
library(liana)
library(tidyverse)


# Read adata
adata <- readH5AD('data/prc/sc/raw.h5ad')
adata@colData$leiden <- as.character(adata@colData$leiden)

# Filter by lesion types
msk <-  adata@colData$lesion_type %in% c('Chronic Active')
adata <- adata[, msk]

# Filter by cell types
msk <-  adata@colData$leiden %in% c('Endothelia', 'Macrophages_f', 'Astros', 'Astros_c')
adata <- adata[, msk]


# Update states
states <- readH5AD('data/prc/sc/astros.h5ad')
msk <-  states@colData$lesion_type %in% c('Chronic Active')
states <- states[, msk]
state <- '2'
msk <- states@colData$leiden == state
barcodes <- rownames(states@colData)[msk]
adata@colData[barcodes, 'leiden'] <- state
state <- '5'
msk <- states@colData$leiden == state
barcodes <- rownames(states@colData)[msk]
adata@colData[barcodes, 'leiden'] <- state

# Subset
msk <-  adata@colData$leiden %in% c('2', '5', 'Macrophages_f', 'Endothelia')
adata <- adata[, msk]

# Rename assays
names(adata@assays) <- c("counts")

# Normalize
adata <- scuttle::logNormCounts(adata)

# Assign label
SingleCellExperiment::colLabels(adata) <- adata@colData$leiden

# Run liana
res <- liana_aggregate(liana_wrap(adata)) %>%
    filter(aggregate_rank < 0.05) %>%
    select(source, target, ligand.complex, receptor.complex,
           aggregate_rank, natmi.edge_specificity, sca.LRscore)

# Write
res_path <- 'figures/manuscript/fig4/'
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf('%s/astros_inters.csv', res_path), row.names = FALSE)
