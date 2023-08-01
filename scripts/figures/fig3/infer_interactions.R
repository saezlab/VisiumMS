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
msk <-  adata@colData$leiden %in% c('Oligos', 'Microglia', 'Macrophages_f')
adata <- adata[, msk]


# Update states
states <- readH5AD('data/prc/sc/microglia.h5ad')
msk <-  states@colData$lesion_type %in% c('Chronic Active')
states <- states[, msk]
state <- '2'
msk <- states@colData$leiden == state
barcodes <- rownames(states@colData)[msk]
adata@colData[barcodes, 'leiden'] <- state
state <- '3'
msk <- states@colData$leiden == state
barcodes <- rownames(states@colData)[msk]
adata@colData[barcodes, 'leiden'] <- state

# Subset
msk <-  adata@colData$leiden %in% c('Oligos', '2', '3', 'Macrophages_f')
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
res_path <- 'figures/manuscript/fig3/'
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf('%s/microglia_inters.csv', res_path), row.names = FALSE)
