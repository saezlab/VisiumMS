<<<<<<< HEAD
library(mistyR)
future::plan(future::multisession)
library(tidyverse)


pca_misty_importances <- function(lesions, view_name){
    # Read meta and filter by lesions
    meta <- read.csv('data/metadata.csv')
    rownames(meta) <- meta$sample_id
    sample_ids <- meta$sample_id[meta$lesion_type %in% lesions]

    # Find result directories and load them
    res <- unlist(lapply(sample_ids, function(x){
        dirs <- sprintf("data/prc/visium/%s/props/", x)
    }))
    res <- collect_results(res)

    # Filter importances
    df <- res$importances
    df$sample_id <- unlist(lapply(df$sample, function(x){strsplit(x, "/")[[1]][9]}))
    df <- df %>% 
        filter(!is.na(Importance), .data$view == view_name) %>%
        mutate(feature=paste0(Predictor, "|",Target)) %>%
        pivot_wider(id_cols=sample_id, names_from=feature, values_from=Importance) %>%
        column_to_rownames('sample_id')

    # Compute PCA
    pca <- prcomp(df)
    pcs <- as.data.frame(pca$x[,1:2])
    pcs['Lesion Type'] <- meta[rownames(pcs), 'lesion_type']
    pcs['Sample ID'] <- rownames(pcs)

    # Extract loadings
    ldns <- as.data.frame(pca$rotation[,1:2]) %>%
        filter(abs(PC1) > 0.15 |abs(PC2) > 0.15)

    # Write
    dir_path <- "data/prc/pca_misty/props_props"
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')
    write.csv(pcs, sprintf('%s/pcs_%s_%s.csv', dir_path, lesions, view_name), row.names=FALSE)
    write.csv(ldns, sprintf('%s/ldns_%s_%s.csv', dir_path, lesions, view_name))
}

lesions <- c('Chronic Active', 'Control')
pca_misty_importances(lesions, 'intra')
pca_misty_importances(lesions, 'jux_5')
pca_misty_importances(lesions, 'pra_15')

=======
# MISTy cell abundances model
library(Seurat)
library(mistyR)
future::plan(future::multisession)
library(tidyverse)
source('scripts/process/misty_utils.R')


read_slide <- function(sample_id){

    # Create paths
    slide_path <- sprintf('data/raw/visium/%s/outs/', sample_id)
    props_path <- sprintf('data/prc/visium/%s/cell_props.csv', sample_id)

    # Load data
    slide <- Load10X_Spatial(slide_path)
    slide$sample_id <- sample_id
    props <- t(compositions::clr(compositions::clo(read.csv(props_path, row.names=1))))

    # Intersect
    inter <- intersect(colnames(props), colnames(slide))
    slide <- slide[, inter]
    props <- props[, inter]

    # Create assays
    slide[["props"]] <- CreateAssayObject(counts = props)
    DefaultAssay(object = slide) <- "props"
    slide[['Spatial']] <- NULL

    return(slide)
}


run_misty_props <- function(slide){
    
    
    # Define assay of each view ---------------
    view_assays <- list("main" = 'props',
                        "jux" = 'props',
                        "pra" = 'props'
                       )
    # Define features of each view ------------
    view_features <- list("main" = NULL,
                          "jux" = NULL,
                          "pra" = NULL
                         )
    # Define spatial context of each view -----
    view_types <- list("main" = "intra", 
                       "jux" = "juxta",
                       "pra" = "para"
                      )
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "jux" = 5,
                        "pra" = 15)

    # Define name
    sample_id <- unique(slide$sample_id)
    misty_out <- sprintf('data/prc/visium/%s/props/', sample_id)

    run_misty_seurat(visium.slide = slide,
                     # Seurat object with spatial transcriptomics data.
                     view.assays = view_assays,
                     # Named list of assays for each view.
                     view.features = view_features,
                     # Named list of features/markers to use.
                     # Use all by default.
                     view.types = view_types,
                     # Named list of the type of view to construct
                     # from the assay.
                     view.params = view_params,
                     # Named list with parameters (NULL or value)
                     # for each view.
                     spot.ids = NULL,
                     # spot IDs to use. Use all by default.
                     out.alias = misty_out,
                     # bypass.intra
                     bypass.intra = FALSE
                    )
}

# Get name slides
meta <- read.csv('data/metadata.csv')
sample_ids <- meta$sample_id

for (sample_id in sample_ids){

    # Read slide
    slide <- read_slide(sample_id)

    # Run misty model
    run_misty_props(slide)
}
>>>>>>> Added plotting scripts
