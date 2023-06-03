# MISTy
library(Seurat)
library(mistyR)
future::plan(future::multisession)
library(tidyverse)
source('scripts/process/misty_utils.R')


read_slide <- function(sample_id){

    # Create paths
    slide_path <- sprintf('data/raw/visium/%s/outs/', sample_id)
    props_path <- sprintf('data/prc/visium/%s/cell_props.csv', sample_id)
    abnds_path <- sprintf('data/prc/visium/%s/cell_abunds.csv', sample_id)
    sign_pos_path <- sprintf('data/prc/visium/%s/cell_pos_sign.csv', sample_id)
    sign_neg_path <- sprintf('data/prc/visium/%s/cell_neg_sign.csv', sample_id)

    # Load data
    slide <- Load10X_Spatial(slide_path)
    slide$sample_id <- sample_id
    props <- t(compositions::clr(compositions::clo(read.csv(props_path, row.names=1))))
    abnds <- t(read.csv(abnds_path, row.names=1))
    sign_pos <- t(read.csv(sign_pos_path, row.names=1))
    sign_neg <- t(read.csv(sign_neg_path, row.names=1))

    # Intersect
    inter <- intersect(colnames(props), colnames(slide))
    inter <- intersect(inter, colnames(abnds))
    inter <- intersect(inter, colnames(sign_pos))
    inter <- intersect(inter, colnames(sign_neg))
    slide <- slide[, inter]
    props <- props[, inter]
    abnds <- abnds[, inter]
    sign_pos <- sign_pos[, inter]
    sign_neg <- sign_neg[, inter]

    # Create assays
    slide[["props"]] <- CreateAssayObject(counts = props)
    slide[["abnds"]] <- CreateAssayObject(counts = abnds)
    DefaultAssay(object = slide) <- "props"
    slide[['Spatial']] <- NULL
    for (sign in rownames(sign_pos)){
        tmp <- sign_pos[sign,, drop=F]
        slide[[paste0('pos.', sign)]] <- CreateAssayObject(counts = tmp)
    }
    for (sign in rownames(sign_pos)){
        tmp <- sign_neg[sign,, drop=F]
        slide[[paste0('neg.', sign)]] <- CreateAssayObject(counts = tmp)
    }

    return(slide)
}


run_misty_sign <- function(slide, sign){

    # Extract sign
    sign_type <- substr(sign, 1, 3)
    cell_name <- substr(sign, 5, 100)

    # Subset by proporitons
    min_prop <- 1 / nrow(slide@assays$props@data)
    msk <- slide@assays$props@data[cell_name,] > min_prop
    spot.ids <- colnames(slide@assays$props@data)[msk]
    
    
    # Define assay of each view ---------------
    view_assays <- list("main" = sign,
                        "int" = 'props',
                        "jux" = 'props',
                        "pra" = 'props'
                       )
    # Define features of each view ------------
    view_features <- list("main" = NULL,
                          "int" = NULL,
                          "jux" = NULL,
                          "pra" = NULL
                         )
    # Define spatial context of each view -----
    view_types <- list("main" = "intra", 
                       "int" = "intra",
                       "jux" = "juxta",
                       "pra" = "para"
                      )
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "int" = NULL,
                        "jux" = 5,
                        "pra" = 15)

    # Define name
    sample_id <- unique(slide$sample_id)
    misty_out <- sprintf('data/prc/visium/%s/signs/%s/%s',
                         sample_id, cell_name, sign_type)

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
                     spot.ids = spot.ids,
                     # spot IDs to use. Use all by default.
                     out.alias = misty_out,
                     # bypass.intra
                     bypass.intra = TRUE
                    )
}

# Get name slides
meta <- read.csv('data/metadata.csv')
sample_ids <- meta$sample_id


for (sample_id in sample_ids){

    # Read slide
    slide <- read_slide(sample_id)

    # Extract assay names
    assay_names <- names(slide@assays)
    msk <- grepl('^pos', assay_names) | grepl('^neg', assay_names)

    # Run misty models
    for (sign in assay_names[msk]){
        print(paste0(sample_id, ': ', sign))
        run_misty_sign(slide, sign)
    }
}
<<<<<<< HEAD

=======
>>>>>>> Added plotting scripts
