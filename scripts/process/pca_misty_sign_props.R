library(mistyR)
future::plan(future::multisession)
library(tidyverse)


pca_misty_importances <- function(lesions, view_name, sign_type){
    # Read meta and filter by lesions
    meta <- read.csv('data/metadata.csv')
    rownames(meta) <- meta$sample_id
    sample_ids <- meta$sample_id[meta$lesion_type %in% lesions]

    # Find result directories and load them
    res <- unlist(lapply(sample_ids, function(x){
        dirs <- list.dirs(path = sprintf("data/prc/visium/%s/signs/", x),
                          full.names = TRUE, recursive = TRUE)
        dirs <- dirs[grepl(sprintf('%s$', sign_type), dirs)]
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
    dir_path <- "data/prc/pca_misty/sign_props"
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')
    write.csv(pcs, sprintf('%s/pcs_%s_%s_%s.csv', dir_path, lesions,
                           sign_type, view_name), row.names=FALSE)
    write.csv(ldns, sprintf('%s/ldns_%s_%s_%s.csv', dir_path, lesions,
                            sign_type, view_name))
}

lesions <- c('Chronic Active', 'Control')
pca_misty_importances(lesions, 'int', 'pos')
pca_misty_importances(lesions, 'jux_5', 'pos')
pca_misty_importances(lesions, 'pra_15', 'pos')
pca_misty_importances(lesions, 'int', 'neg')
pca_misty_importances(lesions, 'jux_5', 'neg')
pca_misty_importances(lesions, 'pra_15', 'neg')
