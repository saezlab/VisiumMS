library(mistyR)
future::plan(future::multisession)
library(tidyverse)


read_importances <- function(lesions, type){
    
    # Read meta and filter by lesions
    meta <- read.csv('data/metadata.csv')
    rownames(meta) <- meta$sample_id
    sample_ids <- meta$sample_id[meta$lesion_type %in% lesions]
    
    if (type == 'props'){
        # Find result directories and load them
        res <- unlist(lapply(sample_ids, function(x){
            dirs <- sprintf("data/prc/visium/%s/props/", x)
        }))
        res <- collect_results(res)

        # Filter importances
        df <- res$importances
        df$sample_id <- unlist(lapply(df$sample, function(x){strsplit(x, "/")[[1]][9]}))
        df <- df %>% 
            filter(!is.na(Importance)) %>%
            mutate(view = if_else(view == 'intra', 'int', view)) %>%
            mutate(feature=paste0("prop", "|", view, "|", Predictor, "|", Target)) %>%
            pivot_wider(id_cols=sample_id, names_from=feature, values_from=Importance) %>%
            column_to_rownames('sample_id') %>%
            as_tibble(rownames='sample_id')
    } else if (type == 'sign'){
        # pos
        # Find result directories and load them
        res <- unlist(lapply(sample_ids, function(x){
            dirs <- list.dirs(path = sprintf("data/prc/visium/%s/signs/", x),
                              full.names = TRUE, recursive = TRUE)
            dirs <- dirs[grepl(sprintf('%s$', 'pos'), dirs)]
        }))
        res <- collect_results(res)

        # Filter importances
        pos <- res$importances
        pos$sample_id <- unlist(lapply(pos$sample, function(x){strsplit(x, "/")[[1]][9]}))
        pos <- pos %>% 
            filter(!is.na(Importance), view != 'intra') %>%
            mutate(feature=paste0("pos", "|", view, "|", Predictor, "|", Target)) %>%
            pivot_wider(id_cols=sample_id, names_from=feature, values_from=Importance) %>%
            column_to_rownames('sample_id') %>%
            as_tibble(rownames='sample_id')
        
        # neg
        # Find result directories and load them
        res <- unlist(lapply(sample_ids, function(x){
            dirs <- list.dirs(path = sprintf("data/prc/visium/%s/signs/", x),
                              full.names = TRUE, recursive = TRUE)
            dirs <- dirs[grepl(sprintf('%s$', 'neg'), dirs)]
        }))
        res <- collect_results(res)

        # Filter importances
        neg <- res$importances
        neg$sample_id <- unlist(lapply(neg$sample, function(x){strsplit(x, "/")[[1]][9]}))
        neg <- neg %>% 
            filter(!is.na(Importance), view != 'intra') %>%
            mutate(feature=paste0("neg", "|", view, "|", Predictor, "|", Target)) %>%
            pivot_wider(id_cols=sample_id, names_from=feature, values_from=Importance) %>%
            column_to_rownames('sample_id') %>%
            as_tibble(rownames='sample_id')
        df <- inner_join(pos, neg)
    }
    return(df)
}


# Read meta and filter by lesions
meta <- read.csv('data/metadata.csv')
rownames(meta) <- meta$sample_id

# Read importances
lesions <- c('Chronic Active', 'Control')
props <- read_importances(lesions, 'props')
signs <- read_importances(lesions, 'sign')
df <- inner_join(props, signs)

# Compute PCA
pca <- prcomp(df %>% column_to_rownames('sample_id'))
pcs <- as.data.frame(pca$x[,1:2])
pcs['Lesion Type'] <- meta[rownames(pcs), 'lesion_type']
pcs['Sample ID'] <- rownames(pcs)

# Extract loadings
ldns <- as.data.frame(pca$rotation[,1:2]) %>%
    filter(abs(PC1) > 0.10)

# Write
dir_path <- "data/prc/pca_misty/joint"
dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')
write.csv(pcs, sprintf('%s/pcs_%s_%s.csv', dir_path, lesions, 'joint'), row.names=FALSE)
write.csv(ldns, sprintf('%s/ldns_%s_%s.csv', dir_path, lesions, 'joint'))
