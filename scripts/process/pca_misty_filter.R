# Extract interactions that separate conditions in both PC spaces
library(tidyverse)


filter_pcs <- function(lesions, view_name, sign_type){

    # Read PCS and loadings
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')

    if (view_name == 'int'){
        prop_name = 'intra'

    } else {
        prop_name = view_name
    }

    # Load
    pcs_props_props <- read.csv(sprintf("data/prc/pca_misty/props_props/pcs_%s_%s.csv",
                            lesions, prop_name))
    pcs_sign_props <- read.csv(sprintf("data/prc/pca_misty/sign_props/pcs_%s_%s_%s.csv",
                            lesions, sign_type, view_name))
    ldns_props_props <- read.csv(sprintf("data/prc/pca_misty/props_props/ldns_%s_%s.csv",
                            lesions, prop_name), row.names = 1)
    ldns_sign_props <- read.csv(sprintf("data/prc/pca_misty/sign_props/ldns_%s_%s_%s.csv",
                            lesions, sign_type, view_name), row.names = 1)

    # ANOVA PCs
    pval_props_props <- summary(aov(PC1 ~ Lesion.Type,
                                    data = pcs_props_props))[[1]][["Pr(>F)"]][1]
    pval_sign_props <- summary(aov(PC1 ~ Lesion.Type,
                                   data = pcs_sign_props))[[1]][["Pr(>F)"]][1]

    # Extract name
    corr_name <- glue::glue('corr_{view_name}_{sign_type}_props_props')
    pval_name <- glue::glue('pval_{view_name}_{sign_type}_props_props')
    props_props <- ldns_props_props[, 'PC1', drop=F] %>%
        mutate(!!pval_name := pval_props_props) %>%
        rename(!!corr_name := PC1) %>%
        as_tibble(rownames="CC")
    corr_name <- glue::glue('corr_{view_name}_{sign_type}_sign_props')
    pval_name <- glue::glue('pval_{view_name}_{sign_type}_sign_props')
    sign_props <- ldns_sign_props[, 'PC1', drop=F] %>%
        mutate(!!pval_name := pval_sign_props) %>%
        rename(!!corr_name := PC1) %>%
        as_tibble(rownames='CC')
    cc <- inner_join(props_props, sign_props, by='CC')

    props_props <- tibble(ints = props_props) %>%
        separate(ints, into=c('A','B'), sep='[|]') %>%
        rowwise() %>%
        mutate(inter = paste(sort(c(A, B)), collapse = "|")) %>%
        ungroup() %>%
        pull(inter) %>%
        unique()

    sign_props <- tibble(ints = sign_props) %>%
        separate(ints, into=c('A','B'), sep='[|]') %>%
        rowwise() %>%
        mutate(inter = paste(sort(c(A, B)), collapse = "|")) %>%
        ungroup() %>%
        pull(inter) %>%
        unique()

    inter <- intersect(props_props, sign_props)
    inter <- tibble(ints = inter) %>%
        separate(ints, into=c('A','B'), sep='[|]')

    res_path <- 'data/prc/pca_misty/inters'
    dir.create(res_path, showWarnings = FALSE, recursive = TRUE)
    write.csv(inter, sprintf('%s/%s_%s_%s.csv', res_path, lesions,
                             sign_type, view_name), row.names = F)
}


lesions <- c('Chronic Active', 'Control')

filter_pcs(lesions, 'int', 'pos')
filter_pcs(lesions, 'jux_5', 'pos')
filter_pcs(lesions, 'pra_15', 'pos')
filter_pcs(lesions, 'int', 'neg')
filter_pcs(lesions, 'jux_5', 'neg')
filter_pcs(lesions, 'pra_15', 'neg')
