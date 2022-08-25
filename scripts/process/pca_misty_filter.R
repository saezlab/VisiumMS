# Extract interactions that separate conditions in both PC spaces
library(tidyverse)


filter_pcs <- function(lesions, view_name, sign_type){
    # Read PCS and loadings
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')

    if (view_name == 'intra'){
        prop_name = 'int'
    } else {
        prop_name = view_name
    }

    # Load
    ldns_props_props <- read.csv(sprintf("data/prc/pca_misty/props_props/ldns_%s_%s.csv",
                            lesions, view_name), row.names = 1)
    ldns_sign_props <- read.csv(sprintf("data/prc/pca_misty/sign_props/ldns_%s_%s_%s.csv",
                            lesions, sign_type, prop_name), row.names = 1)

    msk <- ldns_props_props$PC1 > 0.15
    props_props <- rownames(ldns_props_props[msk,])

    msk <- ldns_sign_props$PC1 > 0.15
    sign_props <- rownames(ldns_sign_props[msk,])


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
filter_pcs(lesions, 'intra', 'pos')
filter_pcs(lesions, 'jux_5', 'pos')
filter_pcs(lesions, 'pra_15', 'pos')
filter_pcs(lesions, 'intra', 'neg')
filter_pcs(lesions, 'jux_5', 'neg')
filter_pcs(lesions, 'pra_15', 'neg')

