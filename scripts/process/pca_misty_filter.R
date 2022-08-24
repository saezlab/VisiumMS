# Extract interactions that separate conditions in both PC spaces

library(tidyverse)

lesions <- c('Chronic Active', 'Control')
view_name <- 'intra'
sign_type <- 'pos'

# Read PCS and loadings
lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')


ldns_props_props <- read.csv(sprintf("data/prc/pca_misty/props_props/ldns_%s_%s.csv",
                        lesions, view_name), row.names = 1)
ldns_sign_props <- read.csv(sprintf("data/prc/pca_misty/sign_props/ldns_%s_%s_%s.csv",
                        lesions, sign_type, 'int'), row.names = 1)

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

write.csv(inter, 'data/prc/pca_misty/inters.csv', row.names = F)

