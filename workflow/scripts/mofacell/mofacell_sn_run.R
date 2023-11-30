# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we infer multicellular programs
#' using MOFAcell

library(MOFAcellulaR)
library(tidyverse)
library(compositions)

# Parse args
args <- commandArgs(trailingOnly = F)
pb_data_file <- args[6]
coldata_file <- args[7]
meta_file <- args[8]
mrkrs_file <- args[9]
factors_file <- args[10]
loadings_file <- args[11]
mofamodel_file <- args[12]
hmap_file <- args[13]
dist_file <- args[14]
bplots_file <- args[15]
qc_file <- args[16]

# Importing pb data
pb_data <- read_csv(pb_data_file, show_col_types = FALSE) %>%
  as.data.frame()
rownames(pb_data) <- gsub("_c", "c", pb_data[,1])
pb_data <- pb_data[, -1]
pb_data <- as.matrix(pb_data) %>% t()

# Importing coldata of the matrices
coldat <- read_csv(coldata_file,
                   show_col_types = FALSE)[,-1] %>%
  dplyr::mutate(colname = gsub("_c", "c", colname)) %>%
  dplyr::mutate(cell_types = gsub("_c", "c", leiden)) %>%
  column_to_rownames("colname") %>%
  dplyr::rename(cell_counts = "counts") %>%
  dplyr::rename(donor_id = "Sample id")

pb_data <- pb_data[,rownames(coldat)]

meta <- read_csv(meta_file)[, -1]

# Add mrker genes for background
mrkr_genes <- read_csv(mrkrs_file)

mrkr_genes <- mrkr_genes %>% #dplyr::filter(!name %in% exclude_ct) %>%
  dplyr::filter(FDR < 0.01, logFC > 1) %>%
  dplyr::select(name, gene) %>%
  dplyr::rename("lineage" = name) %>%
  group_by(lineage) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

# 
cts <- coldat$leiden %>% unique()

# Do all processing of counts matrices
pb_dat_df <- MOFAcellulaR::create_init_exp(counts = pb_data[,rownames(coldat)],  
                                           coldata = coldat) %>%
  MOFAcellulaR::filt_profiles(pb_dat = .,
                              cts = cts,
                              ncells = 10,
                              counts_col = "cell_counts", 
                              ct_col = "cell_types") %>%
  MOFAcellulaR::filt_views_bysamples(pb_dat_list = .,
                                     nsamples = 10) %>%
  MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                min.count = 100,
                                min.prop = 0.25) %>%
  MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                   ngenes = 50) %>%
  MOFAcellulaR::filt_samples_bycov(pb_dat_list = ., # Filtering of low quality samples
                                   prop_coverage = 0.90) %>%
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>% 
  MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = .,
                                   prior_mrks = mrkr_genes) %>%
  MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                   ngenes = 50) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .)


# Run MOFA

MOFAobject <- create_mofa(pb_dat_df)

data_opts <- get_default_data_options(MOFAobject)

model_opts <- get_default_model_options(MOFAobject)

model_opts$num_factors <- 4

model_opts$spikeslab_weights <- FALSE

train_opts <- get_default_training_options(MOFAobject)

# Prepare MOFA model:
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train model:
model <- run_mofa(MOFAobject, mofamodel_file)

# Model Factors
factors <- MOFAcellulaR::get_tidy_factors(model, 
                                          factor = "all", 
                                          metadata = meta, 
                                          sample_id_column = "Sample id") %>%
  dplyr::rename("Sample id" = "sample")

write.csv(factors, factors_file)


# Model loadings

factor_loadings <- get_factors(model)[[1]] %>% colnames() %>%
  set_names() %>%
  map(., ~ get_geneweights(model = model, factor = .x)) %>%
  enframe(name = "Factor") %>%
  unnest()

write.csv(factor_loadings, loadings_file)

# Check completeness (plot)
complete_plt <- pb_dat_df %>% 
  dplyr::select(view, sample) %>% 
  unique() %>% 
  group_by(sample) %>% 
  summarize(n_views = n()) %>%
  ggplot(aes(x = n_views, y = sample)) +
  geom_bar(stat="identity") +
  theme_classic()

pdf(qc_file, height = 3, width = 3)

plot(complete_plt)

dev.off()

# Test the associations
# First do data manipulation
meta <- meta %>%
  dplyr::rename("Lesion" = "Lesion type",
                "Batch" = "Batch sn")


expl_var_dis <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "Sample id",
                                               test_variable = "Condition",
                                               test_type = "categorical",
                                               categorical_type = "parametric",
                                               group = FALSE)

expl_var_lession <- MOFAcellulaR::get_associations(model = model,
                                                   metadata = meta,
                                                   sample_id_column = "Sample id",
                                                   test_variable = "Lesion",
                                                   test_type = "categorical",
                                                   categorical_type = "parametric",
                                                   group = FALSE)

expl_var_sex <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "Sample id",
                                               test_variable = "Sex",
                                               test_type = "categorical",
                                               categorical_type = "parametric",
                                               group = FALSE)

expl_var_age <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "Sample id",
                                               test_variable = "Age",
                                               test_type = "continous",
                                               group = FALSE)

expl_var_rin <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "Sample id",
                                               test_variable = "RIN",
                                               test_type = "continous",
                                               group = FALSE)

expl_var_batch <- MOFAcellulaR::get_associations(model = model,
                                                 metadata = meta,
                                                 sample_id_column = "Sample id",
                                                 test_variable = "Batch",
                                                 test_type = "categorical",
                                                 categorical_type = "parametric",
                                                 group = FALSE)

# Make complex heatmap

assoc_list <- list(disease = expl_var_dis, 
                   lession = expl_var_lession,
                   sex = expl_var_sex,
                   age = expl_var_age,
                   rin = expl_var_rin,
                   batch = expl_var_batch)

col_list <- list(Lesion = c("CA" = "#FF6666", 
                            "CI"= "#3CB371",
                            "Ctrl" = "pink"),
                 Condition = c("Control" = "black",
                               "MS" = "darkgrey"))

scores_hmap <- MOFAcellulaR::plot_MOFA_hmap(model = model,
                                            group = FALSE,
                                            metadata = meta,
                                            sample_id_column = "Sample id",
                                            sample_anns = c("Lesion", "Condition"),
                                            assoc_list = assoc_list,
                                            col_rows = col_list)


pdf(hmap_file, height = 6.5, width = 4.0)

draw(scores_hmap)

dev.off()

# Make a 2D plot with MOFA space
plt_2d <- MOFAcellulaR::plot_sample_2D(model, 
                             method = "MDS", 
                             metadata = meta, 
                             sample_id_column = "Sample id",
                             color_by = "Lesion")

plt_2d <- last_plot()

pdf(dist_file, width = 3.3, height = 2.4)

plot(plt_2d)

dev.off()

# Boxplots

bp1 <- MOFAcellulaR::get_tidy_factors(model, meta,factor = "Factor1",sample_id_column = "Sample id") %>%
  ggplot(aes(x = Lesion, y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10)) +
  ylab("Factor1")

bp4 <- MOFAcellulaR::get_tidy_factors(model, meta,factor = "Factor4",sample_id_column = "Sample id") %>%
  ggplot(aes(x = Lesion, y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10)) +
  ylab("Factor4")

pdf(bplots_file, width = 4.5, height = 1.7)

cowplot::plot_grid(bp1, bp4,
                   nrow = 1, ncol = 3)

dev.off()



