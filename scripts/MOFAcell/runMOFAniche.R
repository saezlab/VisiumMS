# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we infer multiniche programs
#' using MOFAcell
library(MOFAcellulaR)
library(tidyverse)
library(compositions)

reticulate::use_condaenv(condaenv = "/Users/ricardoramirez/mambaforge/envs/sc_base")
pb_data_file <- "./submission/data/MOFAniche_pb_data.csv"
coldata_file <- "./submission/data/MOFAniche_pb_coldata.csv"
factors_file <- "./submission/data/MOFAniche_factors.csv"
loadings_file <- "./submission/data/MOFAniche_loadings.csv"
mofamodel_file <- file.path("./submission/data/", "mofaniche_complete.hdf5")
hmap_file <- "./submission/data/MOFAniche_hmap.pdf"
bplots_file <- "./submission/data/MOFAniche_bplots.pdf"
dist_file <- "./submission/data/MOFAniche_fspace.pdf"

# Importing pb data
pb_data <- read_csv(pb_data_file, show_col_types = FALSE) %>%
  as.data.frame()
rownames(pb_data) <- pb_data[,1]
pb_data <- pb_data[, -1]
pb_data <- as.matrix(pb_data) %>% t()

# Importing coldata of the matrices - Only MS
coldat <- read_csv(coldata_file,
                   show_col_types = FALSE)[,-1] %>%
  column_to_rownames("colname") %>%
  dplyr::rename(donor_id = "Sample id") %>%
  dplyr::filter(Condition == "MS")

pb_data <- pb_data[,rownames(coldat)]

meta <- coldat %>%
  dplyr::select_at(c("donor_id", "Condition", 
                     "Lesion type", "Age", "Sex", "Duration (years)")) %>%
  unique()

cts <- coldat$niches %>% unique()

# Do all processing of counts matrices
pb_dat_df <- MOFAcellulaR::create_init_exp(counts = pb_data[,rownames(coldat)],  
                                           coldata = coldat) %>%
  MOFAcellulaR::filt_profiles(pb_dat = .,
                              cts = cts,
                              ncells = 10,
                              counts_col = "psbulk_n_cells", 
                              ct_col = "niches") %>%
  MOFAcellulaR::filt_views_bysamples(pb_dat_list = .,
                                     nsamples = 9) %>%
  MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                min.count = 100,
                                min.prop = 0.25) %>%
  MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                   ngenes = 50) %>%
  MOFAcellulaR::filt_samples_bycov(pb_dat_list = ., # Filtering of low quality samples
                                   prop_coverage = 0.90) %>%
  MOFAcellulaR::tmm_trns(pb_dat_list = .,
                         scale_factor = 1000000) %>% 
  MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                   ngenes = 50) %>%
  MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                            sample_column = "donor_id")
# Run MOFA
MOFAobject <- create_mofa(pb_dat_df)

data_opts <- get_default_data_options(MOFAobject)

model_opts <- get_default_model_options(MOFAobject)

model_opts$num_factors <- 3

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
                                          sample_id_column = "donor_id") %>%
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
  dplyr::rename("Lesion" = "Lesion type") %>%
  dplyr::rename(duration = "Duration (years)")

expl_var_lession <- MOFAcellulaR::get_associations(model = model,
                                                   metadata = meta,
                                                   sample_id_column = "donor_id",
                                                   test_variable = "Lesion",
                                                   test_type = "categorical",
                                                   categorical_type = "parametric",
                                                   group = FALSE)

expl_var_sex <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "donor_id",
                                               test_variable = "Sex",
                                               test_type = "categorical",
                                               categorical_type = "parametric",
                                               group = FALSE)

expl_var_age <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "donor_id",
                                               test_variable = "Age",
                                               test_type = "continous",
                                               group = FALSE)

expl_var_duration <- MOFAcellulaR::get_associations(model = model,
                                               metadata = meta,
                                               sample_id_column = "donor_id",
                                               test_variable = "duration",
                                               test_type = "continous",
                                               group = FALSE)


assoc_list <- list(lession = expl_var_lession,
                   sex = expl_var_sex,
                   age = expl_var_age,
                   duration = expl_var_duration)

col_list <- list(Lesion = c("CA" = "#FF6666", 
                            "CI"= "#3CB371"))

scores_hmap <- MOFAcellulaR::plot_MOFA_hmap(model = model,
                                            group = FALSE,
                                            metadata = meta,
                                            sample_id_column = "donor_id",
                                            sample_anns = c("Lesion"),
                                            assoc_list = assoc_list,
                                            col_rows = col_list)


pdf(hmap_file, height = 4.5, width = 3.5)

draw(scores_hmap)

dev.off()

# Make a 2D plot with MOFA space
plt_2d <- MOFAcellulaR::plot_sample_2D(model, 
                                       method = "MDS", 
                                       metadata = meta, 
                                       sample_id_column = "donor_id",
                                       color_by = "Lesion")

plt_2d <- last_plot()

pdf(dist_file, width = 3.3, height = 2.4)

plot(plt_2d)

dev.off()

# Boxplots

bp1 <- MOFAcellulaR::get_tidy_factors(model, meta,factor = "Factor1",sample_id_column = "donor_id") %>%
  ggplot(aes(x = Age, y = value)) +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10)) +
  ylab("Factor1")


bp3 <- MOFAcellulaR::get_tidy_factors(model, meta,factor = "Factor3",sample_id_column = "donor_id") %>%
  ggplot(aes(x = Lesion, y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10)) +
  ylab("Factor3")

pdf(bplots_file, width = 3.5, height = 1.7)

cowplot::plot_grid(bp1, bp3,
                   nrow = 1, ncol = 2)

dev.off()





