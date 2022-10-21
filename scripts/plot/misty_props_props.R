# Plot misty props props results
library(mistyR)
future::plan(future::multisession)
library(tidyverse)

plot_misty_res <- function(lesion_type){

    meta <- read.csv('data/metadata.csv')
    sample_ids <- meta$sample_id[meta$lesion_type == lesion_type]

    # Find result directories
    res <- unlist(lapply(sample_ids, function(x){
        dirs <- sprintf("data/prc/visium/%s/props/", x)
    }))

    # Aggregate (mean) across samples
    res <- collect_results(res)

    # Plot results
    lesion_type <- stringr::str_replace_all(lesion_type, " ", "")
    res_dir <- sprintf("figures/misty_props_props/%s", lesion_type)
    dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

    # Plot gain R2
    plot_path <- sprintf("%s/R2.pdf", res_dir)
    pdf(file = plot_path, width = 4, height = 4)
    plot_improvement_stats(res, "multi.R2")
    dev.off()

    # Plot contributions
    plot_path <- sprintf("%s/contributions.pdf", res_dir)
    pdf(file = plot_path, width = 4, height = 4)
    plot_view_contributions(res)
    dev.off()

    # Plot int importances
    plot_path <- sprintf("%s/imp_int.pdf", res_dir)
    pdf(file = plot_path, width = 4, height = 4)
    plot_interaction_heatmap(res, view = 'intra', cutoff = 0)
    dev.off()

    # Plot juxta importances
    plot_path <- sprintf("%s/imp_jux.pdf", res_dir)
    pdf(file = plot_path, width = 4, height = 4)
    plot_interaction_heatmap(res, view = 'jux_5', cutoff = 0)
    dev.off()

    # Plot para importances
    plot_path <- sprintf("%s/imp_pra.pdf", res_dir)
    pdf(file = plot_path, width = 4, height = 4)
    plot_interaction_heatmap(res, view = 'pra_15', cutoff = 0)
    dev.off()
}

plot_misty_res('Chronic Active')
plot_misty_res('Control')
<<<<<<< HEAD

=======
>>>>>>> Added plotting scripts
