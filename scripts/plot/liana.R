brary(liana)
library(tidyverse)


plot_liana <- function(lesion_type, n = 5){
    # Read liana results
    lesion_type <- stringr::str_replace_all(lesion_type, " ", "")
    res_path <- 'data/prc/liana'
    res <- read.csv(sprintf("%s/%s.csv", res_path, lesion_type))

    # Select top interacitons per cluster
    top_inters <- res %>%
        group_by(source, target) %>%
        slice_min(aggregate_rank, n=5) %>%
        ungroup() %>%
        distinct_at(c('ligand.complex', 'receptor.complex'))

    # Plot
    plt <- liana_dotplot(inner_join(top_inters, res)) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    # Save
    plot_path <- 'figures/liana'
    dir.create(plot_path, showWarnings = FALSE, recursive = TRUE)
    pdf(file = sprintf('%s/%s.pdf', plot_path, lesion_type), width = 12, height = 12)
    plot(plt)
    dev.off()
}

plot_liana('Chronic Active', n=5)
plot_liana('Control', n=5)

