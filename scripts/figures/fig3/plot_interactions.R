library(tidyverse)
library(liana)

res <- as.tibble(read.csv('figures/manuscript/fig3/microglia_inters.csv')) %>%
    unite(inter, c("ligand.complex", "receptor.complex"), sep = "|", remove = FALSE)
inters <- c('HLA-A|APLP2', 'APP|APLP2', 'APOE|LRP1', 'SPP1|ITGA5', 'APOE|TREM2', 'C3|LRP1')
res <- res %>%
    filter(inter %in% inters)

pdf(file = "figures/manuscript/fig3/interactions.pdf",
    width = 16,
    height = 8)
res %>%
  liana_dotplot(source_groups = c("Macrophages_f", "2", "3", "Oligos"),
                target_groups = c("Macrophages_f", "2", "3", "Oligos"),
                ntop = 5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
