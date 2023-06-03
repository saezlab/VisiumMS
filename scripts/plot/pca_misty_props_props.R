library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(patchwork)


plot_misty_importances <- function(lesions, view_name){

    # Read PCS and loadings
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')
    pcs <- read.csv(sprintf("data/prc/pca_misty/props_props/pcs_%s_%s.csv",
                            lesions, view_name))
    ldns <- read.csv(sprintf("data/prc/pca_misty/props_props/ldns_%s_%s.csv",
                            lesions, view_name), row.names = 1)

    # ANOVA PCs
    aov(PC1 ~ Lesion.Type, data = pcs)
    pc1_pval <- -log10(summary(aov(PC1 ~ Lesion.Type, data = pcs))[[1]][["Pr(>F)"]][1])
    pc2_pval <- -log10(summary(aov(PC2 ~ Lesion.Type, data = pcs))[[1]][["Pr(>F)"]][1])
    pvals <- data.frame(row.names='Assoc', PC1=pc1_pval, PC2=pc2_pval)

    # Scatter
    pca_scatter <- ggplot(pcs, 
                          aes(x=PC1, y=PC2, label=Sample.ID, color=Lesion.Type)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, max.time=5, max.iter=1000000) +
    theme_bw()

    # Heatmaps
    h_pcs <- Heatmap(pcs[c('PC1', 'PC2')], cluster_columns = FALSE, name = 'PC scores',
                    right_annotation = rowAnnotation(Condition = pcs$Lesion.Type))
    h_ldns <- Heatmap(ldns, cluster_columns = FALSE, name = 'PC Loadings')
    h_pvals <- Heatmap(pvals, cluster_columns=F, name='P-value',
                      col = circlize::colorRamp2(c(0, 4), c("white", "orange3")))
    h_list <- h_ldns %v% h_pvals %v% h_pcs

    # Plot
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')
    dir.create("figures/pca_misty_props_props/", showWarnings = FALSE, recursive = TRUE)
    plot_path <- sprintf("figures/pca_misty_props_props/%s_%s.pdf", lesions, view_name)
    layout <- "
    ##BBBB
    AABBBB
    ##BBBB
    "
    plt <- pca_scatter + grid.grabExpr(draw(h_list)) +
    plot_layout(design = layout) + plot_annotation(
      title = view_name,
      theme = theme(plot.title = element_text(size = 18))
    )
    
    # Save
    pdf(file = plot_path, width = 12, height = 8)
    plot(plt)
    dev.off()
}

lesions <- c('Chronic Active', 'Control')
plot_misty_importances(lesions, 'intra')
plot_misty_importances(lesions, 'jux_5')
plot_misty_importances(lesions, 'pra_15')
