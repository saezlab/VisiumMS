library(mistyR)
future::plan(future::multisession)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(patchwork)




plot_misty_importances <- function(lesions, view_name){
    # Read meta and filter by lesions
    meta <- read.csv('data/metadata.csv')
    rownames(meta) <- meta$sample_id
    sample_ids <- meta$sample_id[meta$lesion_type %in% lesions]

    # Find result directories and load them
    res <- unlist(lapply(sample_ids, function(x){
        dirs <- sprintf("data/prc/visium/%s/props/", x)
    }))
    res <- collect_results(res)

    # Filter importances
    df <- res$importances
    df$sample_id <- unlist(lapply(df$sample, function(x){strsplit(x, "/")[[1]][9]}))
    df <- df %>% 
        filter(!is.na(Importance), .data$view == view_name) %>%
        mutate(feature=paste0(Predictor, "|",Target)) %>%
        pivot_wider(id_cols=sample_id, names_from=feature, values_from=Importance) %>%
        column_to_rownames('sample_id')

    # Compute PCA
    pca <- prcomp(df)
    pcs <- as.data.frame(pca$x[,1:2])
    pcs['Lesion Type'] <- meta[rownames(pcs), 'lesion_type']
    pcs['Sample ID'] <- rownames(pcs)

    # ANOVA PCs
    aov(PC1 ~ `Lesion Type`, data = pcs)
    pc1_pval <- -log10(summary(aov(PC1 ~ `Lesion Type`, data = pcs))[[1]][["Pr(>F)"]][1])
    pc2_pval <- -log10(summary(aov(PC2 ~ `Lesion Type`, data = pcs))[[1]][["Pr(>F)"]][1])
    pvals <- data.frame(row.names='Assoc', PC1=pc1_pval, PC2=pc2_pval)

    # Scatter
    pca_scatter <- ggplot(pcs, 
                          aes(x=PC1, y=PC2, label=`Sample ID`, color=`Lesion Type`)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, max.time=5, max.iter=1000000) +
    theme_bw()

    # Heatmaps
    h_pcs <- Heatmap(pcs[c('PC1', 'PC2')], cluster_columns = FALSE, name = 'PC scores',
                    right_annotation = rowAnnotation(Condition = pcs$`Lesion Type`))
    ldns <- as.data.frame(pca$rotation[,1:2]) %>%
        rownames_to_column('Interaction') %>%
        gather(key='PC', value='Loading', PC1, PC2) %>%
        filter(abs(Loading) > 0.15)
    ldns <- as.data.frame(pca$rotation[,1:2]) %>%
        filter(abs(PC1) > 0.15 |abs(PC2) > 0.15)
    h_ldns <- Heatmap(ldns, cluster_columns = FALSE, name = 'PC Loadings')

    h_pvals <- Heatmap(pvals, cluster_columns=F, name='P-value',
                      col = circlize::colorRamp2(c(0, 4), c("white", "orange3")))
    h_list <- h_ldns %v% h_pvals %v% h_pcs

    # Plot
    lesions <- stringr::str_c(stringr::str_replace_all(lesions, " ", ""), collapse='|')
    dir.create("figures/pca_misty/", showWarnings = FALSE, recursive = TRUE)
    plot_path <- sprintf("figures/pca_misty/%s_%s.pdf", lesions, view_name)
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
    theme = theme(plot.title = element_text(size = 18))
    
    # Save
    pdf(file = plot_path, width = 12, height = 8)
    plot(plt)
    dev.off()
}

lesions <- c('Chronic Active', 'Control')
plot_misty_importances(lesions, 'intra')
plot_misty_importances(lesions, 'jux_5')
plot_misty_importances(lesions, 'pra_15')

