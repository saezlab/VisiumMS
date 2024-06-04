library(rhdf5)
library(lme4)
library(MuMIn)
library(tidyverse)

read_factor_data <- function(x){
    codes <- c(x$codes)
    cats <- c(x$categories)
    return(as.character(factor(codes, levels = 0:(length(cats)-1), labels = cats)))
}

read_rna <- function(ctype, genes){
    path_data <- paste0('data/prc/ctypes/', ctype, '.h5ad')
    indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
    # RNA
    rna_X <- indata$X
    rna_X <- Matrix::sparseMatrix(
        i=indata$X$indices,
        p=indata$X$indptr,
        x=as.numeric(indata$X$data),
        index1 = FALSE
    )
    colnames(rna_X) <- indata$obs$`_index`
    rownames(rna_X) <- indata$var$`_index`
    rna_X <- as.data.frame(t(as.matrix(rna_X[genes, , drop=FALSE])))
    rna_X$condition <- read_factor_data(indata$obs$`Lesion type`)
    rna_X$patient <- read_factor_data(indata$obs$`Patient id`)
    h5closeAll()
    return(rna_X)
}

fit_lmm <- function(df){
    # Fit model
    model <- lmerTest::lmer(measurement ~ condition + (1 | patient), data = df)

    # Get fixed effect R2
    FE_R2 <- unname(r.squaredGLMM(model)[,"R2m"])

    # Get random effect R2
    model_res <- summary(model)
    RE_R2 <- model_res$varcor %>%
      as.data.frame() %>%
      dplyr::select(grp, vcov) %>%
      pivot_wider(names_from = grp, values_from = vcov) %>%
      dplyr::mutate(perc_studyvar = patient/(patient + Residual)) %>% # change study for 
      pull(perc_studyvar)

    # Get pvalue
    pval <- model_res$coefficients[2, 5]
    return(c(FE_R2, RE_R2, pval))
}

test_sign <- function(data, genes, contrast){
    res <- map(genes, function(gene){
        print(gene)
        df <- data[data$condition %in% contrast, c(gene, 'condition', 'patient')] %>% rename('measurement'=gene)
        fit_lmm(df)
    })
    res <- as.data.frame(do.call(rbind, res))
    colnames(res) <- c('fe', 're', 'pval')
    res$padj <- p.adjust(res$pval, method='BH')
    res$genes <- genes
    res$
    return(res)
}


contrast <- c('Ctrl', 'CA')

res = list()
genes <- c('HMGB1', 'ITGB1')
res$AS <- test_sign(read_rna('AS', genes=genes), genes=genes, contrast=contrast)
genes <- c('CD163', 'CD14', 'ITGB1', 'TLR2')
res$MG <- test_sign(read_rna('MG', genes=genes), genes=genes, contrast=contrast)
genes <- c('ITGB1')
res$EC <- test_sign(read_rna('EC', genes=genes), genes=genes, contrast=contrast)
res <- as.data.frame(do.call(rbind, res))
print(res)
