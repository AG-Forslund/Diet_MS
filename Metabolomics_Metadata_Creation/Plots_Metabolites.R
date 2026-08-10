make.correlation.plot_p.Value <- function(corr_data.frame, filter = TRUE, text_size = 7) {
  corr_df <- corr_data.frame %>%
    cor_gather() %>%
    mutate(star = ifelse(cor <= 0.05, "*", ""))
  
  if(filter) {
    sig.cor <- corr_df %>%
      filter(cor < 0.1)
    corr_df <- corr_df %>%
      filter(var1 %in% sig.cor$var1) %>%
      filter(var2 %in% sig.cor$var2)
  }
  
  p <- ggplot(data = corr_df, aes(var1, var2, fill = cor)) +
    geom_tile() +
    geom_text(aes(label=star), size = text_size) +
    scale_fill_distiller(palette = "PuBuGn", name = "P-Values", limits = c(0,1)) +
    theme_classic() +
    ggtitle("Correlation p-Values", subtitle="FDR-adjusted p-values: < 0.1 = * ") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = text_size),
          axis.text.y = element_text(size = text_size),
          axis.line = element_blank(),
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  return(p)
}



make.correlation.plot_cor.estimate <- function(corr_data.frame, cor.lim = 0.3, text_size = 7, filter = TRUE, scale_limits = c(-0.5,0.5)) {
  corr_df <- corr_data.frame %>%
    cor_gather() %>%
    mutate(star = ifelse(abs(cor) >= cor.lim, "+", ""))
  
  
  if(filter) {
    sig.cor <- corr_df %>%
      filter(star == "+")
    corr_df <- corr_df %>%
      filter(var1 %in% sig.cor$var1) %>%
      filter(var2 %in% sig.cor$var2)
  }
  
  corr_df %>%
    ggplot(aes(var2,var1, fill = cor)) +
    geom_tile() +
    geom_text(aes(label=star), size = 3) +
    scale_fill_distiller(palette = "PRGn", name = "Correlation\nEffect Size", limits = scale_limits) +
    theme_classic() +
    ggtitle("Correlation Estimates", subtitle=paste("Correlation", paste(paste(">", cor.lim)," = + ", sep = ""))) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = text_size),
          axis.text.y = element_text(size = text_size),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}



make_model_boxplot <- function(pdata, metadata, var, plot.kind = "by.sample") {
  col.class <- check_col_kind(metadata[,var])
  
  pdata_merged <- pdata %>%
    merge(., metadata[c(var, meta.merge)], by = meta.merge)

  if(col.class == "Binary") {
    metadata[,var] <- as.factor(metadata[,var])
  }
  
  if(col.class == "Binary" | col.class == "Factor") {
    levelsize <-  length(unique(metadata[,var]))
    if(plot.kind == "by.variable" | plot.kind == "by.sample") {
      if (plot.kind == "by.variable") {
        p <- pdata_merged %>%
          ggplot(aes_string(x = var, y = metric, color = var))
      } else if (plot.kind == "by.sample") {
        p <- pdata_merged %>%
          ggplot(aes_string(x = paste("fct_reorder(",paste(meta.merge, metric, sep = ", "), ", na.rm=TRUE)", sep = ""), y = metric, group = meta.merge, color = var))
      }
      plot <- p +
        geom_boxplot(outlier.size = 0.5) +
        theme_classic() +
        xlab(var) +
        theme(axis.text.x = element_blank()) +
        facet_grid(as.formula(paste("~", var)), scale="free", space = "free")
    } else if (plot.kind == "by.metabolite") {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = comp.name, y = metric, color = var)) +
        geom_boxplot(position=position_dodge(.9), outlier.size = 0.5) +
        geom_point(position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),alpha=0.3) +
        theme_classic() +
        scale_color_manual(values = safe_colorblind_palette[1:levelsize]) +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    }
  } else {
    if (plot.kind == "by.variable" | plot.kind == "by.sample") {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = paste("reorder(",paste(meta.merge, var, sep = ","), ")", sep = ""), y = metric, group = meta.merge, color = var)) +
        geom_boxplot(outlier.size = 0.5) +
        xlab(var) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))
    } else {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = comp.name, y = metric, group = comp.name)) +
        geom_boxplot(outlier.size = 0.5) +
        geom_point(aes_string(color = var)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    }
  }
    return(plot)
}



make_variable_boxplot <- function(pdata, metadata, var, metabolite, plot.kind = "by.sample") {
  col.class <- check_col_kind(metadata[,var])
  
  if(col.class == "Binary") {
    metadata[,var] <- as.factor(metadata[,var])
  }
  
  pdata_merged <- pdata %>%
    filter(Compound == metabolite) %>%
    merge(., metadata[c(var, meta.merge)], by = meta.merge)
  
  if(col.class == "Binary" | col.class == "Factor") {
    levelsize <- length(unique(metadata[,var]))
    if (plot.kind == "by.sample") {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = paste("reorder(",paste(meta.merge, var, sep = ","), ")", sep = ""), y = metric, group = meta.merge, color = var)) +
        geom_point() +
        xlab(var) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))
    } else if (plot.kind == "by.variable") {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = var, y = metric)) +
        geom_boxplot(aes_string(color = var), outlier.size = 0.5) +
        geom_point(aes_string(color = var)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    }
  } else {
    if (plot.kind == "by.sample") {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = paste("reorder(",paste(meta.merge, var, sep = ","), ")", sep = ""), y = metric, group = meta.merge, color = var)) +
        geom_point() +
        xlab(var) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))
    } else if (plot.kind == "by.variable") {
      plot <- pdata_merged %>%
        ggplot(aes_string(x = var, y = metric)) +
        geom_smooth(method = "lm") +
        geom_point(aes_string(color = var)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    }
  }
  
  return(plot)
}



post_hoc_effectsize <- function(d_results, cor.lim = 0.33, text_size = 7){
  d_results %>%
    gather(Feature, CD, -Comparison) %>%
    mutate(star = ifelse(abs(CD) >= cor.lim, "+", "")) %>%
    ggplot(aes(Comparison, Feature, fill = CD)) +
    geom_tile() +
    geom_text(aes(label=star), size = 3) +
    scale_fill_distiller(palette = "PRGn", name = "Cliff's Delta", limits = c(-1, 1)) +
    theme_classic() +
    ggtitle("Group post hoc effect sizes", subtitle=paste("Cliff's Delta", paste(paste(">", cor.lim)," = + ", sep = ""))) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = text_size),
          axis.text.y = element_text(size = text_size),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}






