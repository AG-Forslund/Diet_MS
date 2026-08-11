beta_div_plot <- function(data, limits = c(-0.12, 0.12), kind = "per_visit"){
  if(kind == "per_visit"){
    ggplot (data, aes (x = X1, y = X2, colour = as.factor(Visit), fill = as.factor(Visit))) +
    geom_point (aes (shape = as.factor(Visit)), size = 1) +
    scale_colour_manual(values = c("#882255", "#999933")) +
    geom_density_2d() +
    theme_classic() +
    scale_x_continuous(limits = limits) +
    scale_y_continuous(limits = limits) +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0))
  } else if (kind == "per_group") {
    ggplot (data, aes (x = X1, y = X2, colour = as.factor(Group), fill = as.factor(Group))) +
    geom_point (aes (shape = as.factor(Visit)), size = 1) +
    scale_colour_manual(values = group_palette) +
    geom_density_2d() +
    theme_classic() +
    scale_x_continuous(limits = limits) +
    scale_y_continuous(limits = limits) +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0))
  }
  
  
}


beta_div_plot_2 <- function(pcoa_data){
  plot <- pcoa_data %>%
    ggplot(aes(X1,X2, color = as.factor(Group), shape = Group)) +
    geom_point(size = 3) +
    geom_density_2d(contour_var = "density") +
    scale_colour_manual(values = group_palette) +
    scale_fill_manual(values = group_palette) +
    theme_classic() +
    scale_x_continuous(limits = c(min(pcoa_data$X1) - 0.1, max(pcoa_data$X1) + 0.1)) +
    scale_y_continuous(limits = c(min(pcoa_data$X2) - 0.1, max(pcoa_data$X2) + 0.1)) +
    labs(x = "PCOA 1", y = "PCOA2") +
    facet_wrap(~Visit)
  
  return(plot)
  
}


star_plot_Sofia <- function(data){
  dataV1G1 <- colMeans (data %>% filter(Group == "Control", Visit == "V0") %>% dplyr::select(c("X1", "X2")))
  dataV1G2 <- colMeans (data %>% filter(Group == "ketogenic", Visit == "V0") %>% dplyr::select(c("X1", "X2")))
  dataV1G3 <- colMeans (data %>% filter(Group == "fasting", Visit == "V0") %>% dplyr::select(c("X1", "X2")))
  dataV3G1 <- colMeans (data %>% filter(Group == "Control", Visit == "V3") %>% dplyr::select(c("X1", "X2")))
  dataV3G2 <- colMeans (data %>% filter(Group == "ketogenic", Visit == "V3") %>% dplyr::select(c("X1", "X2")))
  dataV3G3 <- colMeans (data %>% filter(Group == "fasting", Visit == "V3") %>% dplyr::select(c("X1", "X2")))
  
  means_df <- data.frame(
    Visit = rep(c("V0", "V3"), each = 3),
    Group = as.factor(rep(c("Control", "ketogenic", "fasting"), 2)),
    X1 = c(dataV1G1[1], dataV1G2[1], dataV1G3[1], dataV3G1[1], dataV3G2[1], dataV3G3[1]),
    X2 = c(dataV1G1[2], dataV1G2[2], dataV1G3[2], dataV3G1[2], dataV3G2[2], dataV3G3[2])
  )
  
  ggplot(data %>% mutate(Group = as.factor(Group)), aes(x = X1, y = X2)) +
    theme_classic() +
    geom_segment(aes(xend = 0, yend = 0, color = Group)) +
    #annotate("point", color = "#117733", x = dataV1G1 [1], y = dataV1G1 [2], size = 40, alpha = 0.3) +
    #annotate("point", color = "#661100", x = dataV1G2 [1], y = dataV1G2 [2], size = 40, alpha = 0.3) +
    #annotate("point", color = "#332288", x = dataG3 [1], y = dataV1G3 [2], size = 40, alpha = 0.3) +
    geom_point(data = means_df, aes(x = X1, y = X2, color = as.factor(Group)), size = 20, alpha = 0.3) +
    geom_label(aes(label = PatientID, color = Group), size = 2, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = group_palette) +
    scale_fill_manual(values = group_palette) +
    facet_wrap(~Visit)+
    coord_cartesian (xlim = c (min(data$X1), max(data$X2)), ylim = c (min(data$X2), max(data$X2))) +
    theme(legend.position = "bottom")
  
}


star_plot_diffs <- function(data){
  centered_data <- data %>%
    group_by(PatientID) %>%
    mutate(X1_ref = if(any(Visit == "V0")) X1[Visit == "V0"] else first(X1), X2_ref = if(any(Visit == "V0")) X2[Visit == "V0"] else first(X2)) %>%
    mutate(X1_centered = X1 - X1_ref, X2_centered = X2 - X2_ref)
  
  plot <- centered_data %>%
    ggplot(aes(x = X1_centered, y = X2_centered, color = as.factor(Group), group = PatientID, label = PatientID)) +
    theme_classic() +
    #geom_point() +
    geom_line(aes(group = PatientID)) +
    geom_label(data = centered_data %>% filter(Visit == "V3"), size = 2) +
    scale_colour_manual(values = group_palette) +
    theme(legend.position = "bottom")
  
  return(plot)
}


alpha_div_plot <- function(data, alpha.div.metric){
  ggplot(data, aes(Group,diversity_shannon, group = Group, col = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.05) +
    scale_colour_manual(values = group_palette) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0)) +
    ylab(paste(alpha.div.metric, "[%]", sep = " "))
}


MetadeconfoundR_Heatmap2 <- function(metaDeconfOutput, q_cutoff = 0.1, star_size = 7, x_axis_size = 7, y_axis_size = 7) {
  status <- data.frame(metaDeconfOutput$status) %>%
    #filter(if_any(everything(), ~ . != "NS")) %>%
    #mutate(across(everything(), .fns = function(x){na_if(x, "NS")})) %>%
    #dplyr::select(where(~!all(is.na(.x)))) %>%
    rownames_to_column("Metabolites") %>%
    gather("Variable", "Confound", -Metabolites)
  
  Ds <- data.frame(metaDeconfOutput$Ds) %>%
    dplyr::select(unique(status$Variable)) %>%
    rownames_to_column("Metabolites") %>%
    filter(Metabolites %in% status$Metabolites) %>%
    gather("Variable", "EffectSize", -Metabolites)
  
  Qs <- data.frame(metaDeconfOutput$Qs) %>%
    dplyr::select(unique(status$Variable)) %>%
    rownames_to_column("Metabolites") %>%
    filter(Metabolites %in% status$Metabolites) %>%
    gather("Variable", "qs", -Metabolites)
  
  merged.frames <- merge(status, Ds, by = c("Metabolites", "Variable"))
  plotframe <- merge(merged.frames, Qs, by = c("Metabolites", "Variable")) %>%
    #filter(qs <= q_cutoff) %>%
    mutate(Stars = ifelse((Confound %in% c("OK_sd", "OK_d", "OK_nc", "AD")) & qs <= q_cutoff/100, "***", 
                          ifelse((Confound %in% c("OK_sd", "OK_d", "OK_nc", "AD")) & qs <= q_cutoff/10 & qs > q_cutoff/100, "**",
                                 ifelse((Confound %in% c("OK_sd", "OK_d", "OK_nc", "AD")) & qs <= q_cutoff & qs > q_cutoff/10, "*",
                                        ifelse(!is.na(Confound) & qs <= q_cutoff/10 & qs > q_cutoff/100, "+++",
                                               ifelse(!is.na(Confound) & qs <= q_cutoff/100, "++",
                                                      ifelse(!is.na(Confound) & qs <= q_cutoff & qs > q_cutoff/10, "+","")))))))
  
  p <- plotframe %>%
    ggplot(aes(Variable, Metabolites, fill = EffectSize)) +
    geom_tile() +
    geom_text(aes(label = Stars), size = star_size) +
    scale_fill_distiller(palette = "PRGn", name = "Correlation\nEffect Size", limits = c(-1,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = x_axis_size, angle = 45, hjust = 1),
          axis.text.y = element_text(size = y_axis_size, angle = 0, hjust = 1, vjust = 0)) +
    labs(title="MetaDeconfoundR summarizing heatmap",
         subtitle=paste("FDR-corrected p-values:", paste(paste(paste(q_cutoff/100," = ***"), paste(q_cutoff/10," = **"), sep = ", "), paste(q_cutoff," = *"), sep = ", ")), x = "Metadata variables", y = "Metabolites")
  
  return(p)
}


make_metavariables_plot <- function(variables = variables, visit_var = visit_var, matrix_var = matrix_var, color_by = color_var) {
  
  
  if(matrix_var == "stool" & visit_var == "percent") {
    matrix_df <- stool_metabolomics_metadata_diffs
  } else if (matrix_var == "serum" & visit_var == "percent") {
    matrix_df <- serum_metabolomics_metadata_diffs
  } else if(matrix_var == "serum" & visit_var != "percent") {
    matrix_df <- serum_metabolomics_metadata %>%
      filter(Visit == visit_var)
  } else if (matrix_var == "stool" & visit_var != "percent") {
    matrix_df <- stool_metabolomics_metadata %>%
      filter(Visit == visit_var)
  }
  
  y.class <- check_col_kind(matrix_df[,variables[1]])
  z.class <- check_col_kind(matrix_df[,variables[2]])
  col.class <- check_col_kind(matrix_df[,color_by])
  
  if(y.class == "Factor" | y.class == "Binary" ) {
    cor_res <- format(kruskal.test(matrix_df[,variables[2]], matrix_df[,variables[1]])$p.value,scientific = T)
  } else if ((z.class == "Factor" | z.class == "Binary" )){
    cor_res <- format(kruskal.test(matrix_df[,variables[1]], matrix_df[,variables[2]])$p.value,scientific = T)
  } else {
    cor_res <- format(cor.test(matrix_df[,variables[1]], matrix_df[,variables[2]], method = "spearman")$p.value, scientific = T)
  }
  
  which.scale.colour <- function(plot){
    if(col.class == "Ranked" | col.class == "Normally.Distributed") {
      return(scale_color_viridis_c(option = "plasma"))
    } else {
      if(color_by == "Group") {
        color_pal <- c("#117733", "#332288", "#661100")
      } else {
        color_pal <- safe_colorblind_palette[1:levelsize]
      }
      return(scale_color_manual(values = color_pal))
    }
  }
  
  if(col.class == "Binary") {
    matrix_df <- matrix_df %>%
      mutate(across(all_of(color_by), as.factor))
  }
  
  if (y.class == "Binary" & z.class == "Binary") {
    p <- matrix_df %>%
      #dplyr::select(variables) %>%
      mutate(Sum = as.factor(as.numeric(matrix_df[,variables[1]]) + as.numeric(matrix_df[,variables[2]]))) %>%
      mutate(Match = ifelse(Sum == 1, FALSE, TRUE)) %>%
      ggplot(aes(x = Match, fill = Sum)) +
      geom_bar(position = "stack") +
      ylab("Absolute Count") +
      annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
      theme_linedraw()
    
  } else if (y.class == "Factor" & z.class == "Factor") {
    
    p <- matrix_df %>%
      mutate(F1 = as.factor(matrix_df[,variables[1]])) %>%
      mutate(F2 = as.factor(matrix_df[,variables[2]])) %>%
      ggplot(aes(F1, F2)) +
      geom_jitter(width = 0.2, height = 0.2, aes_string(color = color_by)) +
      which.scale.colour(.) +
      xlab(variables[1]) +
      ylab(variables[2]) +
      annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
      theme_linedraw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else if(y.class == "Normally.Distributed" & z.class == "Normally.Distributed") {
    
    p <- ggplot(matrix_df, aes_string(variables[1], variables[2])) +
      geom_point(aes_string(color = color_by)) +
      annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
      theme_linedraw()
    
  } else if(y.class == "Ranked" & z.class == "Ranked") {
    
    p <- ggplot(matrix_df, aes_string(variables[1], variables[2])) +
      geom_point(aes_string(color = color_by)) +
      which.scale.colour(.) +
      annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
      theme_linedraw()
    
  } else if(y.class != z.class) {
    
    if ((y.class == "Normally.Distributed" | y.class == "Ranked") & (z.class == "Ranked" | z.class == "Normally.Distributed")){
      p <- ggplot(matrix_df, aes_string(variables[1], variables[2])) +
        geom_point(aes_string(color = color_by)) +
        which.scale.colour(.) +
        annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
        theme_linedraw()
      
    } else if ((y.class == "Normally.Distributed" | y.class == "Ranked" | y.class == "Binary") & (z.class == "Binary" | z.class == "Ranked" | z.class == "Normally.Distributed")) {
      
      a <- ifelse(y.class == "Binary", 1, 2)
      p <- matrix_df %>%
        mutate(Bin = as.factor(matrix_df[,variables[a]])) %>%
        ggplot(aes_string(y = variables[-a], color = color_by)) +
        geom_boxplot(aes(x = Bin), position = position_dodge(0.8), outlier.shape = NA) +
        geom_point(aes(x = Bin), position = position_dodge(0.8), size = 1) +
        which.scale.colour(.) +
        annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
        xlab(variables[a]) +
        theme_linedraw()
      
    } else if (y.class == "Factor" | z.class == "Factor") {
      
      if (z.class == "Binary" | y.class == "Binary") {
        
        a <- ifelse(y.class == "Binary", 1, 2)
        p <- matrix_df %>%
          mutate(Bin = as.factor(matrix_df[,variables[a]])) %>%
          ggplot(aes_string(fill = variables[-a])) +
          geom_bar(aes(x = Bin), position = position_fill()) +
          which.scale.colour(.) +
          annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
          xlab(variables[a]) +
          ylab("Count in %") +
          theme_linedraw()
        
      } else if (z.class == "Normally.Distributed" | z.class == "Ranked" | y.class == "Normally.Distributed" | y.class == "Ranked") {
        
        a <- ifelse(y.class == "Factor", 1, 2)
        p <- matrix_df %>%
          ggplot(aes_string(x = variables[a], y = variables[-a])) +
          geom_boxplot() +
          geom_point(aes_string(color = color_by)) +
          annotate("text", x = Inf, y = Inf, label = paste("p-Value =", cor_res, sep = ""), hjust = 1.2, vjust = 1.2) +
          which.scale.colour(.) +
          xlab(variables[a]) +
          theme_linedraw()
        
      }
    } 
  }
  return(p)
}


interaction_heatmap <- function(results, filter = T, limits = c(0,0.3), q.limit = 0.1) {
  
  if (filter == T) {
    results_dataQ <- results$Qs %>%
      filter(if_any(!Feature, ~. < q.limit)) %>%
      dplyr::select(Feature, where(~ any(. < q.limit)))
  } else {
    results_dataQ <- results$Qs
  }
  
  results_dataQ <- results_dataQ %>%
    gather(Taxon, Q, -Feature) %>%
    mutate(Star = if_else(Q < q.limit/10, "***", if_else(Q < q.limit/2, "**", if_else(Q < q.limit, "*", NA))))
  
  results_dataD <- results$Ds %>%
    filter(Feature %in% results_dataQ$Feature) %>%
    dplyr::select(Feature, all_of(results_dataQ$Taxon)) %>%
    gather(Taxon, D, -Feature)
  
  results_data <- merge(results_dataQ, results_dataD, by = c("Feature", "Taxon"))
  
  interaction_pPlot <- ggplot(results_data, aes(Feature, Taxon, fill = D, label = Star)) +
    geom_tile() +
    geom_text() +
    scale_fill_viridis_c(option = "plasma", name = "Partial R2", begin = 0.15, end = 1, limits = limits) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0))
  
  return(interaction_pPlot)
}


Interaction_Heatmap2 <- function(result, p_limit = 0.05, q_cutoff = 0.1, star_size = 7, x_axis_size = 7, y_axis_size = 7, cor_limits = c(0, 0.1), p_adj_method = NULL) {
  relevant_covariates <- result$Qs %>%
    filter(Group <= p_limit) %>%
    dplyr::select(Group, where(~ any(. <= p_limit, na.rm = T))) %>%
    names()

  relevant_features <- result$Qs %>%
    filter(Group <= p_limit) %>%
    pull("Feature")
  
  result_filtered <- lapply(result, function(df) {
    df <- df %>%
      filter (Feature %in% relevant_features) %>%
      dplyr::select(Feature, all_of(relevant_covariates)) %>%
      gather(Covariate,Value,-Feature)
  })

  
  merged.frames <- merge(merge(result_filtered$Qs, result_filtered$Ds, by = c("Feature", "Covariate"), suffixes = c("_Significance", "_EffSize", "_Status")), result_filtered$Status, by = c("Feature", "Covariate"))
  
  confounded_covariates <- merged.frames %>%
    filter(Value == "Group") %>%
    pull(Covariate)
  
  plotframe <- merged.frames %>%
    filter(!Covariate %in% confounded_covariates) %>%
    mutate(Stars = ifelse((Value %in% c("Both", "None", "Deconfounded")) & Value_Significance <= q_cutoff/100, "***", 
                          ifelse((Value %in% c("Both", "None", "Deconfounded")) & Value_Significance <= q_cutoff/10 & Value_Significance > q_cutoff/100, "**",
                                 ifelse((Value %in% c("Both", "None", "Deconfounded")) & Value_Significance <= q_cutoff & Value_Significance > q_cutoff/10, "*",
                                        ifelse((Value %in% c("Confounding")) & Value_Significance <= q_cutoff/100, "°°°", 
                                               ifelse((Value %in% c("Confounding")) & Value_Significance <= q_cutoff/10 & Value_Significance > q_cutoff/100, "°°",
                                                      ifelse((Value %in% c("Confounding")) & Value_Significance <= q_cutoff & Value_Significance > q_cutoff/10, "°",
                                                          ifelse(!is.na(Value) & Value_Significance <= q_cutoff/10 & Value_Significance > q_cutoff/100, "+++",
                                                                 ifelse(!is.na(Value) & Value_Significance <= q_cutoff/100, "++",
                                                                        ifelse(!is.na(Value) & Value_Significance <= q_cutoff & Value_Significance > q_cutoff/10, "+",""))))))))))
  if(p_adj_method == "none") {
    header <- "Raw p-values:"
  } else {
    header <- "FDR-corrected q-values:"
  }
  
  subtitle <- paste(plotframe %>% filter(Covariate == "Group") %>% pull(Value), collapse = ", ")
   
  p <- plotframe %>%
    ggplot(aes(Feature, Covariate, fill = abs(Value_EffSize))) +
    geom_tile() +
    geom_text(aes(label = Stars), size = star_size) +
    scale_fill_viridis_c(option = "plasma", name = "Partial R2", begin = 0.15, end = 1, limits = cor_limits) +
    #scale_fill_distiller(palette = "PRGn", name = "Partial\nR2", limits = cor_limits) +
    theme_classic() +
    theme(axis.text.x = element_text(size = x_axis_size, angle = 45, hjust = 1),
          axis.text.y = element_text(size = y_axis_size, angle = 0, hjust = 1, vjust = 0)) +
    labs(title="Interaction calculation summarizing heatmap",
         subtitle=paste(header, paste(paste(paste(q_cutoff/100," = ***"), paste(q_cutoff/10," = **"), sep = ", "), paste(q_cutoff," = *"), sep = ", ")), x = "Metadata variables", y = "Features")
  
  return(p)
}


interaction_heatmap3 <- function(results, filter = T, limits = c(-0.3,0.3), triangle_size = 3, star_size = 5) {
  
  if (filter == T) {
    results_dataQ <- results$Qs %>%
      filter(if_any(!Feature, ~. < 0.1)) %>%
      dplyr::select(Feature, where(~ any(. < 0.1)))
  } else {
    results_dataQ <- results$Qs
  }
  
  results_dataQ <- results_dataQ %>%
    gather(Taxon, Q, -Feature) %>%
    mutate(Star = if_else(Q < 0.01, "***", if_else(Q < 0.05, "**", if_else(Q < 0.1, "*", NA))))
  
  results_dataD <- results$Ds %>%
    filter(Feature %in% results_dataQ$Feature) %>%
    dplyr::select(Feature, all_of(results_dataQ$Taxon)) %>%
    gather(Taxon, D, -Feature)
  
  results_data <- merge(results_dataQ, results_dataD, by = c("Feature", "Taxon"))
  
  interaction_pPlot <- ggplot(results_data, aes(Feature, Taxon)) +
    geom_point(aes(color = ifelse(Q > 0.1, "gray70", "black"), shape = ifelse(D > 0, 24, 25), fill = D), size = triangle_size) +
    geom_text(aes(label = Star), size = star_size) +
    #scale_color_viridis_c(option = "plasma", name = "Partial R2", begin = 0.15, end = 1, limits = limits) +
    #scale_fill_viridis_c(option = "plasma", name = "Partial R2", begin = 0.15, end = 1, limits = limits) +
    scale_fill_distiller(palette = "PRGn", name = "Partial\nR2", limits = limits) +
    scale_color_identity() +
    scale_shape_identity() +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0))
  
  return(interaction_pPlot)
}


Interaction_posthoc_plot <- function(result_list, type = NULL) {
  if(type == "metabolomics") {
    if (Trp_ratio_calc == T) {
      features_data <- merge(Trp_ratios %>%
        mutate(Sample.Identification = paste("S", Sample.Identification, sep = "")) %>%
        column_to_rownames("Sample.Identification") %>%
        dplyr::select(-starts_with("Batch"), -starts_with("Sample.Type")) %>%
          rownames_to_column("Tube_ID"),
        SampleIDs, by = "Tube_ID")
    } else {
      features_data <- merge(metabolomics_data %>%
                               mutate(Sample.Identification = paste("S", Sample.Identification, sep = "")) %>%
                               column_to_rownames("Sample.Identification") %>%
                               dplyr::select(-starts_with("Batch"), -starts_with("Sample.Type")) %>%
                               rownames_to_column("Tube_ID"),
                             SampleIDs, by = "Tube_ID")
    }
  } else if (type == "microbiome") {
    features_data <- merge(data.frame(t(microbiome_data)) %>% rownames_to_column("Tube_ID"), SampleIDs, by = "Tube_ID")
  }
  
  sig.features <- result_list$Status %>%
    filter(!Group == "NS") %>%
    pull(Feature)
  
  results_df <- bind_rows(lapply(sig.features, function(f) feature_posthoc_tests(feature = f, countdata = features_data, type = type))) #%>%
    #pivot_wider(names_from = Feature, values_from = P_Value)
  
  plot <- results_df %>%
    #gather(Feature, q_Value, -Comparison) %>%
    filter(P_Value < 0.1) %>%
    mutate(Comparison = as.factor(Comparison)) %>%
    mutate(Estimate_direction = factor(ifelse(Estimate > 0, "positive", "negative"))) %>%
    #ggplot(aes(x = Feature, y = Comparison, fill = P_Value, size = -log10(P_Value), shape = Estimate_direction)) +
    ggplot(aes(x = Feature, y = Comparison, fill = P_Value, size = abs(Estimate), shape = Estimate_direction)) +
    #geom_point(shape = 21, alpha = 0.7 ) +
    geom_point() +
    scale_shape_manual(values = c("positive" = 24, "negative" = 25)) +
    #scale_size_continuous(range = c(5, 15), name = "P-value") +
    scale_size_continuous(range = c(5, 15), name = "abs(Estimate)") +
    scale_fill_viridis(option = "mako", begin = 1, end = 0, name = "q-Value") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0)) 
    
  list <- list(result_df = results_df, result_plot = plot)
  
  return(list)
}



Interaction_posthoc_plot_outcome <- function(result, features_data, covariates_data, type = NULL, outcome_name = NULL) {
  
  sig.features <- result$Qs %>%
    dplyr::select(Feature, outcome_name) %>%
    filter(if_any(outcome_name, ~ . < 1)) %>%
    pull(Feature)

  combined_df <- merge(features_data %>% dplyr::select(all_of(sig.features)), covariates_data %>% dplyr::select(all_of(outcome_name), Group, Status, PatientID), by = 0) %>%
    column_to_rownames("Row.names") %>%
    rename(Visit = Status)
  
  results_df <- bind_rows(lapply(sig.features, function(f) feature_posthoc_tests(f, combined_df, outcome = T, outcome_name = outcome_name, type = type))) %>%
    pivot_wider(names_from = Feature, values_from = P_Value)
  
  plot <- results_df %>%
    gather(Feature, q_Value, -Comparison) %>%
    filter(q_Value < 0.1) %>%
    mutate(Comparison = as.factor(Comparison)) %>%
    ggplot(aes(Feature, Comparison, fill = q_Value, size = -log10(q_Value))) +
    geom_point(shape = 21, alpha = 0.7 ) +
    scale_size_continuous(range = c(1, 15), name = "P-value") +
    scale_fill_viridis(option = "mako", begin = 1, end = 0, name = "q-Value") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0)) 
  
  list <- list(result_df = results_df, result_plot = plot)
  
  return(list)
}



future_assoc_bars <- function(result_future, filter = T, q.limit = 0.1, bar.plot = T, y.limits = c(-0.2, 0.2)) {
  
  limits = c(0,q.limit)
  
  group_names <- names(result_future)[3:4]
  
  if (filter == T) {
    results_dataQ <- result_future$Qs %>%
      filter(if_any(!Feature, ~. < q.limit)) %>%
      dplyr::select(Feature, where(~ any(. < q.limit))) %>%
      gather(Outcome, Q, -Feature)

  } else {
    results_dataQ <- result_future$Qs %>%
      gather(Outcome, Q, -Feature)
  }
  
  results_dataD1 <- result_future[[group_names[1]]] %>%
    filter(Feature %in% results_dataQ$Feature) %>%
    dplyr::select(Feature, all_of(results_dataQ$Outcome)) %>%
    gather(Outcome, D, -Feature) %>%
    mutate(Comparison = group_names[1])
  
  results_dataD2 <- result_future[[group_names[2]]] %>%
    filter(Feature %in% results_dataQ$Feature) %>%
    dplyr::select(Feature, all_of(results_dataQ$Outcome)) %>%
    gather(Outcome, D, -Feature) %>%
    mutate(Comparison = group_names[2])
  
  results_data <- merge(results_dataQ, rbind(results_dataD1, results_dataD2), by = c("Feature", "Outcome")) %>%
    mutate(Outcome2 = Outcome) %>%
    unite(Combination, c("Feature", "Outcome"))
  
  if (filter == T) {
    results_data <- results_data %>%
      filter(Q < q.limit)
  }
  
  if (bar.plot) {
  future_Plot <- ggplot(results_data, aes(x = Combination, y = D, fill = Q, group = Comparison)) +
    geom_col(position = position_dodge(width = 0.9)) +
    scale_fill_viridis(option = "rocket", name = "Q-Value", limits = limits, direction = -1, begin = 0.2, end = 0.9) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    ylim(y.limits) +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0))+
    facet_wrap(~Outcome2)
  } else {
    future_Plot <- ggplot(results_data, aes(y = interaction(Combination, Comparison, sep = "_"))) +
      geom_segment(aes(x = 0, xend = D, color = Q), linewidth = 1, arrow = arrow(length = unit(2, "cm"), type = "closed", angle = 25)) +
      scale_color_viridis(option = "cividis", name = "Q-Value", limits = limits, direction = -1, begin = 0, end = 1) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      ylim(y.limits) +
      labs(x = "Effect Size (D)", y = "Feature - Group", title = "Effect Size Directionality by Feature Group") +
      theme_minimal()
  }
  
  return(future_Plot)
}


tree_func <- function(final_model, tree_num) {
  
  # get tree by index
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  
  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
  
  # plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  
  print(plot)
}





