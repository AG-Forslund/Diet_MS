check_col_kind <- function(x){
  if (class(x) == "factor"){
    col.class <- "Factor"
  } else if(all(na.omit(x) %in% 0:1) == TRUE){
    col.class <- "Binary"
  } else if(class(x) == "integer" | class(x) == "numeric"){
    if(shapiro.test(x)$p.value >= 0.05){
      #if(which(check_distribution(na.omit(x))$p_Vector == max(check_distribution(na.omit(x))$p_Vector)) == 12) {
      col.class <- "Normally.Distributed"
    } else {
      col.class <- "Ranked"
    }
  }
  return(col.class)
}



calc.metabolite.correlations <- function(pdata, mdata, p.adj.method = "BH"){
  
  res.data.frame <- data.frame(matrix(nrow = length(metabolites), ncol = length(mdata.names)))
  res.data.frame_cor = res.data.frame
  
  for (i in 1:length(metabolites)) {
    metabolite <- metabolites[i]
    metabolite.areas <- subset(pdata, pdata[comp.name] == metabolite) %>%
      select(all_of(c(comp.name,metric,meta.merge)))
    
    metabolite.data <- merge(mdata, metabolite.areas, by = meta.merge) %>%
      drop_na(all_of(metric))
    
    metabolite.area_final <- metabolite.data[,metric]
    #print(i)
    for (j in 1:length(mdata.names)){
      mdata.name <- metabolite.data[,mdata.names[j]]
      
      if(length(na.omit(mdata.name)) <= 4){
        p <- NA
        cor <- NA
      } else {
        mdata.name_class <- check_col_kind(mdata.name)
        
        if (mdata.name_class == "Factor") {
          levelsize <- length(unique(mdata.name))
        }
        
        #print(j)
        
        if (mdata.name_class == "Binary"){
          p <- cor.test(mdata.name, metabolite.area_final, method = "pearson")$p.value
          cor <- cor.test(mdata.name, metabolite.area_final, method = "pearson")$estimate
        } else if (mdata.name_class == "Factor") {
          if(length(levels(mdata.name)) > length(unique(mdata.name)) | length(levels(mdata.name)) == 1){
            p <- NA
            cor <- NA
          } else {
            if (levelsize == 2) {
              p <- wilcox.test(metabolite.area_final ~ mdata.name, exact = F)$p.value # Man Whitney U test
              cor <- wilcoxonR(metabolite.area_final, mdata.name)
            } else {
              p <- kruskal.test(metabolite.area_final ~ mdata.name)$p.value # Kruskal Wallis test
              cor <- epsilonSquared(metabolite.area_final, mdata.name)
            }
            #one.way.anova <- aov(metabolite.area_final ~ mdata.name)
            #p <- summary(one.way.anova)[[1]][["Pr(>F)"]][1]
            #cor <- eta_squared(one.way.anova)
          }
        } else if (mdata.name_class == "Normally.Distributed" | mdata.name_class == "Ranked") {
          p <- cor.test(mdata.name, metabolite.area_final, method = "spearman")$p.value
          cor <- cor.test(mdata.name, metabolite.area_final, method = "spearman")$estimate
        }
      }
      res.data.frame[i,j] <- p
      res.data.frame_cor[i,j] <- cor
    }
  }
  res.data.frame_FDR.adjusted <- data.frame(lapply(res.data.frame, function(x) p.adjust(x, method = p.adj.method)))
  result.list <- list(p.Value = res.data.frame, adjusted = res.data.frame_FDR.adjusted, cor.estimate = res.data.frame_cor) #adjust over metabolites
  result.list <- lapply(result.list, function(x){
    rownames(x) <- metabolites
    names(x) <- mdata.names
    return(x)
  })
  return(result.list)
}



calc.nas.per.group <- function(merged_file, groups){
  result_list <- vector("list", length = length(groups))
  for (i in 1:length(groups)){
    group <- groups[i]
    nas <- merged_file %>%
      mutate(across(!!group, ~ as.factor(.))) %>%
      rename(NAGroup = !!group) %>%
      group_by(NAGroup) %>%
      summarise_each(~ 100*mean(is.na(.)))
      #summarise_all(~ (sum(is.na(.)/n())*100))

    result_list[i] <- list(nas)
  }
  
  names(result_list) <- groups
  
  return(result_list)
}



calc.RSDs <- function(file, QC_regex = "Pooled QC 1"){
  
    result_list <- vector("list", length = length(unique(file$Batch)))
    
    # For QCs
    for(i in 1:length(unique(file$Batch))) {
      QC_Batch_n <- paste("Batch", i, sep = "")
      
      rsds_QC <- file %>%
        filter(grepl(QC_regex, Sample.Type)) %>%
        filter(Batch == QC_Batch_n) %>%
        select(-c(Sequence.Position, Batch, Sample.Type)) %>%
        column_to_rownames("Sample.Identification") %>%
        summarize(across(everything(),~ sd(.x, na.rm = TRUE)/mean(.x, na.rm = TRUE)*100))
      
      result_list[i] <- list(rsds_QC)
      names(result_list)[i] <- QC_Batch_n
    }
    
    m_QCs <- bind_rows(result_list)
    rownames(m_QCs) <- paste(names(result_list), "QC", sep=".")
    
    # Now for samples
    for(i in 1:length(unique(file$Batch))) {
      QC_Batch_n <- paste("Batch", i, sep = "")
      
      rsds_sample <- file %>%
        filter(grepl("Sample", Sample.Type)) %>%
        filter(Batch == QC_Batch_n) %>%
        select(-c(Sequence.Position, Batch, Sample.Type)) %>%
        column_to_rownames("Sample.Identification") %>%
        summarize(across(everything(),~ sd(.x, na.rm = TRUE)/mean(.x, na.rm = TRUE)*100))
      
      result_list[i] <- list(rsds_sample)
      names(result_list)[i] <- QC_Batch_n
    }

    m_samples <- bind_rows(result_list)
    rownames(m_samples) <- paste(names(result_list), "Samples", sep=".")
    
    RSD_result <- rbind(m_samples, m_QCs) %>%
      rownames_to_column("Kind") %>%
      mutate(Batch = str_extract(Kind, "Batch\\d+")) %>%
      mutate(Kind = str_remove(Kind, "Batch\\d+."))
  
  return(RSD_result)
}



batch_diff_sig <- function(file, compounds){
  table <- tibble(compounds) %>%
    mutate(MCT = sapply(compounds, function(c) {
      #summary(aov(file[,c] ~ file[,"Batch"]))[[1]][["Pr(>F)"]][1]
      kruskal.test(file[,c] ~ file[,"Batch"])$p.value})) %>%
    mutate(BH = p.adjust(MCT)) %>%
    mutate(B1B2 = unlist(sapply(compounds, function(c) {
      dunn_test(file, as.formula(paste(c, "Batch", sep = "~")))[1,"p.adj"]}))) %>%
    mutate(B1B3 = unlist(sapply(compounds, function(c) {
      dunn_test(file, as.formula(paste(c, "Batch", sep = "~")))[2,"p.adj"]}))) %>%
    mutate(B2B3 = unlist(sapply(compounds, function(c) {
      dunn_test(file, as.formula(paste(c, "Batch", sep = "~")))[3,"p.adj"]})))
  
  return(table)
}