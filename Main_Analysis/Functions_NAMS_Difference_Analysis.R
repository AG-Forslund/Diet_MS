feature_association <- function(feature_data, covariate_data, covariate_name, type = type, analysis_type = analysis_type) {
  # Merge feature data with covariates
  data_general <- data.frame(feature_data, covariate_data, check.names = FALSE) %>%
    rename_with(~ "Feature", .cols = contains("_Diff")) %>%
    rename_with(~ "Base", .cols = contains("_Base"))
  
  if(type == "metabolomics" & analysis_type == "biological") {
    data_general <- data_general %>%
      mutate(Base = log(Base))
  } else if(type == "microbiome" & analysis_type == "biological") {
     data_general <- data_general %>%
      mutate(Base = rank(Base)) %>%
      mutate(Feature = rank(Feature))
  } else if (analysis_type == "outcomes") {
    data_general <- data_general %>%
      mutate(Base = rank(Base))
  }
  
  # Is the group effect on percentage change depended on baseline values?
  general_model <- glm(Feature ~ data_general[[covariate_name]]*Base, data = data_general)
  general_model_no_interaction <- glm(Feature ~ Base + data_general[[covariate_name]], data = data_general)
  general_test <- anova(general_model, general_model_no_interaction)

  
  #Calculate effect size
  r2_full <- r.squaredGLMM(general_model)[1]
  r2_reduced <- r.squaredGLMM(general_model_no_interaction)[1]
  f_squared <- (r2_full - r2_reduced) / (1 - r2_full)
  
  if(f_squared < 0) {
    f_squared = 0
  }
  
  if(covariate_name == "Group") {
    data_general$Group <- relevel(as.factor(data_general$Group), ref = "ketogenic diet")
    posthoc_model <- glm(Feature ~ Group*Base, data = data_general)
    slopes <- emtrends(posthoc_model, specs = "Group", var = "Base")
    pairwise_contrasts <- pairs(slopes)
    print(summary(pairwise_contrasts, infer = c(TRUE, TRUE)))
  }
  
  
  #Summarize results
  result_list <- list(pValue = general_test$"Pr(>F)"[2], EffSize = f_squared)
  
  return(result_list)
}


analyze_feature <- function(feature_data, covariate_name, covariate_data, p_limit = 0.1, type = type, analysis_type = analysis_type) {
  #create dataframe specific for covariate (without its NA rows for example)
  data <- data.frame(feature_data, covariate_data, check.names = FALSE) %>%
    rename_with(~ "Feature", .cols = contains("_Diff")) %>%
    rename_with(~ "Base", .cols = contains("_Base"))
  
  if(type == "metabolomics" & analysis_type == "biological") {
    data <- data %>%
      mutate(Base = log(Base)) %>%
      drop_na(all_of(covariate_name))
  } else if(type == "microbiome" & analysis_type == "biological") {
    data <- data %>%
      mutate(Base = rank(Base)) %>%
      mutate(Feature = rank(Feature)) %>%
      drop_na(all_of(covariate_name))
  } else if (analysis_type == "outcomes") {
    data <- data %>%
      mutate(Base = rank(Base)) %>%
      drop_na(all_of(covariate_name))
  }
  
  if(length(unique(data[[covariate_name]])) <=1) {
    return("NS")
  } else {
    # Fit the model with that potentially trimmed data with and without the covariate term as interaction
    full_model <- glm(Feature ~ Base + Group + Base:Group + data[[covariate_name]], data = data)
    no_group_interaction_model <- glm(Feature ~ Base + Group + data[[covariate_name]], data = data)
    no_covariate_interaction_model <- glm(Feature ~ Base + Group + Base:Group, data = data)
      
    # Compare models
    interaction_test <- anova(full_model, no_group_interaction_model)
    covariate_test <- anova(full_model, no_covariate_interaction_model)
      
    # Return results
    GroupTest = interaction_test$"Pr(>F)"[2]
    CovariateTest = covariate_test$"Pr(>F)"[2]
      
    if (GroupTest < p_limit & CovariateTest < p_limit){
        return("Both")
    } else if (GroupTest > p_limit & CovariateTest > p_limit){
        return("None")
    } else if (GroupTest < p_limit & CovariateTest > p_limit) {
        return("Group_Base_Interaction")
    } else if(GroupTest > p_limit & CovariateTest < p_limit) {
        return("Confounding")
    }
  } 
}


interaction_deconf <- function(features_data, covariates_data, p_limit, p_adj_method = "fdr", features, type = type, analysis_type = analysis_type) {
  
  patterns <- c("Both", "None", "Group", "Confounding")
  covariates <- setdiff(names(covariates_data), c("PatientID", "Group"))
  results_df <- data.frame(matrix(ncol = length(covariates) + 2, nrow = length(features)))
  colnames(results_df) <- c("Feature", "Group", covariates)
  results_df$Feature <- features
  
  results <- list(Ps = results_df, Qs = results_df, Ds = results_df, Status = results_df)
  
  for (feature in features) {
    #print(paste("Feature:", feature))
    feature_diff <- paste(feature, "_Diff", sep = "")
    feature_base <- paste(feature, "_Base", sep = "")
    feature_data <- features_data[, c(feature_diff, feature_base)]
    
    for (covariate in c("Group", covariates)) {
      #print(paste("Covariate:", covariate))
      if (str_detect(covariate, feature, negate = T)) {
        association_res <- feature_association(feature_data, covariates_data, covariate, type = type, analysis_type = analysis_type)
        results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
        results$Ds[[covariate]][results$Ds$Feature == feature] <- association_res$EffSize
      }
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  
  #Calculate if it's the covariate or the Interaction term that is significant for all features which are associated with the interaction term
  for(feature in features) {
    #print(paste("Feature:", feature))
    feature_diff <- paste(feature, "_Diff", sep = "")
    feature_base <- paste(feature, "_Base", sep = "")
    feature_data <- features_data[, c(feature_diff, feature_base)]
    
    if (results$Qs$Group[results$Qs$Feature == feature] < p_limit) {
      for (covariate in covariates) {
        #print(paste("Covariate:", covariate))
        if (str_detect(covariate, feature, negate = T)) {
          if (results$Qs[[covariate]][results$Qs$Feature == feature] < p_limit) {
            results$Status[[covariate]][results$Status$Feature == feature] <- analyze_feature(feature_data, covariate, covariates_data, p_limit, type = type, analysis_type = analysis_type)
          } else {
            results$Status[[covariate]][results$Status$Feature == feature] <- "NS"
          }
          }
      }
    
    
      Group_conf <- results$Status %>%
        filter(Feature == feature) %>%
        rowwise() %>%
        mutate(has_confounded = any(c_across(-Feature) == "Confounding", na.rm = TRUE)) %>%
        ungroup() %>%
        pull(has_confounded) %>%
        any(na.rm = TRUE)
      
      Conf_variables <- results$Status %>%
        filter(Feature == feature) %>%
        summarise(across(-Feature, ~. == "Confounding")) %>%
        dplyr::select(where(~any(. == TRUE, na.rm = T))) %>%
        names()
      
      if (results$Qs$Group[results$Qs$Feature == feature] > p_limit) {
        results$Status$Group[results$Status$Feature == feature] <- "NS"
      } else if(Group_conf == T) {
        results$Status$Group[results$Status$Feature == feature] <- paste(Conf_variables, collapse = ", ")
      } else if (Group_conf == F) {
        results$Status$Group[results$Status$Feature == feature] <- "Deconfounded"
      }
    }
  }
  
  return(results)
}
