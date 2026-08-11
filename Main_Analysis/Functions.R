feature_association <- function(feature_data, covariate_data, covariate_name, type = NULL, visit = F) {
  # Merge feature data with covariates
  if(type == "metabolomics") {
      data_general <- data.frame(Feature = log(feature_data), covariate_data, check.names = FALSE) %>% mutate(Status = as.numeric(Status))
  } else if(type == "microbiome") {
      data_general <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% mutate(Status = as.numeric(Status))
  } else if (type == "outcome") {
      #data_general <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% mutate("{covariate_name}" := rank(!!sym(covariate_name)))
      data_general <- data.frame(Feature = feature_data, covariate_data, check.names = FALSE) %>% mutate(Status = as.numeric(Status))
  }
  
  #Test if interaction is generally associated with this feature regardless of covariate
  if (visit) {
    general_model <- lmer(Feature ~ Status + (1 | PatientID) , data = data_general, REML = F)
    general_model_no_interaction <- lmer(Feature ~ 1 + (1 | PatientID), data = data_general, REML = F)
  } else {
    general_model <- lmer(Feature ~ Status + data_general[[covariate_name]] + Status:data_general[[covariate_name]] + (1 | PatientID), data = data_general, REML = F)
    general_model_no_interaction <- lmer(Feature ~ Status + data_general[[covariate_name]] + (1 | PatientID), data = data_general, REML = F)
  }
  
  general_test <- anova(general_model, general_model_no_interaction)
  
  #Calculate effect size
  r2_full <- r.squaredGLMM(general_model)[2]
  r2_reduced <- r.squaredGLMM(general_model_no_interaction)[2]
  f_squared <- (r2_full - r2_reduced) / (1 - r2_full)
  
  if(f_squared < 0) {
    f_squared = 0
  }
  #Summarize results
  result_list <- list(pValue = general_test$"Pr(>Chisq)"[2], EffSize = f_squared)
  
  return(result_list)
}


### This function is called within the interaction_deconf_both() function
feature_association_both_multivariate <- function(feature_data, covariate_data, covariate_name, type = NULL) {
  # Merge feature data with covariates
  if(type == "microbiome") {
    data_general <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  } else if(type == "metabolomics") {
    data_general <- data.frame(Feature = log(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  }
  
  #Test if the association between the microbial feature and the metabolite is dependent on the group
  general_model <- lmer(Feature ~ Status + data_general[[covariate_name]] + Group + data_general[[covariate_name]] : Group + (1 | PatientID), data = data_general, REML = F)
  general_model_no_interaction <- lmer(Feature ~ Status + data_general[[covariate_name]] + Group + (1 | PatientID), data = data_general, REML = F)
  general_test <- anova(general_model, general_model_no_interaction)
  
  #Calculate effect size
  r2_full <- r.squaredGLMM(general_model)[2]
  r2_reduced <- r.squaredGLMM(general_model_no_interaction)[2]
  f_squared <- (r2_full - r2_reduced) / (1 - r2_full)
  
  if(f_squared < 0) {
    f_squared = 0
  }
  
  #Summarize results
  result_list <- list(pValue = general_test$"Pr(>Chisq)"[2], EffSize = f_squared)
  
  return(result_list)
}


### This function is called within the interaction_deconf_both() function
feature_association_both_univariate <- function(feature_data, covariate_data, covariate_name, type = NULL) {
  # Merge feature data with covariates
  if(type == "microbiome") {
    data_general <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  } else if(type == "metabolomics") {
    data_general <- data.frame(Feature = log(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  }
  
  #Test if the association between the microbial feature and the metabolite is dependent on a feature, while group is accounted for
  general_model <- lmer(Feature ~ Status + data_general[[covariate_name]] + Group + (1 | PatientID), data = data_general, REML = F)
  general_model_no_interaction <- lmer(Feature ~ Status + Group + (1 | PatientID), data = data_general, REML = F)
  general_test <- anova(general_model, general_model_no_interaction)
  
  #Calculate effect size
  r2_full <- r.squaredGLMM(general_model)[2]
  r2_reduced <- r.squaredGLMM(general_model_no_interaction)[2]
  f_squared <- (r2_full - r2_reduced) / (1 - r2_full)
  
  estimate_sign <- sign(summary(general_model)$coefficients[,"Estimate"]["data_general[[covariate_name]]"])
  
  if(f_squared < 0) {
    f_squared = 0
  } else if (estimate_sign < 0) {
    f_squared = -f_squared
  } else {
    f_squared = f_squared
  }
  
  #Summarize results
  result_list <- list(pValue = general_test$"Pr(>Chisq)"[2], EffSize = f_squared)
  
  return(result_list)
}


feature_association_future <- function(feature_data, covariate_data_future, covariate_name, type = NULL) {
  
  # Merge feature data with covariates
  if(type == "metabolomics") {
      data_general <- data.frame(Feature = feature_data, covariate_data_future, check.names = FALSE) %>%
        filter(!is.na(!!sym(covariate_name))) %>%
        group_by(PatientID) %>%
        reframe(across(Feature, ~ .x[Status == 1]/.x[Status == 0]), !!sym(covariate_name)) %>%
        distinct()
  } else if(type == "microbiome") {
      data_general <- data.frame(Feature = feature_data, covariate_data_future, check.names = FALSE) %>%
        filter(!is.na(!!sym(covariate_name))) %>%
        group_by(PatientID) %>%
        reframe(across(Feature, ~ .x[Status == 1]-.x[Status == 0]), !!sym(covariate_name)) %>%
        distinct()
  }
  
  #Test if a change in feature over the intervention period is associated with this future outcome while accounting for baseline values
  p <- cor.test(data_general$Feature, data_general[[covariate_name]], method = "spearman")$p.value

  #Calculate effect size
  rho <- cor.test(data_general$Feature, data_general[[covariate_name]], method = "spearman")$estimate
  
  #Summarize results
  result_list <- list(pValue = p, EffSize = rho)
  
  return(result_list)
}


feature_association_outcomes <- function(feature_data, covariate_data, covariate_name, feature_type = NULL, cov_type = NULL, reference_level = NA) {
  # Merge feature data with covariates
  if(feature_type == "metabolomics") {
    data_general <- data.frame(Feature = log(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  } else if (feature_type == "microbiome"){
    data_general <- data.frame(Feature = feature_data, covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  }
  
  if(cov_type == "metabolomics") {
    data_general[[covariate_name]] <- log(data_general[[covariate_name]])
  } else if (cov_type == "microbiome") {
    data_general[[covariate_name]] <- rank(data_general[[covariate_name]])
  }
  
  data_general$Group <- relevel(data_general$Group, ref = reference_level)
  level_name1 <- levels(covariate_data$Group)[levels(covariate_data$Group) != reference_level][1]
  level_name2 <- levels(covariate_data$Group)[levels(covariate_data$Group) != reference_level][2]
  
  opposite_level1 <- paste("Status:Feature:Group", level_name1, sep="")
  opposite_level2 <- paste("Status:Feature:Group", level_name2, sep="")
  
  general_model <- lmer(data_general[[covariate_name]] ~ Status * Feature * Group + (1 | PatientID), data = data_general, REML = F)
  general_model_no_interaction <- lmer(data_general[[covariate_name]] ~ Status + Group + Feature + Status : Group + Status : Feature + Group : Feature + (1 | PatientID), data = data_general, REML = F)
  general_test <- anova(general_model, general_model_no_interaction, test = "LRT")
  
  
  #Calculate effect size
  r2_full <- r.squaredGLMM(general_model)[2]
  r2_reduced <- r.squaredGLMM(general_model_no_interaction)[2]
  f_squared <- (r2_full - r2_reduced) / (1 - r2_full)
  
  estimate_sign1 <- sign(summary(general_model)$coefficients[,"Estimate"][opposite_level1])
  estimate_sign2 <- sign(summary(general_model)$coefficients[,"Estimate"][opposite_level2])
  
  if (is.na(estimate_sign1)) {
    f_squared1 = NA
  } else if (estimate_sign1 < 0) {
    f_squared1 = -abs(f_squared)
  } else if (estimate_sign1 > 0) {
    f_squared1 = abs(f_squared)
  }
  
  if (is.na(estimate_sign2)) {
    f_squared2 = NA
  } else if (estimate_sign2 < 0) {
    f_squared2 = -abs(f_squared)
  } else {
    f_squared2 = abs(f_squared)
  }
  
  #Summarize results
  result_list <- list(pValue = general_test$"Pr(>Chisq)"[2], level_name1 = f_squared1, level_name2 = f_squared2)
  
  return(result_list)
}


feature_mediation <- function(feature_data, covariate_data, covariate_name, type = NULL, reference_level = reference_level) {
  # Merge feature data with covariates
  if(type == "microbiome") {
    data_general <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  } else if(type == "metabolomics") {
    data_general <- data.frame(Feature = log(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate(Status = as.numeric(Status))
  }
  
  data_general$Group <- relevel(data_general$Group, ref = reference_level)
  level_name1 <- levels(data_general$Group)[levels(data_general$Group) != reference_level][1]
  level_name2 <- levels(data_general$Group)[levels(data_general$Group) != reference_level][2]
  
  result_list <- list()
  
  for(comparison in c(level_name1, level_name2)) {
    data_general_analysis <- data_general %>% filter(Group %in% c(comparison, reference_level))

      model.m <- lmer(Feature ~ Status * Group + (1|PatientID), data = data_general_analysis, REML = F)
      model.y <- lmer(data_general_analysis[[covariate_name]] ~ Feature + Status * Group + (1|PatientID), data = data_general_analysis, REML = F)
      mediation_result <- mediate(model.m , model.y, treat = "Group", mediator = "Feature", time = "Status", control.value = comparison, treat.value = reference_level)
      result_list[[paste("p.Val",comparison, sep = "_")]] <- mediation_result$d0.p
      result_list[[paste("Proportion",comparison, sep = "_")]] <- abs(mediation_result$d1) / (abs(mediation_result$d1) + abs(mediation_result$z1))
  }
  
  data_general_analysis <- data_general %>% filter(Group %in% c(level_name1, level_name2))
  
  model.m <- lmer(Feature ~ Status * Group + (1|PatientID), data = data_general_analysis, REML = F)
  model.y <- lmer(data_general_analysis[[covariate_name]] ~ Feature + Status * Group + (1|PatientID), data = data_general_analysis, REML = F)
  mediation_result <- mediate(model.m , model.y, treat = "Group", mediator = "Feature", time = "Status", control.value = level_name1, treat.value = level_name2)
  result_list[[paste("p.Val",paste(level_name1, level_name2, sep = "_"), sep = "_")]] <- mediation_result$d0.p
  result_list[[paste("Proportion",paste(level_name1, level_name2, sep = "_"), sep = "_")]] <- abs(mediation_result$d1) / (abs(mediation_result$d1) + abs(mediation_result$z1))
  

  return(result_list)
}




analyze_feature <- function(feature_data, covariate_name, covariate_data, p_limit = 0.05, visit = F) {
  #create dataframe specific for covariate (without its NA rows for example)
  if(type == "microbiome") {
    data <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name))
  } else if(type == "metabolomics") {
    data <- data.frame(Feature = log(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name))
  } else if (type == "outcome") {
    #data <- data.frame(Feature = rank(feature_data), covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name)) %>% mutate("{covariate_name}" := rank(!!sym(covariate_name)))
    data <- data.frame(Feature = feature_data, covariate_data, check.names = FALSE) %>% drop_na(all_of(covariate_name))
  }
  
  if(length(unique(data[[covariate_name]])) <=1) {
    return("NS")
  } else {
    # Fit the model with that potentially trimmed data with and without the covariate term as interaction
    if (visit) {
      full_model <- lmer(Feature ~ Status + (1 | PatientID) + data[[covariate_name]] + Status:data[[covariate_name]], data = data, REML = F)
      no_group_interaction_model <- lmer(Feature ~ (1 | PatientID) + data[[covariate_name]], data = data, REML = F)
      no_covariate_interaction_model <- lmer(Feature ~ Status + (1 | PatientID), data = data, REML = F)
    } else {
      full_model <- lmer(Feature ~ Status + Group + Status:Group + (1 | PatientID) + data[[covariate_name]] + Status:data[[covariate_name]], data = data, REML = F)
      no_group_interaction_model <- lmer(Feature ~ Status + Group + (1 | PatientID) + data[[covariate_name]] + Status:data[[covariate_name]], data = data, REML = F)
      no_covariate_interaction_model <- lmer(Feature ~ Status + Group + Status:Group + (1 | PatientID) + data[[covariate_name]], data = data, REML = F)
    }

    # Compare models
    interaction_test <- anova(full_model, no_group_interaction_model)
    covariate_test <- anova(full_model, no_covariate_interaction_model)
      
    # Return results
    GroupTest = interaction_test$"Pr(>Chisq)"[2]
    CovariateTest = covariate_test$"Pr(>Chisq)"[2]
      
    if (GroupTest < p_limit & CovariateTest < p_limit){
        return("Both")
    } else if (GroupTest > p_limit & CovariateTest > p_limit){
        return("None")
    } else if (GroupTest < p_limit & CovariateTest > p_limit) {
      if(visit) {
        return("Visit")
      } else {
        return("Group")
      }
    } else if(GroupTest > p_limit & CovariateTest < p_limit) {
        return("Confounding")
    }
  } 
}


feature_posthoc_tests <- function(feature, countdata, outcome = F, outcome_name = NULL, type = NULL) {
  if (outcome) {
    if(type == "metabolomics") {
      model_function <- as.formula(paste(outcome_name, paste(paste("log(", feature, ")", sep = ""), "Visit*Group + (1 | PatientID)", sep = " * "), sep = " ~ "))
      lmer_model <- lmer(model_function, data = countdata, REML = FALSE)
      specs_formula <- as.formula(paste("pairwise ~ Group | Visit", paste("log(", feature, ")", sep = ""), sep = "*"))
      emmeans_results <- emmeans(lmer_model, specs = specs_formula, adjust = "mvt")
      comparisons <- contrast(contrast(emmeans_results, "revpairwise", by = "Group"), "revpairwise", by = NULL)
      s <- summary(comparisons)
    } else if (type == "microbiome") {
      model_function <- as.formula(paste(outcome_name, paste(feature, "Visit*Group + (1 | PatientID)", sep = " * "), sep = " ~ "))
      lmer_model <- lmer(model_function, data = countdata, REML = FALSE)
      specs_formula <- as.formula(paste("pairwise ~ Group | Visit", feature, sep = "*"))
      emmeans_results <- emmeans(lmer_model, specs = specs_formula, adjust = "dunn")
      comparisons <- contrast(contrast(emmeans_results, "revpairwise", by = "Group"), "revpairwise", by = NULL)
      s <- summary(comparisons)
    }
  } else {
    if(type == "metabolomics") {
      model_function <- as.formula(paste(paste(paste("log(", feature, sep = ""), ")", sep = ""), "Visit*Group + (1 | PatientID)", sep = " ~ "))
      lmer_model <- lmer(model_function, data = countdata, REML = FALSE)
      emmeans_results <- emmeans(lmer_model, specs = ~ Visit | Group, adjust = "mvt")
      comparisons <- contrast(contrast(emmeans_results, "revpairwise", by = "Group"), "revpairwise", by = NULL)
      s <- summary(comparisons)
    } else if (type == "microbiome") {
      model_function <- as.formula(paste(paste(paste("rank(", feature, sep = ""), ")", sep = ""), "Visit*Group + (1 | PatientID)", sep = " ~ "))
      lmer_model <- lmer(model_function, data = countdata, REML = FALSE)
      emmeans_results <- emmeans(lmer_model, specs = ~ Visit | Group, adjust = "dunn")
      comparisons <- contrast(contrast(emmeans_results, "revpairwise", by = "Group"), "revpairwise", by = NULL)
      s <- summary(comparisons)
    }
  }
  
  # Extract p-values
  p_values <- s$p.value
  estimates <- data.frame(eff_size(comparisons, sigma = sigma(lmer_model), edf = df.residual(lmer_model), method = "identity"))$effect.size
  
  # Create a dataframe with feature name and p-values
  result <- data.frame(
    Feature = feature,
    Comparison = data.frame(s)$contrast,
    P_Value = p_values,
    Estimate = estimates
  )
  
  return(result)
}


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


calculate_permanovas <- function(metabolite_df, sampleIDs){
  interaction_metab <- metabolite_df %>%
    filter(Sample.Identification %in% sampleIDs) %>%
    column_to_rownames("Sample.Identification") %>%
    mutate(across(everything(), ~ as.numeric(.)))
  
  interaction_micro <- data.frame(t(microbiome_data))[rownames(interaction_metab),] #only include the samples that were included in metabolomics data
  
  micro_dist <- vegdist(interaction_micro, method = "bray") #calculate bray curtis distance between samples based on microbiome counts
  metab_dist <- vegdist(interaction_metab, method = "mahalanobis") #calculate mahalanobis distance based on metabolite AUCs 
  
  micro_pcoa <- data.frame(cmdscale (micro_dist, k = 2)) #calculate principle components based on distances
  metab_pcoa <- data.frame(cmdscale (metab_dist, k = 2)) #calculate principle components based on distances
  
  loadings_micro <- data.frame(cor(interaction_micro, cmdscale (micro_dist, k = 2, eig = T)$points, method = "spearman")) %>%
    rownames_to_column("Feature") #calculate the correlation of the counts and the PCs to see which counts contribute most to the PCs
  loadings_metab <- data.frame(cor(interaction_metab, cmdscale (metab_dist, k = 2, eig = T)$points, method = "spearman")) %>%
    rownames_to_column("Feature") #calculate the correlation of the counts and the PCs to see which counts contribute most to the PCs
  
  return(list(interaction_metab = interaction_metab,
              interaction_micro = interaction_micro,
              micro_dist = micro_dist,
              metab_dist = metab_dist,
              micro_pcoa = micro_pcoa,
              metab_pcoa = metab_pcoa,
              loadings_micro = loadings_micro,
              loadings_metab = loadings_metab))
}


interaction_deconf <- function(features_data, covariates_data, q_limit, p_limit, p_adj_method = "fdr", type = NULL) {
  
  patterns <- c("Both", "None", "Group", "Confounding")
  covariates <- setdiff(names(covariates_data), c("PatientID", "Group", "Status"))
  results_df <- data.frame(matrix(ncol = length(covariates) + 2, nrow = length(names(features_data))))
  colnames(results_df) <- c("Feature", "Group", covariates)
  results_df$Feature <- names(features_data)
  
  results <- list(Ps = results_df, Qs = results_df, Ds = results_df, Status = results_df)
  
  for (feature in colnames(features_data)) {
    feature_data <- features_data[, feature]
    #print(feature)
    
    for (covariate in c("Group", covariates)) {
      #print(covariate)
      if (covariate!=feature) {
        association_res <- feature_association(feature_data, covariates_data, covariate, type = type, visit = F)
        results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
        results$Ds[[covariate]][results$Ds$Feature == feature] <- association_res$EffSize
      }
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  
  #Calculate if it's the covariate or the Interaction term that is significant for all features which are associated with the interaction term
  for(feature in colnames(features_data)) {

    feature_data <- features_data[, feature]
    
    if (results$Qs$Group[results$Qs$Feature == feature] < q_limit) {
      for (covariate in covariates) {
        if (covariate!=feature) {
          if (results$Qs[[covariate]][results$Qs$Feature == feature] < q_limit) {
            results$Status[[covariate]][results$Status$Feature == feature] <- analyze_feature(feature_data, covariate, covariates_data, p_limit, visit = F)
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
      
      if (results$Qs$Group[results$Qs$Feature == feature] > q_limit) {
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


interaction_deconf_both <- function(features_data, covariates_data, p_adj_method = "fdr", feature_type = NULL, univariate = F){

  features_data <- features_data %>% dplyr::select(-c("PatientID", "Status", "Group"))
  covariates <- setdiff(names(covariates_data), c("PatientID", "Status", "Group"))
  results_df <- data.frame(matrix(ncol = length(covariates) + 1, nrow = length(names(features_data))))
  colnames(results_df) <- c("Feature", covariates)
  results_df$Feature <- names(features_data)
  
  results <- list(Ps = results_df, Qs = results_df, Ds = results_df)
  
  for (feature in colnames(features_data)) {
    feature_data <- features_data[, feature]
    
    for (covariate in covariates) {
      if (univariate) {
        association_res <- feature_association_both_univariate(feature_data, covariates_data, covariate, type = feature_type)
      } else {
        association_res <- feature_association_both_multivariate(feature_data, covariates_data, covariate, type = feature_type)
      }
      results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
      results$Ds[[covariate]][results$Ds$Feature == feature] <- association_res$EffSize
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  return(results)
}


interaction_deconf_future <- function(features_data_posthoc, covariates_data_posthoc, p_adj_method = "fdr", feature_type = NULL){
  
  covariates <- covariates_data_posthoc %>% dplyr::select(-Status, -PatientID) %>% names()
  results_df <- data.frame(matrix(ncol = length(covariates) + 1, nrow = length(names(features_data_posthoc))))
  colnames(results_df) <- c("Feature", covariates)
  results_df$Feature <- names(features_data_posthoc)
  
  results <- list(Ps = results_df, Qs = results_df, Ds = results_df)
  
  for (feature in colnames(features_data_posthoc)) {
    feature_data <- features_data_posthoc[, feature]
    
    for (covariate in c(covariates)) {
      association_res <- feature_association_future(feature_data, covariates_data_posthoc, covariate, type = feature_type)
      results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
      results$Ds[[covariate]][results$Ds$Feature == feature] <- association_res$EffSize
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  return(results)
}


interaction_outcome_check <- function(multivariate_result_list, features_data, covariates_data, relevant_outcomes, sig.level = 0.1, p.adjust.methods) {
  relevant_metabolites <- multivariate_result_list$Qs %>%
    filter(if_any(!Feature, ~. < sig.level)) %>%
    dplyr::select(Feature, where(~ any(. < sig.level))) %>%
    pull(Feature)
  
  relevant_taxa <- multivariate_result_list$Qs %>%
    filter(if_any(!Feature, ~. < sig.level)) %>%
    dplyr::select(where(~ any(. < sig.level))) %>%
    names()
  
  merged_data <- merge(merge(features_data, covariates_data, by = 0) %>% rename(Tube_ID = Row.names), merge(SampleIDs %>% dplyr::select(PatientID, Tube_ID, Visit), metadata, by = c("PatientID", "Visit")), by = "Tube_ID") %>% dplyr::select(all_of(c(relevant_metabolites, relevant_taxa, relevant_outcomes)), PatientID, Visit, Group)
  
  results_df <- data.frame(matrix(ncol = length(relevant_metabolites), nrow = length(relevant_taxa))) %>% mutate(Taxon = relevant_taxa)
  colnames(results_df) <- c(relevant_metabolites, "Taxon")
  
  result_list_p <- setNames(replicate(length(relevant_outcomes), results_df, simplify = FALSE), relevant_outcomes)
  result_list_d <- setNames(replicate(length(relevant_outcomes), results_df, simplify = FALSE), relevant_outcomes)
  
  for (metabolite in relevant_metabolites) {
    for (taxon in relevant_taxa) {
      
      q <- multivariate_result_list$Qs %>%
        filter(Feature == metabolite) %>%
        pull(taxon) # Is their interaction significant
      
      if(q < 0.1) {
        for (outcome in relevant_outcomes) {
          #print(outcome)
          model_data <- merged_data %>% drop_na(all_of(outcome)) %>% mutate(PatientID = as.factor(PatientID)) #%>% mutate(across(all_of(outcome), ~sample(.)))
          f_term <- as.formula(paste0(outcome, "~rank(", taxon, ")*log(", metabolite, ")*", "Group + Visit + (1|PatientID)"))
          f_term_wo <- as.formula(paste0(outcome, "~rank(", taxon, ")+log(", metabolite, ")+rank(",taxon, "):Group+log(", metabolite, "):Group+rank(", taxon, "):log(", metabolite,  ")+Group + Visit + (1|PatientID)"))
          
          model_base <- lmer(f_term, data = model_data, REML = F)
          model_wo <- lmer(f_term_wo, data = model_data, REML = F)
          
          p <- anova(model_base, model_wo)$"Pr(>Chisq)"[2]
          
          r2_full <- r.squaredGLMM(model_base)[1]
          r2_reduced <- r.squaredGLMM(model_wo)[1]
          f_squared <- (r2_full - r2_reduced) / (1 - r2_full)
          
          if(f_squared < 0) {
            f_squared = 0
          }
          
          #print(p)
          result_list_p[[outcome]][[metabolite]][result_list_p[[outcome]]$Taxon == taxon] <- p
          result_list_d[[outcome]][[metabolite]][result_list_d[[outcome]]$Taxon == taxon] <- f_squared
        }
      }
    }
  }
  
  final_result_p <- bind_rows(lapply(relevant_outcomes, function(outcome){
    result_list_p[[outcome]] %>%
      pivot_longer(cols = -Taxon, names_to = "Metabolite", values_to = "PVal") %>%
      drop_na(PVal) %>%
      mutate(PValadj = p.adjust(PVal, method = p.adjust.methods)) %>%
      filter(PValadj <= sig.level) %>%
      mutate(Outcome = outcome)
    }))
  
  final_result_d <- bind_rows(lapply(relevant_outcomes, function(outcome){
    result_list_d[[outcome]] %>%
      pivot_longer(cols = -Taxon, names_to = "Metabolite", values_to = "DVal") %>%
      mutate(Outcome = outcome)}))
  
  final_result <- merge(final_result_p, final_result_d, by = c("Taxon", "Metabolite", "Outcome"))
  
  return(final_result)
}


interaction_deconf_diffs <- function(features_data, covariates_data, q_limit, p_limit, p_adj_method = "fdr", type = NULL) {
  
  patterns <- c("Both", "None", "Group", "Confounding")
  covariates <- setdiff(names(covariates_data), c("PatientID", "Group"))
  results_df <- data.frame(matrix(ncol = length(covariates) + 2, nrow = length(names(features_data))))
  colnames(results_df) <- c("Feature", "Group", covariates)
  results_df$Feature <- names(features_data)
  
  results <- list(Ps = results_df, Qs = results_df, Ds = results_df, Status = results_df)
  
  for (feature in colnames(features_data)) {
    feature_data <- features_data[, feature]
    
    for (covariate in c("Group", covariates)) {
      #print(covariate)
      association_res <- feature_association(feature_data, covariates_data, covariate, type = type)
      results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
      results$Ds[[covariate]][results$Ds$Feature == feature] <- association_res$EffSize
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  
  #Calculate if it's the covariate or the Interaction term that is significant for all features which are associated with the interaction term
  for(feature in colnames(features_data)) {
    
    feature_data <- features_data[, feature]
    
    if (results$Qs$Group[results$Qs$Feature == feature] < q_limit) {
      for (covariate in covariates) {
        if (results$Qs[[covariate]][results$Qs$Feature == feature] < q_limit) {
          results$Status[[covariate]][results$Status$Feature == feature] <- analyze_feature(feature_data, covariate, covariates_data, p_limit)
        } else {
          results$Status[[covariate]][results$Status$Feature == feature] <- "NS"
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


interaction_deconf_visit <- function(features_data, covariates_data, p_limit, p_adj_method = "fdr", type = NULL) {
  
  patterns <- c("Both", "None", "Group", "Confounding")
  covariates <- setdiff(names(covariates_data), c("PatientID"))
  results_df <- data.frame(matrix(ncol = length(covariates), nrow = length(names(features_data))))
  colnames(results_df) <- covariates
  results_df$Feature <- names(features_data)
  
  results <- list(Ps = results_df, Qs = results_df, Ds = results_df, Status = results_df)
  
  for (feature in colnames(features_data)) {
    feature_data <- features_data[, feature]
    #print(feature)
    association_res <- feature_association(feature_data, covariates_data, "Status", type = type, visit = T)
    results$Ps[["Status"]][results$Ps$Feature == feature] <- association_res$pValue
    results$Ds[["Status"]][results$Ds$Feature == feature] <- association_res$EffSize
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  
  #Calculate if it's the covariate or the Interaction term that is significant for all features which are associated with the interaction term
  for(feature in colnames(features_data)) {
    
    feature_data <- features_data[, feature]
    
    if (results$Qs$Status[results$Qs$Feature == feature] < p_limit) {
      for (covariate in covariates[covariates != "Status"]) {
        if (covariate!=feature) {
          results$Status[[covariate]][results$Status$Feature == feature] <- analyze_feature(feature_data, covariate, covariates_data, p_limit, visit = T)
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
      
      if (results$Qs$Status[results$Qs$Feature == feature] > p_limit) {
        results$Status$Status[results$Status$Feature == feature] <- "NS"
      } else if(Group_conf == T) {
        results$Status$Status[results$Status$Feature == feature] <- paste(Conf_variables, collapse = ", ")
      } else if (Group_conf == F) {
        results$Status$Status[results$Status$Feature == feature] <- "Deconfounded"
      }
    }
  }
  
  return(results)
}


interaction_outcome_assoc <- function(features_data, covariates_data, p_adj_method = "fdr", feature_type = NULL, cov_type = NULL, reference_level = NA) {
  
  covariates <- setdiff(names(covariates_data), c("PatientID", "Status", "Group"))
  results_df <- data.frame(matrix(ncol = length(covariates) + 1, nrow = length(names(features_data))))
  colnames(results_df) <- c("Feature", covariates)
  results_df$Feature <- names(features_data)
  level_name1 <- levels(covariates_data$Group)[levels(covariates_data$Group) != reference_level][1]
  level_name2 <- levels(covariates_data$Group)[levels(covariates_data$Group) != reference_level][2]
  
  results <- list(Ps = results_df, Qs = results_df, level_name1 = results_df, level_name2 = results_df)
  
  for (feature in colnames(features_data)) {
    feature_data <- features_data[, feature]
    print(feature)
    for (covariate in covariates) {
      print(covariate)
      association_res <- feature_association_outcomes(feature_data, covariates_data, covariate, feature_type = feature_type, cov_type = cov_type, reference_level = reference_level)
      results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
      results$level_name1[[covariate]][results$level_name1$Feature == feature] <- association_res$level_name1
      results$level_name2[[covariate]][results$level_name2$Feature == feature] <- association_res$level_name2
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  names(results)[3] <- level_name1
  names(results)[4] <- level_name2
  
  return(results)
}


mediation_analysis <- function(features_data, covariates_data, p_adj_method = "fdr", feature_type = NULL, reference_level = NULL){
  
  features_data <- features_data %>% dplyr::select(-c("PatientID", "Status", "Group"))
  covariates <- setdiff(names(covariates_data), c("PatientID", "Status", "Group"))
  
  results_df <- data.frame(matrix(ncol = length(covariates) + 1, nrow = length(names(features_data))))
  colnames(results_df) <- c("Feature", covariates)
  results_df$Feature <- names(features_data)
  
  results_1 <- list(Ps = results_df, Qs = results_df, Ds = results_df)
  results_2 <- results_1
  results_3 <- results_1
  
  for (feature in colnames(features_data)) {
    feature_data <- features_data[, feature]
    
    for (covariate in covariates) {
      association_res <- feature_mediation(feature_data, covariates_data, covariate, type = feature_type, reference_level = reference_level)
      results_1$Ps[[covariate]][results_1$Ps$Feature == feature] <- association_res[[1]]
      results_1$Ds[[covariate]][results_1$Ds$Feature == feature] <- association_res[[2]]
      results_2$Ps[[covariate]][results_2$Ps$Feature == feature] <- association_res[[3]]
      results_2$Ds[[covariate]][results_2$Ds$Feature == feature] <- association_res[[4]]
      results_3$Ps[[covariate]][results_3$Ps$Feature == feature] <- association_res[[5]]
      results_3$Ds[[covariate]][results_3$Ds$Feature == feature] <- association_res[[6]]
    }
  }
  
  #p-adjust results file
  results_1[["Qs"]] <- data.frame(results_1[["Ps"]]) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  results_2[["Qs"]] <- data.frame(results_2[["Ps"]]) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  results_3[["Qs"]] <- data.frame(results_3[["Ps"]]) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  results <- list(results_1, results_2, results_3) %>%
    setNames(sub(".*_", "", names(association_res)[c(1, 3)]))
  
  return(results)
}


# Big deconf function
interaction_outcome_deconf <- function(features_data_outcome, covariates_data_outcome, p_adj_method = "fdr", feature_type = NULL, cov_type = NULL, reference_level = NA, q_limit = 0.1, p_limit = 0.05, other_variables = c()) {
  
  covariates <- setdiff(names(covariates_data_outcome), c("PatientID", "Status", "Group"))
  results_df <- data.frame(matrix(ncol = length(covariates) + 1, nrow = length(names(features_data_outcome))))
  colnames(results_df) <- c("Feature", covariates)
  results_df$Feature <- names(features_data_outcome)
  level_name1 <- levels(covariates_data_outcome$Group)[levels(covariates_data_outcome$Group) != reference_level][1]
  level_name2 <- levels(covariates_data_outcome$Group)[levels(covariates_data_outcome$Group) != reference_level][2]
  
  results <- list(Ps = results_df, Qs = results_df, level_name1 = results_df, level_name2 = results_df)
  
  for (feature in colnames(features_data_outcome)) {
    feature_data <- features_data_outcome[, feature]
    print(feature)
    for (covariate in covariates) {
      #print(covariate)
      association_res <- feature_association_outcomes(feature_data, covariates_data_outcome, covariate, feature_type = feature_type, cov_type = cov_type, reference_level = reference_level)
      results$Ps[[covariate]][results$Ps$Feature == feature] <- association_res$pValue
      results$level_name1[[covariate]][results$level_name1$Feature == feature] <- association_res$level_name1
      results$level_name2[[covariate]][results$level_name2$Feature == feature] <- association_res$level_name2
    }
  }
  
  #p-adjust results file
  results$Qs <- data.frame(results$Ps) %>%
    mutate(across(-Feature,~ p.adjust(., method = p_adj_method)))
  
  names(results)[3] <- level_name1
  names(results)[4] <- level_name2
  
  results$Qs <- data.frame(results$Ps) %>%
    gather(Genus, P, -Feature) %>%
    mutate(across(P, ~ p.adjust(., method = p_adj_method))) %>%
    spread(Genus, P)
  
  sig_outcomes <- names(results$Qs %>%
                          dplyr::select(Feature, where(~ any(. < q_limit))) %>%
                          filter(if_any(everything(), ~ . < q_limit)) %>%
                          select(-Feature))
  
  #Calculate association between outcomes and covariates without interactions
  result_covs <- interaction_deconf(covariates_data %>% dplyr::select(contains(c("_igf_", "NFL", "mri_t2", "lesions", "edss", "bdi2_score_corrected", "sdmt")),-contains(".v5")), covariates_data %>% dplyr::select(-contains(c(".v5", "mri_pbvc", "_igf_", "NFL", "mri_t2", "lesions", "edss", "bdi2_score_corrected", "sdmt"))), p_limit, p_adj_method, type = "outcome")
  result_covs_filtered <- result_covs$Qs %>%
    filter(Feature %in% sig_outcomes) %>%
    dplyr::select(Feature, where(~ any(. < 0.1)))
  potential_confounders <- c(names(result_covs_filtered %>% dplyr::select(-Feature)), other_variables)
  
  #Calculate if it's the covariate or the Interaction term that is significant for all features which are associated with the interaction term
  for(outcome_var in sig_outcomes) {
    
    #Create dataframe to add to the results list for each outcome variable that is associated with a feature
    status_df <- data.frame(matrix(nrow = length(colnames(features_data_outcome)), ncol = length(potential_confounders)+1))
    colnames(status_df) <- c(outcome_var,potential_confounders)
    status_df$Feature <- colnames(features_data_outcome)
    results[[outcome_var]] <- status_df
    
    for(feature in colnames(features_data_outcome)) {
      
      feature_data <- features_data_outcome[, feature]
      
      if (results$Qs[results$Qs$Feature == feature, outcome_var] < q_limit) {
        
        for (confounder in potential_confounders) {
          
          if (result_covs$Qs[result_covs$Qs$Feature == outcome_var, confounder] < q_limit) {
            results[[outcome_var]][results[[outcome_var]]$Feature == feature, confounder] <- analyze_feature_outcome_deconf(feature_data, outcome_var, confounder, covariates_data, p_limit = p_limit, feature_type)
          } else {
            results[[outcome_var]][results[[outcome_var]]$Feature == feature, confounder] <- "NS"
          }
          
        }
        
        
        Group_conf <- results[[outcome_var]] %>%
          filter(Feature == feature) %>%
          rowwise() %>%
          mutate(has_confounded = any(c_across(-Feature) == "Confounding", na.rm = TRUE)) %>%
          ungroup() %>%
          pull(has_confounded) %>%
          any(na.rm = TRUE)
        
        Conf_variables <- results[[outcome_var]] %>%
          filter(Feature == feature) %>%
          summarise(across(-Feature, ~. == "Confounding")) %>%
          dplyr::select(where(~any(. == TRUE, na.rm = T))) %>%
          names()
        
        if (results$Qs[results$Qs$Feature == feature, outcome_var] > q_limit) {
          results[[outcome_var]][[outcome_var]][results[[outcome_var]]$Feature == feature] <- "NS"
        } else if(Group_conf == T) {
          results[[outcome_var]][[outcome_var]][results[[outcome_var]]$Feature == feature] <- paste(Conf_variables, collapse = ", ")
        } else if (Group_conf == F) {
          results[[outcome_var]][[outcome_var]][results[[outcome_var]]$Feature == feature] <- "Deconfounded"
        }
      }
    }
  }
  
  return(results)
}


# Deconf Test
analyze_feature_outcome_deconf <- function(feature_data, outcome_var, confounder, covariates_data, p_limit = 0.05, feature_type = NULL) {
  #create dataframe specific for covariate (without its NA rows for example)
  if(feature_type == "metabolomics") {
    data_deconf <- data.frame(Feature = log(feature_data), covariates_data, check.names = FALSE) %>% drop_na(all_of(outcome_var)) %>% drop_na(all_of(confounder)) %>% mutate(Status = as.numeric(Status))
  } else if (feature_type == "microbiome"){
    data_deconf <- data.frame(Feature = feature_data, covariates_data, check.names = FALSE) %>% drop_na(all_of(outcome_var)) %>% drop_na(all_of(confounder)) %>% mutate(Status = as.numeric(Status))
  }
  
  
  if(length(unique(data_deconf[[outcome_var]])) <=1) {
    return("NS")
  } else {
    # Fit the model with that potentially trimmed data with and without the covariate term as interaction
    data_deconf$Group <- relevel(data_deconf$Group, ref = reference_level)
    level_name1 <- levels(covariates_data$Group)[levels(covariates_data$Group) != reference_level][1]
    level_name2 <- levels(covariates_data$Group)[levels(covariates_data$Group) != reference_level][2]
    
    opposite_level1 <- paste("Status:Feature:Group", level_name1, sep="")
    opposite_level2 <- paste("Status:Feature:Group", level_name2, sep="")
    
    full_model <- lmer(data_deconf[[outcome_var]] ~ Status * Feature * Group + Status*data_deconf[[confounder]] + (1 | PatientID), data = data_deconf, REML = F)
    no_covariate_interaction_model <- lmer(data_deconf[[outcome_var]] ~ Status * Feature * Group + (1 | PatientID), data = data_deconf, REML = F)
    no_group_interaction_model <- lmer(data_deconf[[outcome_var]] ~ Status + Group + Feature + Status : Group + Status : Feature + Group : Feature + Status:data_deconf[[confounder]] + (1 | PatientID), data = data_deconf, REML = F)
    
    # Compare models
    interaction_test <- anova(full_model, no_group_interaction_model, test = "LRT")
    covariate_test <- anova(full_model, no_covariate_interaction_model, test = "LRT")
    
    # Return results
    GroupTest = interaction_test$"Pr(>Chisq)"[2]
    CovariateTest = covariate_test$"Pr(>Chisq)"[2]
    
    if (GroupTest < p_limit & CovariateTest < p_limit){
      return("Both")
    } else if (GroupTest > p_limit & CovariateTest > p_limit){
      return("None")
    } else if (GroupTest < p_limit & CovariateTest > p_limit) {
      return("Outcome")
    } else if(GroupTest > p_limit & CovariateTest < p_limit) {
      return("Confounding")
    }
  } 
}


