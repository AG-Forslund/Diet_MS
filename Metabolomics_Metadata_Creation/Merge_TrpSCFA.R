library(here)
library(dplyr)
library(stringr)

type = "Stool"

if(type == "Serum"){
  SCFA_file <- readRDS(here("Results", "Metabolomics","Serum", "SCFA", "Filtered_Normalized_BatchCorrected_AreaNormISTD_NAMS_SCFA_Serum.RDS"))
  Trp_file <- readRDS(here("Results", "Metabolomics","Serum", "Tryptophan", "Run2", "Filtered_Normalized_BatchCorrected_AreaNormISTD_NAMS_Trp_Serum.RDS"))
  merged_file <- merge(SCFA_file, Trp_file %>% select(-"Sample.Type"), by = "Sample.Identification") %>%
    select(-c(Sequence.Position.x, Sequence.Position.y)) %>%
    rename(Batch.SCFA = Batch.x) %>%
    rename(Batch.Trp = Batch.y)
  write_rds(merged_file, here("Results", "Metabolomics","Serum","NAMS_Serum_AreaNormISTD_data.RDS"))
} else {
  SCFA_file <- readRDS(here("Results", "Metabolomics","Stool", "Filtered_Normalized_BatchCorrected_AreaNormISTD_NAMS_SCFA_Stool.RDS"))
  Trp_file <- readRDS(here("Results", "Metabolomics","Stool", "Filtered_Normalized_BatchCorrected_AreaNormISTD_NAMS_Trp_Stool.RDS"))
  merged_file <- merge(SCFA_file, Trp_file %>% select(-"Sample.Type"), by = "Sample.Identification") %>%
    select(-c(Sequence.Position.x, Sequence.Position.y)) %>%
    rename(Batch.SCFA = Batch.x) %>%
    rename(Batch.Trp = Batch.y)
  write_rds(merged_file, here("Results", "Metabolomics","Stool","NAMS_Stool_AreaNormISTD_data.RDS"))
}


#----------------------------
# For Concentration results

if(type == "Serum"){
  SCFA_file <- readRDS(here("Results", "Metabolomics","Serum", "SCFA", "Filtered_Normalized_Concentration_NAMS_SCFA_Serum.RDS"))
  Trp_file <- readRDS(here("Results", "Metabolomics","Serum", "Tryptophan", "Run2", "Filtered_Normalized_BatchCorrected_Concentration_NAMS_Trp_Serum.RDS"))
  merged_file <- merge(SCFA_file, Trp_file %>% select(-"Sample.Type"), by = "Sample.Identification") %>%
    select(-c(Sequence.Position.x, Sequence.Position.y)) %>%
    rename(Batch.SCFA = Batch.x) %>%
    rename(Batch.Trp = Batch.y)
  write_rds(merged_file, here("Results", "Metabolomics","Serum","NAMS_Serum_Concentration_data.RDS"))
} else {
  SCFA_file <- readRDS(here("Results", "Metabolomics","Stool", "Filtered_Normalized_BatchCorrected_Concentration_NAMS_SCFA_Stool.RDS"))
  Trp_file <- readRDS(here("Results", "Metabolomics","Stool", "Filtered_Normalized_BatchCorrected_Concentration_NAMS_Trp_Stool.RDS"))
  merged_file <- merge(SCFA_file, Trp_file %>% select(-"Sample.Type"), by = "Sample.Identification") %>%
    select(-c(Sequence.Position.x, Sequence.Position.y)) %>%
    rename(Batch.SCFA = Batch.x) %>%
    rename(Batch.Trp = Batch.y)
  write_rds(merged_file, here("Results", "Metabolomics","Stool","NAMS_Stool_Concentration_data.RDS"))
}





