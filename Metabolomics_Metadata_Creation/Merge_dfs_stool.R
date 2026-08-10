#Merge dataframes
library(here)
library(dplyr)
library(stringr)

setwd(here("Results", "Metabolomics", "Stool"))

number_pools <- function(x){
  y <- x %>%
    mutate(Sample.Type = ifelse(grepl("Pool._1",Sample.Identification), "Pooled QC 1", Sample.Type)) %>%
    mutate(Sample.Type = ifelse(grepl("Pool._2",Sample.Identification), "Pooled QC 2", Sample.Type)) %>%
    mutate(Sample.Type = ifelse(grepl("Pool._3",Sample.Identification), "Pooled QC 3", Sample.Type))
  return(y)
}



batch1 <- read.csv("20230808_scFAs_stool_data_long_corr_Batch1.csv") %>%
  mutate(Sample.Type = ifelse(grepl("Pool",Sample.Identification), "Pooled QC 1", Sample.Type))

batch2 <- read.csv("20230808_scFAs_stool_data_long_corr_Batch2.csv") %>%
  mutate(Sample.Identification = ifelse(grepl("Standby|CP|Pool|Equilibration|QC|Process|Solven",Sample.Identification), paste("B2", Sample.Identification, sep = "_"), Sample.Identification)) %>%
  number_pools() %>%
  mutate(Batch = "Batch2")

batch3 <- read.csv("20230808_scFAs_stool_data_long_corr_Batch3.csv") %>%
  mutate(Sample.Identification = ifelse(grepl("Standby|CP|Pool|Equilibration|QC|Process|Solven",Sample.Identification), paste("B3", Sample.Identification, sep = "_"), Sample.Identification)) %>%
  number_pools() %>%
  mutate(Batch = "Batch3")

merged_file <- rbind(batch1, batch2, batch3)

write.csv(merged_file, "NAMS_SCFA_Stool_long_080823.csv")



