#Merge dataframes
library(here)
library(dplyr)
library(stringr)

setwd(here("Results", "Metabolomics", "Serum", "SCFA"))

batch1 <- read.csv("20230308_scFAs_data_long_corr_Batch1.csv")

batch2 <- read.csv("20230308_scFAs_data_long_corr_Batch2.csv") %>%
  mutate(Sample.Identification = ifelse(grepl("Standby|CP|Pool|Equilibration|QC|Process|Solven",Sample.Identification), paste("B2", Sample.Identification, sep = "_"), Sample.Identification)) %>%
  mutate(Batch = "Batch2")

batch3 <- read.csv("20230308_scFAs_data_long_corr_Batch3.csv") %>%
  mutate(Sample.Identification = ifelse(grepl("Standby|CP|Pool|Equilibration|QC|Process|Solven",Sample.Identification), paste("B3", Sample.Identification, sep = "_"), Sample.Identification)) %>%
  mutate(Batch = "Batch3")


merged_file <- rbind(batch1, batch2, batch3)

write.csv(merged_file, "NAMS_SCFA_Serum_long_040823.csv")



