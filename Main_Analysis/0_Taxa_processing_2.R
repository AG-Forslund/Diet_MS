##Project: NAMS project
##Aim: Taxonomy data preparation
##Author: Víctor Hugo Jarquín-Díaz
##Date: 28.11.2023

##Libraries 
##Load packages 
library(tidyverse)
library(phyloseq)
library(microbiome)
library(data.table)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(ggsci)

##large size data path
large.path<- "/fast/AG_Forslund/shared/NAMS/"
repo.path<- "/fast/AG_Forslund/Victor/GitRepos/nams/"

##Load motus table
tax_mOTU<- read_tsv(file = paste0(large.path, "nams_motusv3.1.tsv"))
##Check samples with zero counts
colSums(tax_mOTU[2:ncol(tax_mOTU)])
##Check mOTUs with zero counts
rowSums(tax_mOTU[2:ncol(tax_mOTU)])
##Remove zero counts
tax_mOTU<- tax_mOTU[rowSums(tax_mOTU[2:ncol(tax_mOTU)])>0, ]
##Adjust 1st colum name to keep same structure
tax_mOTU%>%
  dplyr::rename(mOTU_ID= `...1`)->tax_mOTU

##Adjust the mess
tax_mOTU%>%
  separate("mOTU_ID", into = c("Species", "mOTUs_ID"), sep = " \\[")%>%
  dplyr::mutate(mOTUs_ID= str_replace_all(mOTUs_ID, "\\]", ""))%>%
  dplyr::mutate(mOTUs_ID= case_when(is.na(mOTUs_ID) ~ "unassigned",
                                    TRUE ~ mOTUs_ID))%>%
  dplyr::select(!"Species")%>%
  column_to_rownames("mOTUs_ID")->mOTU_table

##Load taxonomy
load(url("https://github.com/AlessioMilanese/motus_taxonomy/blob/master/data/motus_taxonomy_3.0.1.Rdata?raw=true"))
#Clear numbers from tax names
motus3.0_taxonomy[] <- lapply(motus3.0_taxonomy, gsub, pattern='[[:digit:]]+ ', replacement='')

motus3.0_taxonomy %>%
  rownames_to_column("rownames")%>%
  dplyr::select(!c("profiled",  "rownames"))%>%
  mutate(mOTUs_ID = gsub("_v3_", "_v31_", mOTUs_ID)) %>%
  column_to_rownames("mOTUs_ID")->tax_table

# some mOTUs have missing annotation
missing_annotation <- rownames(mOTU_table)[!rownames(mOTU_table) %in% rownames(tax_table)]

##Which species name are assigned at the mOTU table for these without annotation
tax_mOTU%>%
  separate("mOTU_ID", into = c("Species", "mOTUs_ID"), sep = " \\[")%>%
  dplyr::mutate(mOTUs_ID= str_replace_all(mOTUs_ID, "\\]", ""))%>%
  dplyr::mutate(mOTUs_ID= case_when(is.na(mOTUs_ID) ~ "unassigned",
                                    TRUE ~ mOTUs_ID))%>%
  dplyr::filter(mOTUs_ID%in%missing_annotation)%>%
  dplyr::mutate(Species= paste0("NA ", Species))%>%
  dplyr::select(c(mOTUs_ID, Species))->tmp

##Extract the taxonomy from other mOTUs with the same species name
tax_table%>%
  dplyr::filter(Species%in%unique(tmp$Species))%>%
  dplyr::distinct()%>%
  remove_rownames()->tmp2

##Create a taxonomy table for the mOTUs with missing annotation
tmp%>%
  plyr::join(tmp2, by = "Species")%>%
  relocate(Species, .after = last_col())%>%
  column_to_rownames("mOTUs_ID")-> missing_taxa

tax_table%>%
  rbind(., missing_taxa)-> tax_table

##clean the environment
rm(tax_mOTU, motus3.0_taxonomy, tmp, tmp2, missing_taxa, missing_annotation)

##Let's get the total number of reads taxonomically assigned
tax_mapped<- as.data.frame(colSums(mOTU_table))

##Add tax assigned reads
tax_mapped%>%
  rownames_to_column("Sample_ID")%>%
  rename("Tax_mapped"= "colSums(mOTU_table)")%>%
  relocate(Sample_ID)-> mapped_reads_stats

##Load sample data
sample_data<- readxl::read_xlsx(paste0(repo.path, "data/Upload_NAMS_DNA.xlsx"))
metadata<- vroom::vroom(paste0(repo.path, "data/SampleIDs_PerProtocol_V1V3.csv"))

sample_data%>%
  plyr::join(metadata, by="Tube_ID")-> sample_data

##To make Phyloseq object
##A - Taxonomy
##1) Use the mOTUS matrix and transform it to "OTU table" format
mOTUS<- otu_table(mOTU_table, taxa_are_rows = T)
sample_names(mOTUS)
##2) Use sample dataframe and transform it to "sample data" format
sample<- sample_data(sample_data)
sample_names(sample)<- sample_data$Tube_ID
##3) Use taxa matrix and transform it to "tax table" format
tax<-tax_table(as.matrix(tax_table))
sample_names(tax)

PS.mOTUs <- merge_phyloseq(mOTUS, tax, sample)
saveRDS(PS.mOTUs, file=paste0(repo.path, "data/PS.nams1.mOTUS3.rds")) 

##Now run this for the second batch 
##Load motus table
tax_mOTU<- read_tsv(file = paste0(large.path, "nams_2_motusv3.1.tsv"))
##Check samples with zero counts
colSums(tax_mOTU[2:ncol(tax_mOTU)])
##Check mOTUs with zero counts
rowSums(tax_mOTU[2:ncol(tax_mOTU)])
##Remove zero counts
tax_mOTU<- tax_mOTU[rowSums(tax_mOTU[2:ncol(tax_mOTU)])>0, ]
##Adjust 1st colum name to keep same structure
tax_mOTU%>%
  dplyr::rename(mOTU_ID= `...1`)->tax_mOTU

##Adjust the mess
tax_mOTU%>%
  separate("mOTU_ID", into = c("Species", "mOTUs_ID"), sep = " \\[")%>%
  dplyr::mutate(mOTUs_ID= str_replace_all(mOTUs_ID, "\\]", ""))%>%
  dplyr::mutate(mOTUs_ID= case_when(is.na(mOTUs_ID) ~ "unassigned",
                                    TRUE ~ mOTUs_ID))%>%
  dplyr::select(!"Species")%>%
  column_to_rownames("mOTUs_ID")->mOTU_table

##Load taxonomy
load(url("https://github.com/AlessioMilanese/motus_taxonomy/blob/master/data/motus_taxonomy_3.0.1.Rdata?raw=true"))
#Clear numbers from tax names
motus3.0_taxonomy[] <- lapply(motus3.0_taxonomy, gsub, pattern='[[:digit:]]+ ', replacement='')

motus3.0_taxonomy %>%
  rownames_to_column("rownames")%>%
  dplyr::select(!c("profiled",  "rownames"))%>%
  mutate(mOTUs_ID = gsub("_v3_", "_v31_", mOTUs_ID)) %>%
  column_to_rownames("mOTUs_ID")->tax_table

# some mOTUs have missing annotation
missing_annotation <- rownames(mOTU_table)[!rownames(mOTU_table) %in% rownames(tax_table)]

##Which species name are assigned at the mOTU table for these without annotation
tax_mOTU%>%
  separate("mOTU_ID", into = c("Species", "mOTUs_ID"), sep = " \\[")%>%
  dplyr::mutate(mOTUs_ID= str_replace_all(mOTUs_ID, "\\]", ""))%>%
  dplyr::mutate(mOTUs_ID= case_when(is.na(mOTUs_ID) ~ "unassigned",
                                    TRUE ~ mOTUs_ID))%>%
  dplyr::filter(mOTUs_ID%in%missing_annotation)%>%
  dplyr::mutate(Species= paste0("NA ", Species))%>%
  dplyr::select(c(mOTUs_ID, Species))->tmp

##Extract the taxonomy from other mOTUs with the same species name
tax_table%>%
  dplyr::filter(Species%in%unique(tmp$Species))%>%
  dplyr::distinct()%>%
  remove_rownames()->tmp2

##Create a taxonomy table for the mOTUs with missing annotation
tmp%>%
  plyr::join(tmp2, by = "Species")%>%
  relocate(Species, .after = last_col())%>%
  column_to_rownames("mOTUs_ID")-> missing_taxa

tax_table%>%
  rbind(., missing_taxa)-> tax_table

##clean the environment
rm(tax_mOTU, motus3.0_taxonomy, tmp, tmp2, missing_taxa, missing_annotation)

##Let's get the total number of reads taxonomically assigned
tax_mapped<- as.data.frame(colSums(mOTU_table))

##Add tax assigned reads
tax_mapped%>%
  rownames_to_column("Sample_ID")%>%
  rename("Tax_mapped"= "colSums(mOTU_table)")%>%
  relocate(Sample_ID)-> mapped_reads_stats_2

##Load sample data
sample_data_2<- readxl::read_xlsx(paste0(repo.path, "data/Upload_NAMS_DNA_2.xlsx"))
sample_data_2%>%
  plyr::join(metadata, by="Tube_ID")%>%
  dplyr::mutate(Sample_Number = str_extract(Tube_ID, "\\d{2,3}"))%>%
  dplyr::mutate(Visit= case_when(is.na(Visit)~ "V0",
                                 TRUE~Visit))-> sample_data_2

##To make Phyloseq object
##A - Taxonomy
##1) Use the mOTUS matrix and transform it to "OTU table" format
mOTUS<- otu_table(mOTU_table, taxa_are_rows = T)
sample_names(mOTUS)
##2) Use sample dataframe and transform it to "sample data" format
sample<- sample_data(sample_data_2)
sample_names(sample)<- sample_data_2$Tube_ID
##3) Use taxa matrix and transform it to "tax table" format
tax<-tax_table(as.matrix(tax_table))
sample_names(tax)

PS.mOTUs2 <- merge_phyloseq(mOTUS, tax, sample)
saveRDS(PS.mOTUs2, file=paste0(repo.path, "data/PS.nams2.mOTUS3.rds")) 

merge_phyloseq(PS.mOTUs, PS.mOTUs2)-> PS.mOTUS.all
saveRDS(PS.mOTUS.all, file=paste0(repo.path, "data/PS.nams.all.mOTUS3.rds")) 

##Clean environment
# Define the list of names of relevant objects to keep
To_keep <- c("PS.mOTUs", "PS.mOTUs2", "PS.mOTUS.all",
             "sample_data", "sample_data_2",
             "mapped_reads_stats", "mapped_reads_stats_2", 
             "large.path", "repo.path")

# Use the list parameter of rm() to remove all objects except the ones in objects_to_keep
rm(list = setdiff(ls(), To_keep))
