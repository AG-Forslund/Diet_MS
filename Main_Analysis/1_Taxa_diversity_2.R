##Project: NAMS project
##Aim: Taxonomy data preparation
##Author: Víctor Hugo Jarquín-Díaz
##Date: 28.11.2023

##Libraries 
##Load packages 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(ggsci)
library(vegan)

##Repository location
repo.path<- "/fast/AG_Forslund/Victor/GitRepos/nams/"

ReRun<- T
if(ReRun){
  source(paste0(repo.path, "scripts/0_Taxa_processing_2.R"))
}

##Rarefaction
## Rarefy without replacement to the min sequencing depth 
vegan::rarecurve(t(otu_table(PS.mOTUs)), step=20, cex=0.5,  label =T)
abline(v=4000, col="red")

PS.Rare<- prune_samples(!is.na(PS.mOTUs@sam_data$Group), PS.mOTUs)
PS.Rare<- rarefy_even_depth(PS.Rare, rngseed=2020, sample.size=4000, replace=F)

#--------------
### Addition from Friederike
rare_df <- PS.mOTUs %>%
  otu_table() %>%
  `@`(".Data") %>%
  t() %>%
  vegan::rarecurve(step = 20, tidy = TRUE)

plot_rarefy <- ggplot(rare_df, aes(x = Sample, y = Species, group = Site, colour = Site)) +
  geom_line(linewidth = 0.6) +
  labs(x = "Sequencing depth", y = "Number of OTUs") +
  geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
  scale_color_viridis_d(option = "plasma") +
  theme_classic() +
  theme(legend.position = "none")

#ggsave(here("figures", "Rarefaction_plot.pdf"))
#--------------


##Estimate alpha diversity
motus.alpha.div <- microbiome::alpha(PS.mOTUs, index =c("Observed","Chao1", "Shannon", "Simpson")) 
motus.alpha.div<- as.data.frame(motus.alpha.div)
motus.alpha.div%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::left_join(sample_data, by= "Tube_ID")-> motus.alpha.div

#--------------
### Friederike Addition
motus.alpha.div <- microbiome::alpha(PS.mOTUs) 
motus.alpha.div<- as.data.frame(motus.alpha.div)
motus.alpha.div <- motus.alpha.div%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::left_join(SampleIDs, by= "Tube_ID") %>%
  drop_na()

alphadiffs <- motus.alpha.div %>%
  group_by(PatientID) %>%
  summarise(across(all_of("observed"), ~ (.x[Visit == "V3"]/.x[Visit == "V0"]-1)*100)) %>%
  rename(diversity_shannon = observed) %>% # this is just because of the plotting function
  merge(., SampleIDs %>% dplyr::select(PatientID, Group) %>% unique(), by = "PatientID") %>%
  mutate(Group = as.factor(Group))

adp_diffs <- alpha_div_plot(alphadiffs, "Richness")
#ggsave(here("figures", "alpha_diff_raw.pdf"),adp_diffs)
#--------------

tax_mapped%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::rename(mapped= `colSums(mOTU_table)`)%>%
  dplyr::left_join(motus.alpha.div, by= "Tube_ID")-> motus.alpha.div

##Rarefied alpha diversity
rare.alpha.div <- microbiome::alpha(PS.Rare, index =c("Observed","Chao1", "Shannon", "Simpson")) 
rare.alpha.div<- as.data.frame(rare.alpha.div)
rare.alpha.div%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::left_join(sample_data, by= "Tube_ID")-> rare.alpha.div

#--------------
### Friederike Addition
rare.alpha.div <- microbiome::alpha(PS.Rare) 
rare.alpha.div<- as.data.frame(rare.alpha.div)
rare.alpha.div <- rare.alpha.div%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::left_join(SampleIDs, by= "Tube_ID") %>%
  drop_na()

alphadiffs_rare <- rare.alpha.div %>%
  group_by(PatientID) %>%
  summarise(across(all_of("observed"), ~ (.x[Visit == "V3"]/.x[Visit == "V0"]-1)*100)) %>%
  rename(diversity_shannon = observed) %>% # this is just because of the plotting function
  merge(., SampleIDs %>% dplyr::select(PatientID, Group) %>% unique(), by = "PatientID") %>%
  mutate(Group = as.factor(Group))

adp_diffs_rare <- alpha_div_plot(alphadiffs_rare, "Richness")
#ggsave(here("figures", "alpha_diff_raw.pdf"),adp_diffs)

alphadiffs_merged <- merge(alphadiffs %>% dplyr::select(-Group), alphadiffs_rare, by = "PatientID")

richness_comparison_plot <- ggplot(alphadiffs_merged, aes(x = diversity_shannon.x, y = diversity_shannon.y)) +
  geom_point() +
  theme_classic() +
  labs(y= "Change Richness Rarefied Counts [%]", x = "Change Richness Raw Counts [%]")+
  xlim(-40, 90) + 
  ylim(-40, 90) +
  scale_x_continuous(breaks = seq(-40, 100, by = 20)) +
  scale_y_continuous(breaks = seq(-40, 100, by = 20)) +
  #scale_color_manual(values = c("#117733", "#661100", "#332288")) +
  geom_smooth(method = "lm", se = FALSE)
  

#ggsave(here("figures", "richness_comparison.pdf"), richness_comparison_plot)
#--------------

mapped_reads_stats%>%
  rename("Tube_ID"="Sample_ID")%>%
  dplyr::left_join(rare.alpha.div, by= "Tube_ID")%>%
  dplyr::filter(!is.na(observed))-> rare.alpha.div

mapped_reads_stats%>%
  rename("Tube_ID"="Sample_ID")%>%
  dplyr::left_join(motus.alpha.div, by= "Tube_ID")%>%
  dplyr::filter(!is.na(observed))-> motus.alpha.div

##Lets plot diversity taxa by DNA concentration
motus.alpha.div%>%
  ggplot(aes(x=Tax_mapped, y=diversity_shannon))+
  geom_point(aes(fill= Concentration_ng_ul), shape= 21,
             color= "black", size= 3,
             position = position_dodge(0.8))+
  scale_fill_gradient(low="#aec8ce", high="#9A342C", name= "DNA\nConcentration\n(ng/uL)")+
  ylab("Diversity (Shannon Index)")+
  xlab("Read count")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))-> A

##Lets plot beta diversity taxa-args by sample type 
#Bray-Curtis dissimilarity 
### 1) Taxa
bray_dist<- phyloseq::distance(PS.mOTUs, 
                               method="bray", weighted=T)
ordination<- ordinate(PS.mOTUs,
                      method="PCoA", distance="bray")

## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, sample_data$Sample_type, type = "centroid")

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Sample_type","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
motus.alpha.div%>%
  column_to_rownames("Tube_ID")%>%
  dplyr::select(!c(Sample_type))%>%
  cbind(seg.data)-> seg.data

##Adjust it to merge and make a single PCoA for both datasets 
seg.data%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::mutate(data_set= "mOTUS")-> seg.data.taxa

##Just to have an overview
ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Concentration_ng_ul), 
             shape= 21, size=3) +
  scale_fill_gradient(low="#aec8ce", high="#9A342C", name= "DNA\nConcentration\n(ng/uL)")+
  guides(color= "none")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> B

##Group them together 
Fig.0 <- grid.arrange(A,B,
                      ncol = 2,
                      nrow = 2,
                      layout_matrix = rbind(c(1,1),
                                            c(2,2))
)

##Save the map
## Save final plots
ggsave(filename = paste0(repo.path, "figures/Figure_0_2.pdf"), 
       plot = Fig.0,
       width = 15, 
       height = 15,
       units = "cm",
       dpi = 600)

ggsave(filename = paste0(repo.path, "figures/Figure_0_2.png"), 
       plot = Fig.0,
       width = 15, 
       height = 15,
       units = "cm",
       dpi = 600)

##Lets plot diversity taxa by groups
motus.alpha.div%>%
  dplyr::filter(!is.na(Group))%>%
  wilcox_test(diversity_shannon ~ Group)%>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

motus.alpha.div%>%
  dplyr::filter(!is.na(Group))%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot(aes(x=Group, y=diversity_shannon))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, shape=21, aes(fill= Group), color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0.1, step.increase = 0.05, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")-> A

##Lets plot beta diversity taxa-args by sample type 
#Bray-Curtis dissimilarity 
### 1) Taxa
PS.mOTUs<- subset_samples(PS.mOTUs, !is.na(Group))

bray_dist<- phyloseq::distance(PS.mOTUs, 
                               method="bray", weighted=T)
ordination<- ordinate(PS.mOTUs,
                      method="PCoA", distance="bray")

sample_data%>%
  dplyr::filter(!is.na(Group))-> tmp

##Set permutations
perm <- permute::how(nperm = 999)
##Define strata 
permute::setBlocks(perm) <- with(tmp, PatientID)

##Quantitative distance
groups.adonis<- vegan::adonis2(bray_dist~Group + Visit, 
                               permutations = perm, data = tmp, 
                               na.action = na.fail, by="margin")


## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, PS.mOTUs@sam_data$Group, type = "centroid")

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
motus.alpha.div%>%
  dplyr::filter(!is.na(Group))%>%
  column_to_rownames("Tube_ID")%>%
  dplyr::select(!c(Group))%>%
  cbind(seg.data)-> seg.data

##Just to have an overview
ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Group, shape= Visit), size=3) +
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), 
                    labels = c("Control", "Ketogenic", "Fasting"))+
  scale_shape_manual(values = c(21, 23),
                     labels= c("V0", "V3"))+
  guides(color= "none", fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Group), linetype = 2)+
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps , group=grps ), size=4, shape= 4) + 
  scale_color_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> B

##Group them together 
Fig.00 <- grid.arrange(A,B,
                       ncol = 2,
                       nrow = 2,
                       layout_matrix = rbind(c(1,1),
                                             c(2,2))
)

##Save the map
## Save final plots
ggsave(filename = paste0(repo.path, "figures/Figure_00_2.pdf"), 
       plot = Fig.00,
       width = 15, 
       height = 15,
       units = "cm",
       dpi = 600)

ggsave(filename = paste0(repo.path, "figures/Figure_00_2.png"), 
       plot = Fig.00,
       width = 15, 
       height = 15,
       units = "cm",
       dpi = 600)

##Lets plot diversity taxa by groups with rarefied data 
rare.alpha.div%>%
  wilcox_test(diversity_shannon ~ Group)%>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

rare.alpha.div%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot(aes(x=Group, y=diversity_shannon))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, shape=21, aes(fill= Group), color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0.1, step.increase = 0.05, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")-> A

rare.alpha.div%>%
  group_by(Visit)%>%
  wilcox_test(diversity_shannon ~ Group)%>%
  add_significance()%>%
  add_xy_position(x = "Visit")-> stats.test

rare.alpha.div%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  group_by(Group, Visit)%>%
  ggplot()+
  geom_boxplot(aes(x=Visit, y=diversity_shannon, fill=Group), alpha= 0.5, outlier.shape=NA,
               position = position_dodge(width = 0.75))+
  geom_point(aes(x=Visit, y=diversity_shannon, fill=Group),
             position=position_dodge(width = 0.75), size=3, shape=21, color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0.1, step.increase = 0.05, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")-> A

rare.alpha.div%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  group_by(PatientID) %>%
  filter(all(c("V0", "V3") %in% Visit)) %>%
  ungroup()%>%
  dplyr::group_by(PatientID)%>%
  dplyr::mutate(shannon_difference = diversity_shannon[Visit == "V3"] - diversity_shannon[Visit == "V0"])%>%
  ungroup()%>%
  wilcox_test(shannon_difference ~ Group)%>%
  add_significance()

rare.alpha.div%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  group_by(PatientID) %>%
  filter(all(c("V0", "V3") %in% Visit)) %>%
  ungroup()%>%
  dplyr::group_by(PatientID)%>%
  dplyr::mutate(shannon_difference = diversity_shannon[Visit == "V3"] - diversity_shannon[Visit == "V0"])%>%
  ungroup()%>%
  ggplot(aes(x=Group, y=shannon_difference))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, shape=21, aes(fill= Group), color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("Difference in diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())-> temp.fig.A

ggsave(filename = paste0(repo.path, "figures/Sup_Figure_0_2.pdf"), 
       plot = temp.fig.A,
       width = 20, 
       height = 20,
       units = "cm",
       dpi = 600)

rare.mOTU.counts<- as.data.frame(PS.Rare@otu_table)
rare.mOTU.counts%>%
  rownames_to_column("mOTUs_ID")-> rare.mOTU.counts

rare.mOTU.taxonomy<- as.data.frame(PS.Rare@tax_table)
rare.mOTU.taxonomy%>%
  rownames_to_column("mOTUs_ID")-> rare.mOTU.taxonomy

##DNA concentration
rare.alpha.div%>%
  wilcox_test(Concentration_ng_ul ~ Group)%>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

rare.alpha.div%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot(aes(x=Group, y=Concentration_ng_ul))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, shape=21, aes(fill= Group), color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("DNA concentration (ng/µL)")+
  labs(tag= "B)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())-> B

##Read counts before rarefy
rare.alpha.div%>%
  wilcox_test(Tax_mapped ~ Group)%>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

rare.alpha.div%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot(aes(x=Group, y=Tax_mapped))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, shape=21, aes(fill= Group), color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  scale_y_log10("Raw read count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  labs(tag= "C)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())-> C

##Lets plot beta diversity taxa by sample type 
#Bray-Curtis dissimilarity 
### 1) Taxa
bray_dist<- phyloseq::distance(PS.Rare, 
                               method="bray", weighted=T)
ordination<- ordinate(PS.Rare,
                      method="PCoA", distance="bray")

sample_data%>%
  dplyr::filter(!is.na(Group))-> tmp

##Set permutations
perm <- permute::how(nperm = 999)
##Define strata 
permute::setBlocks(perm) <- with(tmp, PatientID)

##Quantitative distance
groups.adonis<- vegan::adonis2(bray_dist~Group + Visit, 
                               permutations = perm, data = tmp, 
                               na.action = na.fail, by="margin")

## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, PS.Rare@sam_data$Group, type = "centroid")

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
rare.alpha.div%>%
  column_to_rownames("Tube_ID")%>%
  dplyr::select(!c(Group))%>%
  cbind(seg.data)-> seg.data

##Just to have an overview
ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Group, shape= Visit), size=3) +
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), 
                    labels = c("Control", "Ketogenic", "Fasting"))+
  scale_shape_manual(values = c(21, 23),
                     labels= c("V0", "V3"))+
  guides(color= "none", fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Group), linetype = 2)+
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps , group=grps ), size=4, shape= 4) + 
  scale_color_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> D

##Group them together 
Fig.000 <- grid.arrange(A,B,C,D,
                        ncol = 3,
                        nrow = 4,
                        layout_matrix = rbind(c(1,2,3),
                                              c(1,2,3),
                                              c(4,4,4),
                                              c(4,4,4))
)

##Save the map
## Save final plots
ggsave(filename = paste0(repo.path, "figures/Figure_000_2.pdf"), 
       plot = Fig.000,
       width = 30, 
       height = 20,
       units = "cm",
       dpi = 600)

ggsave(filename = paste0(repo.path, "figures/Figure_000_2.png"), 
       plot = Fig.000,
       width = 30, 
       height = 20,
       units = "cm",
       dpi = 600)

#For Friederike
write_csv(rare.alpha.div, paste0(repo.path, "data/Rarefied_alpha_diversity_B1.csv"))
write_csv(rare.mOTU.counts, paste0(repo.path, "data/Rarefied_mOTU_counts_B1.csv"))
write_csv(rare.mOTU.taxonomy, paste0(repo.path, "data/Rarefied_mOTU_taxonomy_B1.csv"))

##Get this for batch 2
PS.Rare.2<- prune_samples(!is.na(PS.mOTUs2@sam_data$Sample_Number), PS.mOTUs2)
PS.Rare.2<- rarefy_even_depth(PS.Rare.2, rngseed=2020, sample.size=4000, replace=F)

##Rarefied alpha diversity
rare.alpha.div.2 <- microbiome::alpha(PS.Rare.2, index =c("Observed","Chao1", "Shannon", "Simpson")) 
rare.alpha.div.2<- as.data.frame(rare.alpha.div.2)
rare.alpha.div.2%>%
  rownames_to_column("Tube_ID")%>%
  dplyr::left_join(sample_data_2, by= "Tube_ID")-> rare.alpha.div.2

mapped_reads_stats_2%>%
  rename("Tube_ID"="Sample_ID")%>%
  dplyr::left_join(rare.alpha.div.2, by= "Tube_ID")%>%
  dplyr::filter(!is.na(observed))-> rare.alpha.div.2

rare.mOTU.counts.2<- as.data.frame(PS.Rare.2@otu_table)
rare.mOTU.counts.2%>%
  rownames_to_column("mOTUs_ID")-> rare.mOTU.counts.2

rare.mOTU.taxonomy.2<- as.data.frame(PS.Rare.2@tax_table)
rare.mOTU.taxonomy.2%>%
  rownames_to_column("mOTUs_ID")-> rare.mOTU.taxonomy.2

#For Friederike
write_csv(rare.alpha.div.2, paste0(repo.path, "data/Rarefied_alpha_diversity_B2.csv"))
write_csv(rare.mOTU.counts.2, paste0(repo.path, "data/Rarefied_mOTU_counts_B2.csv"))
write_csv(rare.mOTU.taxonomy.2, paste0(repo.path, "data/Rarefied_mOTU_taxonomy_B2.csv"))

##Compare samples in both sequencing batches
rare.alpha.div%>%
  dplyr::mutate(Batch= "Batch_1")-> rare.alpha.div

rare.alpha.div.2%>%
  dplyr::mutate(Batch= "Batch_2")%>%
  rbind(rare.alpha.div)-> rare.alpha.div.all

rare.alpha.div.all%>%
  dplyr::filter(Sample_Number%in%c(22,28,30,52,103,120,135))%>%
  wilcox_test(Tax_mapped ~ Batch, paired = T)%>%
  add_significance()-> stats.test

rare.alpha.div.all%>%
  dplyr::filter(Sample_Number%in%c(22,28,30,52,103,120,135))%>%
  wilcox_test(observed ~ Batch, paired = T)%>%
  add_significance()%>%
  rbind(stats.test)-> stats.test

rare.alpha.div.all%>%
  dplyr::filter(Sample_Number%in%c(22,28,30,52,103,120,135))%>%
  wilcox_test(diversity_shannon ~ Batch, paired = T)%>%
  add_significance()%>%
  rbind(stats.test)-> stats.test

rare.alpha.div.all%>%
  dplyr::filter(Sample_Number%in%c(22,28,30,52,103,120,135))%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot()+
  geom_boxplot(aes(x=Batch, y=diversity_shannon), alpha= 0.5, outlier.shape=NA)+
  geom_line(aes(x = as.numeric(as.factor(Batch)) + (as.numeric(as.factor(Group)) - 2) * 0.25, 
                y = diversity_shannon, group = Sample_Number), 
            color = "gray", alpha = 0.5)+
  geom_point(aes(x=Batch, y=diversity_shannon, fill=Group),
             position=position_dodge(width = 0.75), size=3, shape=21, color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())-> A

rare.alpha.div.all%>%
  dplyr::filter(Sample_Number%in%c(22,28,30,52,103,120,135))%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot()+
  geom_boxplot(aes(x=Batch, y=Concentration_ng_ul), alpha= 0.5, outlier.shape=NA)+
  geom_line(aes(x = as.numeric(as.factor(Batch)) + (as.numeric(as.factor(Group)) - 2) * 0.25, 
                y = Concentration_ng_ul, group = Sample_Number), 
            color = "gray", alpha = 0.5)+
  geom_point(aes(x=Batch, y=Concentration_ng_ul, fill=Group),
             position=position_dodge(width = 0.75), size=3, shape=21, color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  ylab("DNA concentration (ng/µL)")+
  labs(tag= "B)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())-> B

rare.alpha.div.all%>%
  dplyr::filter(Sample_Number%in%c(22,28,30,52,103,120,135))%>%
  dplyr::mutate(Group = as.character(Group))%>%
  dplyr::mutate(Group = fct_relevel(Group, 
                                    "1", "2", "3"))%>%
  ggplot()+
  geom_boxplot(aes(x=Batch, y=Tax_mapped), alpha= 0.5, outlier.shape=NA)+
  geom_line(aes(x = as.numeric(as.factor(Batch)) + (as.numeric(as.factor(Group)) - 2) * 0.25, 
                y = Tax_mapped, group = Sample_Number), 
            color = "gray", alpha = 0.5)+
  geom_point(aes(x=Batch, y=Tax_mapped, fill=Group),
             position=position_dodge(width = 0.75), size=3, shape=21, color= "black")+
  scale_fill_manual(values = c("#70b5db", "#86ac55", "#ea9b5e"), labels = c("Control", "Ketogenic", "Fasting"))+
  scale_y_log10("Raw read count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  labs(tag= "C)", fill= "Group")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())-> C

##Group them together 
Fig.0000 <- grid.arrange(A,B,C,
                        ncol = 3,
                        nrow = 1,
                        layout_matrix = rbind(c(1,2,3))
)

##Save the map
## Save final plots
ggsave(filename = paste0(repo.path, "figures/Figure_0000.pdf"), 
       plot = Fig.0000,
       width = 30, 
       height = 10,
       units = "cm",
       dpi = 600)

ggsave(filename = paste0(repo.path, "figures/Figure_0000.png"), 
       plot = Fig.0000,
       width = 30, 
       height = 10,
       units = "cm",
       dpi = 600)

#Bray-Curtis dissimilarity for all the second Batch
### 1) Taxa
bray_dist<- phyloseq::distance(PS.Rare.2, 
                               method="bray", weighted=T)
ordination<- ordinate(PS.Rare.2,
                      method="PCoA", distance="bray")

##Set permutations
perm <- permute::how(nperm = 999)

##Quantitative distance
groups.adonis<- vegan::adonis2(bray_dist~ Visit, 
                               permutations = perm, data = rare.alpha.div.2, 
                               na.action = na.fail, by="margin")

## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, PS.Rare.2@sam_data$Visit, type = "centroid")

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Visit","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
rare.alpha.div.2%>%
  column_to_rownames("Tube_ID")%>%
  dplyr::select(!c(Visit))%>%
  cbind(seg.data)-> seg.data

##Just to have an overview
ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Visit), size=3, shape=21) +
  scale_fill_manual(values = c("#FFD9FF", "#DD5c55"), 
                    labels = c("V0", "V3"))+
  guides(color= "none", fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Visit), linetype = 2)+
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps , group=grps ), size=4, shape= 4) + 
  scale_color_manual(values = c("#FFD9FF", "#DD5c55"))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> D

ggsave(filename = paste0(repo.path, "figures/Figure_0000_2.pdf"), 
       plot = D,
       width = 20, 
       height = 20,
       units = "cm",
       dpi = 600)

ggsave(filename = paste0(repo.path, "figures/Figure_0000_2.png"), 
       plot = D,
       width = 20, 
       height = 20,
       units = "cm",
       dpi = 600)