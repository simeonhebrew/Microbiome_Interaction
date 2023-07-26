library(phyloseq)
library(vegan)
library(ggplot2)
library(plotly.microbiome)
library(DECIPHER)
library(phangorn)
library(microeco)


phyloseq_na_cluster <- readRDS("/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/phyloseq_na_cluster.rds")
phyloseq_na_cluster_its <- readRDS("/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/phyloseq_na_cluster_its.rds")

phyloseq_ms_pos <- readRDS("/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/phyloseq_new_set.rds")


sample_ms <- sample_data(phyloseq_ms_pos)
otu_ms <- otu_table(phyloseq_ms_pos)
otu_ms <- as.data.frame(otu_ms)
tax_table(phyloseq_ms_pos)



write.csv(sample_ms, "/Users/cimi_bioinformatics/Desktop/sample_ms.csv")


phyloseq_na_cluster_its
phyloseq_na_cluster

#Pre_Phyloseq_Object



pre_samples <- sample_data(phyloseq_pre_samples)
fungi_samples <- sample_data(phyloseq_na_cluster_its)
pre_sel <- sample_data(pre_selected_samples)



#Age_gender_file

age_gender <- read.csv("/Users/cimi_bioinformatics/Downloads/Age_Gender_Samples - Sheet1.csv", header = FALSE)


pre_select <- pre_samples %>% 
  filter(str_detect(age_gender$V1, "CER*"))




#Filter samples

na_cluster_filt = prune_samples(sample_sums(phyloseq_na_cluster)>=2500, phyloseq_na_cluster)
na_cluster_filt

#Pre_Phyloseq


phyloseq_pre_samples <- subset_samples(na_cluster_filt,
                                       Study_inclusion == "CTRL_CIS_MS_pre")
phyloseq_pre_samples

sequences <- rownames(seqtab_curated_org)

alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA, processors = 16)

saveRDS(alignment, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/alignment.rds")
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
#assigning distance matrix
dm <- dist.ml(phang.align)
#creating the phylogenetic tree
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0), multithread = TRUE)
saveRDS(fitGTR,"phylogenetic_tree.rds")
#rooting a the tree
set.seed(711)
phy_tree(phyloseq_object) <- root(phy_tree(phyloseq_object), sample(taxa_names(phyloseq_object), 1), resolve.root = TRUE)
is.rooted(phy_tree(phyloseq_object))










#Pos_Phyloseq

phyloseq_pos_samples <- subset_samples(na_cluster_filt,
                                       Study_inclusion %in% c("POS", "POS_100k"))
phyloseq_pos_samples


pre_selected_samples <- subset_samples(phyloseq_pre_samples,
                                       Subject %in% c("CER069", 
                                                      "CER119",   
                                                      "CER202", 
                                                      "CER203",
                                                      "CER206",
                                                      "CER207",  
                                                      "CER210",
                                                      "CER211",  
                                                      "CER112",  
                                                      "CER213",
                                                      "CER214",  
                                                      "CER218",  
                                                      "CER220",
                                                      "CER223",  
                                                      "CER224",  
                                                      "CER226",  
                                                      "ENT001",  
                                                      "ENT002",  
                                                      "ENT003",  
                                                      "ENT004",  
                                                      "ENT005",  
                                                      "ENT006",  
                                                      "ENT007",  
                                                      "ENT010",  
                                                      "ENT012",  
                                                      "ENT014",  
                                                      "ENT016", 
                                                      "ENT020",  
                                                      "ENT023", 
                                                      "ENT024",
                                                      "ENT025", 
                                                      "ENT028",  
                                                      "ENT030",  
                                                      "ENT033",  
                                                      "ENT034",  
                                                      "ENT035",  
                                                      "ENT036",  
                                                      "ENT049",  
                                                      "ENT050",  
                                                      "ENT053",  
                                                      "ENT061", 
                                                      "ENT066",  
                                                      "ENT068",  
                                                      "ENT083",  
                                                      "ENT101",  
                                                      "ENT103", 
                                                      "ENT105", 
                                                      "ENT106", 
                                                      "ENT116",  
                                                      "ENT117"))



fungi_selected_samples <- subset_samples(phyloseq_na_cluster_its,
                                         SampleID %in% c("CER069", 
                                                        "CER119",   
                                                        "CER202", 
                                                        "CER203",
                                                        "CER206",
                                                        "CER207",  
                                                        "CER210",
                                                        "CER211",  
                                                        "CER212",  
                                                        "CER213",
                                                        "CER214",  
                                                        "CER218",  
                                                        "CER220",
                                                        "CER223",  
                                                        "CER224",  
                                                        "CER226",  
                                                        "ENT001",  
                                                        "ENT002",  
                                                        "ENT003",  
                                                        "ENT004",  
                                                        "ENT005",  
                                                        "ENT006",  
                                                        "ENT007",  
                                                        "ENT010",  
                                                        "ENT012",  
                                                        "ENT014",  
                                                        "ENT016", 
                                                        "ENT020",  
                                                        "ENT023", 
                                                        "ENT024",
                                                        "ENT025", 
                                                        "ENT028",  
                                                        "ENT030",  
                                                        "ENT033",  
                                                        "ENT034",  
                                                        "ENT035",  
                                                        "ENT036",  
                                                        "ENT049",  
                                                        "ENT050",  
                                                        "ENT053",  
                                                        "ENT061", 
                                                        "ENT066",  
                                                        "ENT068",  
                                                        "ENT083",  
                                                        "ENT101",  
                                                        "ENT103", 
                                                        "ENT105", 
                                                        "ENT106", 
                                                        "ENT116",  
                                                        "ENT117"))






pre_df <- tabcheck %>% 
  filter(str_detect(rownames(tabcheck), "_pre"))

pos_df <- tabcheck %>% 
  filter(str_detect(rownames(tabcheck), "_pos*"))

tab <- t(tab)

tax_table <- as.vector(tax_table)
tax_table <- t(tax_table)

rare_pre <- rarecurve(tax_table, step=100, lwd=2, ylab="Species_Richness",  label=T)
rare_pre



tax_table <- otu_table(phyloseq_na_250_phy)
class(tax_table)


rarecurve(tax_table, step=50, cex=0.5)

rare_pos <- rarecurve(pos_df, step=100, lwd=2, ylab="Species_Richness",  label=T)
rare_pos

ps.rarefied = rarefy_even_depth(na_cluster_filt, rngseed=1, sample.size=0.9*min(sample_sums(na_cluster_filt)), replace=F)
ps.rarefied_fungi = rarefy_even_depth(fungi_selected_samples, rngseed=1, sample.size=0.9*min(sample_sums(na_cluster_filt)), replace=F)

#Rarecurve
rarecurve(t(otu_table(tab)), step=5, cex=0.5)

f <- otu_table(na_cluster_filt)

f <- as(otu_table(na_cluster_filt), "matrix")






alpha_pre <- subset_samples(pre_selected_samples,
                            Condition_disease %in% c("CTRL", "CIS", "MS"))
                    
alpha_pre <- as.factor(alpha_pre)

#Alpha Diversity

comparisons = list( c("CTRL", "RR"), c("CTRL", "SP"), c("CTRL", "CIS"))

alpha_bacteria_cluster = plot_richness(mphlan_control, x="Type", color="Type", measures=c("Observed")) + 
  geom_boxplot(alpha=0.5, lwd=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() , axis.line = element_line(colour = "black")) +
  scale_x_discrete(limits = c("Pos", "Neg"))
  
alpha_bacteria_cluster

bac_results = estimate_richness(pre_selected_samples, measures = 'Shannon')
d = sample_data(pre_selected_samples)

# calculate t-test
CTRL = bac_results[d[,'Condition_disease'] == 'CTRL',]
MS = bac_results[d[,'Condition_disease'] == 'MS',]
wilcox.test(CTRL, MS) # 0.1993


CTRL = bac_results[d[,'Condition_disease'] == 'CTRL',]
CIS = bac_results[d[,'Condition_disease'] == 'CIS',]
wilcox.test(CTRL, CIS) #0.4706


alpha_fungi_cluster = plot_richness(fungi_selected_samples, x="DiseaseEvolution", color="DiseaseEvolution", measures=c("Observed", "Chao1", "Shannon", "Simpson" )) + 
  geom_boxplot(alpha=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() , axis.line = element_line(colour = "black")) +
  geom_signif(comparisons = comparisons, 
              map_signif_level=TRUE)

alpha_fungi_cluster


#Beta Diversity

otu_table(pre_selected_samples)[is.na(otu_table(pre_selected_samples))] <- 0
otu_table(pre_selected_samples)[otu_table(pre_selected_samples)==0]<-0.1

ordination_bac_pre_cluster = ordinate(bac_sel_dataset_physeq, method="PCoA", distance="bray")

class(ordination_bac_pre_cluster)

ordinate_bray <- as.matrix(ordination_bac_pre_cluster)

ps_bray <- phyloseq::distance(phyloseq_250_sub, method = "bray")


sample_df_250 <- as.data.frame(sample_data(ordination_bac_pre_cluster))


sam <- (sample_data(phyloseq_250_sub))

meta_colors <- c(CIS = "#D95F02", CTRL = "#1B9E77", RR = "#7570B3",
                   SP = "#E7298A")

meta_colors <- c(CTRL = "#1B9E77",CIS = "#D95F02" , MS = "#7570B3")



beta_diversity_3d(ordination_bac_pre_cluster, 
                  sample_data(bac_sel_dataset_physeq), "Condition",
                  color.key = meta_colors)


adonis2(ps_bray)


metadata <- as(sample_data(bac_sel_dataset_physeq), "data.frame")
dist.bc <- phyloseq::distance(bac_sel_dataset_physeq, method = "bray")
dist.bc <- t(as.matrix(dist.bc))
permanova <- adonis(dist ~ Condition, data = metadata, perm=9999)


dist <- unname(dist.bc)
permanova$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Condition  2    0.4050 0.20253  1.3787 0.06161 0.0073 **
#  Residuals 42    6.1697 0.14690         0.93839          
#Total     44    6.5747                 1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


metadata_its <- as(sample_data(fun_sel_dataset_physeq), "data.frame")
dist.bc_its <- phyloseq::distance(fun_sel_dataset_physeq, method = "bray")
dist_its <- unname(dist.bc_its)
permanova_its <- adonis(dist_its ~ Condition, data = metadata_its, perm=9999)

permanova_its$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Condition  2    0.6029 0.30147  2.3076 0.11363 0.0077 **
#  Residuals 36    4.7031 0.13064         0.88637          
#Total     38    5.3060                 1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Pos_ControlCIS

metadata_pos <- as(sample_data(mod_pos_edit), "data.frame")
dist.bc_pos <- phyloseq::distance(mod_pos_edit, method = "bray")
dist_pos <- unname(dist.bc_pos)
permanova_pos <- adonis(dist_pos ~ Condition, data = metadata_pos, perm=9999)

permanova_pos$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Condition  1    0.1622 0.16219  1.2942 0.04922 0.1126
#Residuals 25    3.1330 0.12532         0.95078       
#Total     26    3.2951                 1.00000  


##Neg_ControlCIS
metadata_neg <- as(sample_data(mod_negsamples), "data.frame")
dist.bc_neg <- phyloseq::distance(mod_negsamples, method = "bray")
dist_neg <- unname(dist.bc_neg)
permanova_neg <- adonis(dist_neg ~ Condition, data = metadata_neg, perm=9999)

permanova_neg$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Condition  1    0.1160 0.11598 0.83425 0.03229 0.9612
#Residuals 25    3.4756 0.13902         0.96771       
#Total     26    3.5915                 1.00000



##ShotgunPos_ControlCIS
metadata_poshot <- as(sample_data(mphlan_pos), "data.frame")
dist.bc_poshot <- phyloseq::distance(mphlan_pos, method = "bray")
dist_poshot <- unname(dist.bc_poshot)
permanova_poshot<- adonis(dist_poshot ~ Condition, data = metadata_poshot, perm=9999)

permanova_poshot$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Condition  1    0.3566 0.35658  1.1726 0.04019 0.2232
#Residuals 28    8.5149 0.30410         0.95981       
#Total     29    8.8714                 1.00000




#PosNeg_ControlCIS

metadata_posnegcondition <- as(sample_data(mod_control_cis), "data.frame")
dist.bc_posnegcondition <- phyloseq::distance(mod_control_cis, method = "bray")
permanova_condition <- adonis(dist.bc_posnegcondition ~ Condition, data = metadata_posnegcondition, perm=9999)

permanova_condition$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Condition  1    0.1484 0.14843  1.1323 0.02092 0.1855
#Residuals 53    6.9480 0.13109         0.97908       
#Total     54    7.0964                 1.00000 


#PosNeg_Control

metadata_posnegcontrol <- as(sample_data(mod_control), "data.frame")
dist.bc_posnegcontrol <- phyloseq::distance(mod_control, method = "bray")
dist.bc_posnegcontrol <- unname(dist.bc_posnegcontrol)
permanova_control <- adonis(dist.bc_posnegcontrol ~ Study, data = metadata_posnegcontrol, perm=9999)

permanova_control$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Study      1    0.2941 0.29411  2.3237 0.07924  1e-04 ***
#Residuals 27    3.4174 0.12657         0.92076           
#Total     28    3.7115                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#PosNeg_CIS
metadata_posnegcis <- as(sample_data(mod_cis), "data.frame")
dist.bc_posnegcis <- phyloseq::distance(mod_cis, method = "bray")
dist.bc_posnegcis <- unname(dist.bc_posnegcis)
permanova_cis <- adonis(dist.bc_posnegcis ~ Study, data = metadata_posnegcis, perm=9999)

permanova_cis$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
#Study      1    0.2644 0.26442  2.1352 0.0817  2e-04 ***
#Residuals 24    2.9721 0.12384         0.9183           
#Total     25    3.2366                 1.0000           
#---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Shotgun_CTRLCIS

metadata_mphlancontrolcis <- as(sample_data(mphlan_dataset_physeq), "data.frame")
dist.bc_mphlancontrolcis <- phyloseq::distance(mphlan_dataset_physeq, method = "bray")
permanova_mphlancontrolcis <- adonis(dist.bc_mphlancontrolcis ~ Condition, data = metadata_mphlancontrolcis, perm=9999)

permanova_mphlancontrolcis$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Condition  1    0.4203 0.42032   1.395 0.02349 0.0799 .
#Residuals 58   17.4758 0.30131         0.97651         
#Total     59   17.8961                 1.00000         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Shotgun_CTRL_PosNeg

metadata_mphlancontrol <- as(sample_data(mphlan_control), "data.frame")
dist.bc_mphlancontrol <- phyloseq::distance(mphlan_control, method = "bray")
dist_ctrlposneg <- unname(dist.bc_mphlancontrol)
permanova_mphlancontrol <- adonis(dist_ctrlposneg ~ Type, data = metadata_mphlancontrol, perm=9999)


permanova_mphlancontrol$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Type       1    0.8529 0.85289  3.2545 0.09787  1e-04 ***
#Residuals 30    7.8619 0.26206         0.90213           
#Total     31    8.7148                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Shotgun_CIS_PosNeg

metadata_mphlancis <- as(sample_data(mphlan_cis), "data.frame")
dist.bc_mphlancis <- phyloseq::distance(mphlan_cis, method = "bray")
permanova_mphlancis <- adonis(dist.bc_mphlancis ~ Type, data = metadata_mphlancis, perm=9999)


permanova_mphlancis$aov.tab

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
#Type       1    1.1468 1.14679  3.9159 0.1309  1e-04 ***
#Residuals 26    7.6142 0.29286         0.8691           
#Total     27    8.7610                 1.0000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1






Y <- samp[, c("DCOX", "SPU", "PUM")]
## Perform a one-way NPMANOVA:
adonis(Y ~ samp$Sex, method = "euclidean",
       permutations = 999)


#Fungi
otu_table(fungi_selected_samples)[is.na(otu_table(fungi_selected_samples))] <- 0
otu_table(fungi_selected_samples)[otu_table(fungi_selected_samples)==0]<-0.1

ordination_bac_pre_cluster_fungi = ordinate(fun_sel_dataset_physeq, method="PCoA", distance="bray")



meta_colors <- RColorBrewer::brewer.pal(7, "Dark2")



beta_diversity_3d(ordination_bac_pre_cluster_fungi, 
                  sample_data(fun_sel_dataset_physeq), "Condition",
                  color.key = meta_colors)


brewer.pal(n = 7, name = "Dark2")


# NMDS with UniFrac distance (ordinate won't let you adjust k)

phy_tree(pre_selected_samples, errorIfNULL=TRUE)


pre_unifrac <- phyloseq::distance(pre_selected_samples, "unifrac")
ord_mfiber_prop4a <- vegan::metaMDS(dist_mfiber_4, k = 3)

beta_diversity_3d(ord_mfiber_prop4a,
                  sample_data(ps_mfiber_prop), "TRT", "Label")


pr <- prevalence(pre_selected_samples, detection=0, sort=TRUE, count=TRUE)
pr <- as.data.frame(pr)



#Differential analysis


otu_table(phyloseq_250_sub)[otu_table(phyloseq_250_sub)==0]<-1

otu_table(phyloseq_250_sub_its)[otu_table(phyloseq_250_sub_its)==0]<-1
 

sample_data(phyloseq_250_sub)$Condition <- as.factor(sample_data(phyloseq_250_sub)$Condition)
sample_data(phyloseq_250_sub_its)$Condition <- as.factor(sample_data(phyloseq_250_sub_its)$Condition)




ds_bac = phyloseq_to_deseq2(phyloseq_250_sub, ~ Condition)
ds_fungi = phyloseq_to_deseq2(phyloseq_250_sub_its, ~ Condition)
ds_bac = DESeq(ds_bac)
ds_fungi = DESeq(ds_fungi)


reveal <- results(ds)

ds = phyloseq_to_deseq2(phyloseq_pos, ~ PCR_BY)
ds = DESeq(ds)

ds_bac = estimateSizeFactors(ds_bac, type="iterate")



alpha = 0.005
res = results(ds_bac, contrast=c("Condition", "CTRL", "RR"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig

alpha = 0.01
res = results(ds_fungi, contrast=c("DiseaseEvolution", "CTRL", "SP"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig







pos_log = -5
neg_log = 2


res = res[order(res$log2FoldChange, na.last=NA), ]
res_sig <- res[res$log2FoldChange >= alpha, ]
res_sigg <- res_sig[res_sig$log2FoldChange >= neg_log, ]
summary(res_sig)

res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(phyloseq_250_sub)[rownames(res_sig), ], "matrix"))

ggplot () + geom_boxplot(res_sig, mapping=aes(x=Species, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) #+
#theme(axis.text.x=element_blank())

ggplot(res_sig, aes(x=Species, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() 









##Prevalence

pseq.rel <- microbiome::transform(phylomerge, "compositional")
pseq.rel_fungi <- microbiome::transform(phylomerge_its, "compositional")


ps.m3.rel <- aggregate_taxa(pseq.rel, "Species")
ps.m3.rel_fun <- aggregate_taxa(pseq.rel_fungi, "Species")


pseq.rel_control <- microbiome::transform(mphlan_dataset_physeq, "compositional")
pseq.rel_cis <- microbiome::transform(mphlan_cis, "compositional")``
ps.m3.rel_control <- aggregate_taxa(pseq.rel_control, "Species")

ps.m3.rel_cis <- aggregate_taxa(mphlan_cis, "Species")




bac_ms <- subset_samples(bac_sel_dataset_physeq,
                            Condition == "MS")

fun_cis <- subset_samples(fun_sel_dataset_physeq,
                         Condition == "CIS")

fun_ms <- subset_samples(fun_sel_dataset_physeq,
                          Condition == "MS")

pseqbacsel <- microbiome::transform(bac_ms, "compositional")
pseqbacsel_rel <- aggregate_taxa(pseqbacsel, "Species")

pseqfunsel <- microbiome::transform(fun_ms, "compositional")
pseqfunsel_rel <- aggregate_taxa(pseqfunsel, "Species")

pseqpossel <- microbiome::transform(mphlan_pos, "compositional")
pseqpossel_rel <- aggregate_taxa(pseqpossel, "Species")




saveRDS(ps.m3.rel.sp, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/phyloseq_network_16.rds")
saveRDS(ps.m3.rel.sp_fun, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/phyloseq_network_ITS.rds")

taxa_names(ps.m3.rel.sp_fun)
taxa_names(ps.m3.rel.sp)

prev_mphlan <- prevalence(ps.m3.rel.batch_fun, detection = 90/100, sort = TRUE, prev = .5)
prev <- as.data.frame(prev)

phyloseq_filter_prevalence(
  ps.m3.rel.sp,
  prev.trh = 0.1,
  abund.trh = NULL,
  threshold_condition = "OR",
  abund.type = "total"
)
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(ps.m3.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()



pseq.core_df <- core(ps.m3.rel, detection = 0.0001, prevalence = 10/100)
pseq.core_df


pseq.core_fun <- core(pseqfunsel_rel, detection = 0.001, prevalence = 10/100)
pseq.core_fun


pseq.core_pos <- core(pseqpossel_rel, detection = 0.001, prevalence = 10/100)
pseq.core_pos

bac_sel_dataset_physeq


core_members(ps.m3.rel_control,  prevalence = 50/100)
core_members(ps.m3.rel.sp_fun,  prevalence = 1/100)

core.taxa.standard <- core_members(ps.m3.rel.sp, prevalence = 99/100)
core.taxa.standard 

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(ps.m3.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()






library(RColorBrewer)
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)
#detections <- 10^seq(log10(1e-4), log10(.2), length = 10), 3)


taxa_names(ps.m3.rel.sp) <- gsub("Firmicutes_Clostridia_Lachnospirales_Lachnospiraceae_Blautia_A_",
                                  "", taxa_names(ps.m3.rel.sp))

p1 <- plot_core(pseq.core, 
                plot.type = "heatmap", 
                colours = rev(brewer.pal(5, "RdBu")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = 0.7) +
  xlab("Detection Threshold (Relative Abundance (%))") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size = 10),
        axis.text.y = element_text(face = "italic", size = 15))
p1

p2 <- plot_core(pseq.core_fun, 
                plot.type = "heatmap", 
                colours = rev(brewer.pal(5, "RdBu")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .07) +
  xlab("Detection Threshold (Relative Abundance (%))") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.text.y = element_text(face = "italic", size = 17)) 
  
p2




#Comparison_Samples
comparison <- read.csv("/Users/cimi_bioinformatics/Downloads/Comparison_Samples - Sheet1.csv", sep = ",", header =TRUE)


z <- data.frame(
  X1 = comparison$Sample2[na.omit(match(comparison$Sample, comparison$Sample2))]
)
z



## new function to determine no. of distinct OTUs/genes/features occurrence per sample: 
# it returns a data frame and a ggplot2 plot
count_no_features <- function(physeq) {
  
  # physeq: a phyloseq-class object
  
  require("phyloseq")
  require("ggplot2")
  
  # check if taxa are rows; if not transpose 
  if ( taxa_are_rows(physeq) ) { #returns mtx
    otu_tbl <- as( otu_table(physeq), "matrix" )
  } else { # transpose; returns mtx
    otu_tbl <- t( otu_table(physeq) ) 
  }
  
  # sum no. of distinct OTUs observed by sample (different than 0)
  no_otus <- apply( otu_tbl, 2, function(x) sum(x != 0) )
  
  # save the result into a df
  no_otus_df <- data.frame("Samples" = names(no_otus),
                           "Type_alpha" = rep("Observed", length(no_otus)),
                           "Alpha_div" = no_otus, 
                           row.names = names(no_otus))
  colnames(no_otus_df) <- c("Samples", "Type of alpha-diversity", "Alpha-diversity measure") # change column names
  
  ## plot results
  plot_rich <- ggplot(data = no_otus_df, aes(x = Samples, y = `Alpha-diversity measure`)) + 
    geom_point() + 
    facet_wrap(~ `Type of alpha-diversity`) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() , axis.line = element_line(colour = "black"))
  
  
  ## return a list 
  list_results <- list("data" = no_otus_df, 
                       "plot" = plot_rich)
  
  return(list_results)
}


no_of_otus_new_fun <- count_no_features(physeq = phyloseq_mphlan_batch)


no_of_otus_new_fun$plot



df <- data.frame(Condition=c("Control", "CIS"),
                 Number_of_interactions=c(184, 150))
head(df)


mphlan_neg_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))

df <- df$Condition %<>% factor(., levels = c("CTRL", "CIS"))
df <-as.data.frame(df)


p<-ggplot(data=df, aes(x=Condition, y=Number_of_interactions, fill=Condition, show.legend = FALSE)) +
  geom_bar(stat="identity")  + scale_fill_brewer(palette="Dark2") +
  theme(axis.text.y=element_text(size=20)) + theme(axis.text.x=element_text(size=20)) + 
  theme(axis.title.y=element_text(size=20)) + theme(axis.title.x=element_text(size=20)) +
  scale_x_discrete(labels=c('Control','CIS')) + 
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black") + theme(legend.position = "none")
  )

p


#rownames(phycounts_its_cis)[rownames(phycounts_its_cis) == "ENT126"] <- "ENT049_2nd_E"
#rownames(phycounts_its_cis)[rownames(phycounts_its_cis) == "MUN001"] <- "ENT020_2nd_E"
#rownames(phycounts_its_cis)[rownames(phycounts_its_cis) == "MUN003"] <- "ENT020_1st_FT"
#rownames(phycounts_its_cis)[rownames(phycounts_its_cis) == "MUN022"] <- "ENT050_1st_FT"
#rownames(phycounts_its_cis)[rownames(phycounts_its_cis) == "ENT101"] <- "ENT116_2nd_E"
#rownames(phycounts_its_cis)[rownames(phycounts_its_cis) == "ENT103"] <- "ENT116_1st_FT_1_20"
