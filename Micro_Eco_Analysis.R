# If a package is not installed, it will be installed from CRAN.
# First select the packages of interest
packages <- c("MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "picante", "pheatmap", "igraph", "rgexf", 
              "ggalluvial", "ggh4x", "rcompanion", "FSA", "gridExtra", "aplot", "NST", "GGally")
# Now check or install
for(x in packages){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

#Let's create a microtable object with mor'e information


phyloseq_mphlan_batch <- readRDS("/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/phyloseq_mphlan_batch_modphylo.rds")


library(file2meco)
library(microeco)
library(magrittr)
pre_dataset_physeq <- phyloseq(pre_dataset$otu_table, pre_dataset$tax_table, pre_dataset$sample_table)
fungi_dataset_physeq <- phyloseq(fungi_pre_dataset$otu_table, fungi_pre_dataset$tax_table, fungi_pre_dataset$sample_table)
mphlan_dataset_physeq <- phyloseq(batch_pre_dataset$otu_table, batch_pre_dataset$tax_table, batch_pre_dataset$sample_table)

bac_sel_dataset_physeq <- phyloseq(bac_sel_dataset$otu_table, bac_sel_dataset$tax_table, bac_sel_dataset$sample_table)
fun_sel_dataset_physeq <- phyloseq(fun_sel_dataset$otu_table, fun_sel_dataset$tax_table, fun_sel_dataset$sample_table)

merge_new_dataset_physeq <- phyloseq(merge_new_set_dataset$otu_table, merge_new_set_dataset$tax_table, merge_new_set_dataset$sample_table)
merge_new_dataset_its_physeq <- phyloseq(merge_new_set_dataset_its$otu_table, merge_new_set_dataset_its$tax_table, merge_new_set_dataset_its$sample_table)

pos_dataset_physeq <- phyloseq(pos_dataset$otu_table, pos_dataset$tax_table, pos_dataset$sample_table)

new_dataset_physeq <- phyloseq(new_set_dataset$otu_table, new_set_dataset$tax_table, new_set_dataset$sample_table)

saveRDS(new_dataset_physeq, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/new_dataset_physeq.rds")


mphlan_control <- subset_samples(mphlan_dataset_physeq,
               Condition == "CTRL")
mphlan_cis <- subset_samples(mphlan_dataset_physeq,
                                 Condition == "CIS")

mphlan_pos <- subset_samples(mphlan_dataset_physeq,
                                 Type == "Pos")

mphlan_pos_ctrl <- subset_samples(mphlan_pos,
                             Condition == "CTRL")

mphlan_pos_cis <- subset_samples(mphlan_pos,
                                  Condition == "CIS")

mphlan_neg <- subset_samples(mphlan_dataset_physeq,
                             Type == "Neg")


mphlan_neg <- subset_samples(mphlan_dataset_physeq,
                             Type == "Neg")

mphlan_neg_ctrl <- subset_samples(mphlan_neg,
                                  Condition == "CTRL")
mphlan_neg_cis <- subset_samples(mphlan_neg,
                                  Condition == "CIS")



mod_pos_neg_ctrl_cis <- subset_samples(new_dataset_physeq,
                                      Condition %in% c("CTRL", "CIS"))


mod_possamples  <- subset_samples(new_dataset_physeq,
                              Study == "IgA_pos")

mod_possamples2  <- subset_samples(new_dataset_physeq,
                                  Study == "IgA_pos")


mod_possamples_ctrl <- subset_samples(mod_pos_edit,
                             Condition == "CTRL")

mod_possamples_cis <- subset_samples(mod_pos_edit,
                                      Condition == "CIS")

####
df["A", ] <- df["A", ] + df["C", ]
df[rownames(df) != "C", ]
####


sample_list <- c("CER211-pos500k", "CER233-pos-A4")

mod_possamples_pruned <- prune_samples(!(sample_names(mod_possamples) %in% sample_list), mod_possamples)

mod_pos_edit <- remove_samples(samples = sample_list, mod_possamples)
mod_pos_edited <- phyloseq2meco(mod_pos_edit)


saveRDS(mod_pos_edit, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/mod_pos_edit.rds")

mod_pos_pruned <- phyloseq2meco(mod_possamples_pruned)

mod_otu <- as.data.frame(otu_table(mod_possamples_pruned))

otu_table(mod_possamples)[,"CER211-pos"] <- otu_table(mod_possamples)[,"CER211-pos"] + otu_table(mod_possamples)[,"CER211-pos500k"]
otu_table(mod_possamples)[,"CER233-pos"] <- otu_table(mod_possamples)[,"CER233-pos"] + otu_table(mod_possamples)[,"CER233-pos-A4"]


edit_otu <- as.data.frame(otu_table(mod_pos_edit))


otu_table(mod_possamples) <- otu_table(mod_possamples)[colnames(otu_table(mod_possamples)) != "CER211-pos500k" ]
otu_table(mod_possamples)

otu_table(mod_possamples)


mod_negsamples  <- subset_samples(new_dataset_physeq,
                                  Study == "IgA_neg")

mod_negsamples_dataset <- phyloseq2meco(mod_negsamples)



mod_negsamples_ctrl  <- subset_samples(mod_negsamples,
                                  Condition == "CTRL")

mod_negsamples_cis  <- subset_samples(mod_negsamples,
                                  Condition == "CIS")


metadata_newdataset <- sample_data(new_dataset_physeq)



mod_control <- subset_samples(new_dataset_physeq,
                                 Condition == "CTRL")
mod_cis <- subset_samples(new_dataset_physeq,
                             Condition == "CIS")

mod_control_cis <- subset_samples(new_dataset_physeq,
                          Condition %in% c("CTRL", "CIS"))


bacctrl <- subset_samples(bacteria_250_physeq,
                                  Condition == "CTRL")

baccis <- subset_samples(bacteria_250_physeq,
                          Condition == "CIS")

bacrr <- subset_samples(bacteria_250_physeq,
                          Condition == "RR")

bacsp <- subset_samples(bacteria_250_physeq,
                          Condition == "SP")




bacteria_250_physeq <- phyloseq(bacteria_250_dataset$otu_table, bacteria_250_dataset$tax_table, bacteria_250_dataset$sample_table)

fungi_250_physeq <- phyloseq(fungi_250_dataset$otu_table, fungi_250_dataset$tax_table, fungi_250_dataset$sample_table)


pre_dataset <- phyloseq2meco(pre_dataset_physeq)
fungi_pre_dataset <- phyloseq2meco(fungi_dataset_physeq)
posseq_dataset <- phyloseq2meco(pos_dataset_physeq)
merge_new_dataset <- phyloseq2meco(merge_new_dataset_physeq)
merge_new_dataset_its <- phyloseq2meco(merge_new_dataset_its_physeq)
mphlan_pos_dataset <- phyloseq2meco(mphlan_pos)
mphlan_neg_dataset <- phyloseq2meco(mphlan_neg)


new_dataset_posneg <- phyloseq2meco(new_dataset_physeq)

mod_control_dataset <- phyloseq2meco(mod_control)
mod_cis_dataset <- phyloseq2meco(mod_cis)
mod_control_cis_dataset <- phyloseq2meco(mod_control_cis)


bacteria_250_eco <- phyloseq2meco(bacteria_250_physeq)
fungi_250_eco <- phyloseq2meco(fungi_250_physeq)

mphlan_control_dataset <- phyloseq2meco(mphlan_control)
mphlan_cis_dataset <- phyloseq2meco(mphlan_cis)



bacteria_selected_dataset <- phyloseq2meco(bac_sel_dataset_physeq)
fungi_selected_dataset <- phyloseq2meco(fun_sel_dataset_physeq)


bacteria_250_eco$cal_alphadiv(PD = FALSE)
bacteria_250_eco$cal_betadiv(unifrac = FALSE)
bacteria_250_eco$cal_abund()

bacteria_selected_dataset$cal_alphadiv(PD = FALSE)
bacteria_selected_dataset$cal_betadiv(unifrac = FALSE)
bacteria_selected_dataset$cal_abund()

new_dataset_posneg$cal_alphadiv(PD = FALSE)
new_dataset_posneg$cal_betadiv(unifrac = FALSE)
new_dataset_posneg$cal_abund()


fungi_selected_dataset$cal_alphadiv(PD = FALSE)
fungi_selected_dataset$cal_betadiv(unifrac = FALSE)
fungi_selected_dataset$cal_abund()



posseq_dataset$cal_alphadiv(PD = FALSE)
posseq_dataset$cal_betadiv(unifrac = FALSE)
posseq_dataset$cal_abund()

mod_control_dataset$cal_alphadiv(PD = FALSE)
mod_control_dataset$cal_betadiv(unifrac = FALSE)
mod_control_dataset$cal_abund()

mod_cis_dataset$cal_alphadiv(PD = FALSE)
mod_cis_dataset$cal_betadiv(unifrac = FALSE)
mod_cis_dataset$cal_abund()


mod_pos_pruned$cal_alphadiv(PD = FALSE)
mod_pos_pruned$cal_betadiv(unifrac = FALSE)
mod_pos_pruned$cal_abund()


mphlan_pos_dataset$cal_alphadiv(PD = FALSE)
mphlan_pos_dataset$cal_betadiv(unifrac = FALSE)
mphlan_pos_dataset$cal_abund()

mphlan_neg_dataset$cal_alphadiv(PD = FALSE)
mphlan_neg_dataset$cal_betadiv(unifrac = FALSE)
mphlan_neg_dataset$cal_abund()

mphlan_control_dataset$cal_alphadiv(PD = FALSE)
mphlan_control_dataset$cal_betadiv(unifrac = FALSE)
mphlan_control_dataset$cal_abund()

mphlan_cis_dataset$cal_alphadiv(PD = FALSE)
mphlan_cis_dataset$cal_betadiv(unifrac = FALSE)
mphlan_cis_dataset$cal_abund()



mod_pos_pruned$cal_alphadiv(PD = FALSE)
mod_pos_pruned$cal_betadiv(unifrac = FALSE)
mod_pos_pruned$cal_abund()






mod_pos_edited$cal_alphadiv(PD = FALSE)
mod_pos_edited$cal_betadiv(unifrac = FALSE)
mod_pos_edited$cal_abund()


mod_negsamples_dataset$cal_alphadiv(PD = FALSE)
mod_negsamples_dataset$cal_betadiv(unifrac = FALSE)
mod_negsamples_dataset$cal_abund()



mod_control_cis_dataset$cal_alphadiv(PD = FALSE)
mod_control_cis_dataset$cal_betadiv(unifrac = FALSE)
mod_control_cis_dataset$cal_abund()



merge_new_dataset$cal_alphadiv(PD = FALSE)
merge_new_dataset$cal_betadiv(unifrac = FALSE)
merge_new_dataset$cal_abund()

merge_new_dataset_its$cal_alphadiv(PD = FALSE)
merge_new_dataset_its$cal_betadiv(unifrac = FALSE)
merge_new_dataset_its$cal_abund()




fungi_250_eco$cal_alphadiv(PD = FALSE)
fungi_250_eco$cal_betadiv(unifrac = FALSE)
fungi_250_eco$cal_abund()

batch_dataset$cal_alphadiv(measures = "Observed")
batch_dataset$cal_betadiv(unifrac = FALSE)
batch_dataset$cal_abund()

batch_control_dataset$cal_alphadiv(measures = "Observed")
batch_control_dataset$cal_betadiv(unifrac = FALSE)
batch_control_dataset$cal_abund()

batch_cis_dataset$cal_alphadiv(measures = "Observed")
batch_cis_dataset$cal_betadiv(unifrac = FALSE)
batch_cis_dataset$cal_abund()



batch_dataset <- phyloseq2meco(mphlan_dataset_physeq)
batch_control_dataset <- phyloseq2meco(mphlan_control)
batch_cis_dataset <- phyloseq2meco(mphlan_cis)
pre_dataset$cal_abund()

batch_dataset$cal_abund()
# ok now
pre_abundance <- trans_abund$new(dataset = pre_dataset, taxrank = "Phylum", ntaxa = 10)

pre_abundance$plot_bar(others_color = "grey70", facet = "Condition_disease", xtext_keep = FALSE, legend_text_italic = FALSE)


pre_dataset$cal_abund()
fungi_pre_dataset$cal_abund()


pre_dataset$cal_alphadiv(PD = FALSE)
fungi_pre_dataset$cal_alphadiv(measures = "Chao1")
pre_dataset$cal_betadiv(unifrac = FALSE)
fungi_pre_dataset$cal_betadiv(unifrac = FALSE)

fun_sam <- sample_data(fungi_selected_samples )


sample_data(fungi_selected_samples)[sample_data(fungi_selected_samples)=="RR"]<-"MS"
sample_data(fungi_selected_samples)[sample_data(fungi_selected_samples)=="SP"]<-"MS"


sample_data(fungi_selected_samples)[sample_data(fungi_selected_samples)==0]<-"MS"
sample_data(fungi_selected_samples)[sample_data(fungi_selected_samples)=="R"]<-"MS"


phyloseq_pos_neg <- phyloseq_posneg_sam
phyloseq_set <- phyloseq_new_set
phyloseq_mod <- phyloseq_sorted_mod


#Merging SP to RR and labelling as MS
phylomerge <- phyloseq_250_sub
phylomerge_its <- phyloseq_250_sub_its


sample_data(phylomerge)[sample_data(phylomerge)=="SP"]<-"RR"
sample_data(phylomerge)[sample_data(phylomerge)=="RR"]<-"MS"

sample_data(phylomerge_its)[sample_data(phylomerge_its)=="SP"]<-"RR"
sample_data(phylomerge_its)[sample_data(phylomerge_its)=="RR"]<-"MS"


sample_data(phyloseq_mod)[sample_data(phyloseq_mod)=="RR"]<-"CIS"


otu_table(phyloseq_pos_neg)[otu_table(phyloseq_pos_neg)==0]<-1
otu_table(phyloseq_mod)[otu_table(phyloseq_mod)==0]<-1
otu_table(phyloseq_set)[otu_table(phyloseq_set)==0]<-1
otu_table(phyloseq_pos_seq)[otu_table(phyloseq_pos_seq)==0]<-1
otu_table(phyloseq_pos_seq)[otu_table(phyloseq_pos_seq)==0]<-1

otu_table(bac_selected_phylo)[otu_table(bac_selected_phylo)==0]<-1
otu_table(fun_selected_phylo)[otu_table(fun_selected_phylo)==0]<-1





otu_table(phyloseq_250_sub)[otu_table(phyloseq_250_sub)==0]<-1
otu_table(phyloseq_250_sub_its)[otu_table(phyloseq_250_sub_its)==0]<-1

otu_table(modagg)[otu_table(modagg)==1]<-0

postax <- tax_table(phyloseq_mod)
postax <- as.data.frame(postax)


moddaggtransformed <- transform_sample_counts(modagg, function(OTU) OTU/sum(OTU))

otuaggtr <- t(otu_table(moddaggtransformed))
otuaggtr <- as.data.frame(otuaggtr)



otuposneg <- as.data.frame((otu_table(phyloseq_pos_neg)))
modagg <- aggregate_taxa(phyloseq_mod, "Species", verbose = FALSE)

aggtaxy <- tax_table(modagg)
taxagg <- as.data.frame(aggtaxy)


otuagg <- t(otu_table(modagg))
otuagg <- as.data.frame(otuagg)

otumphlan1 = otumphlan[-c(2),]
otumphlan2 = otumphlan1[-c(10),]
otumphlan3 = otumphlan2[-c(27),]
otumphlan4 = otumphlan3[-c(54),]



otu_table(mphlan_dataset_physeq)[otu_table(mphlan_dataset_physeq)==0.0000000]<-1
otu_table(mphlan_dataset_physeq)[otu_table(mphlan_dataset_physeq)==1]<-0.0000000


phyloseq_250_sub <- subset_samples(phyloseq_na_250,
                                   Condition %in% c("CTRL", "CIS", "RR", "SP"))

phyloseq_250_sub_its <- subset_samples(phyloseq_na_250_its,
                                   Condition %in% c("CTRL", "CIS", "RR", "SP"))




sample_tab <- sample_data(pre_selected_samples)
otu_tab <- t(otu_table(pre_selected_samples))
o <- as.data.frame(otu_tab)
tax_tab <- tax_table(pre_selected_samples)


sample_tab_250 <- sample_data(phyloseq_250_sub)
otu_tab_250 <- t(otu_table(phyloseq_250_sub))
tax_tab_250 <- tax_table(phyloseq_250_sub)

sample_tab_mod <- sample_data(phyloseq_mod)
otu_tab_mod <- t(otu_table(phyloseq_mod))
tax_tab_mod <- tax_table(phyloseq_mod)




sample_tab_pos <- sample_data(phyloseq_pos_neg)
otu_tab_pos <- t(otu_table(phyloseq_pos_neg))
tax_tab_pos <- tax_table(phyloseq_pos_neg)

sample_tab_pos <- sample_data(phyloseq_pos_seq)
otu_tab_pos <- t(otu_table(phyloseq_pos_seq))
tax_tab_pos <- tax_table(phyloseq_pos_seq)





sample_tab_set <- sample_data(phylomerge)
otu_tab_set <- t(otu_table(phylomerge))
tax_tab_set <- tax_table(phylomerge)

sample_tab_sel <- sample_data(bac_selected_phylo)
otu_tab_sel <- t(otu_table(bac_selected_phylo))
tax_tab_sel <- tax_table(bac_selected_phylo)

sample_tab_sel_its <- sample_data(fun_selected_phylo)
otu_tab_sel_its <- t(otu_table(fun_selected_phylo))
tax_tab_sel_its <- tax_table(fun_selected_phylo)




sample_tab_set_its <- sample_data(phylomerge_its)
otu_tab_set_its <- t(otu_table(phylomerge_its))
tax_tab_set_its <- tax_table(phylomerge_its)




sample_tab_250_its <- sample_data(phyloseq_250_sub_its)
otu_tab_250_its <- t(otu_table(phyloseq_250_sub_its))
tax_tab_250_its <- tax_table(phyloseq_250_sub_its)

sample_tab_250_batch <- sample_data(phyloseq_mphlan_batch)
otu_tab_250_batch <- t(otu_table(phyloseq_250_sub_its))
tax_tab_250_batch <- tax_table(phyloseq_250_sub_its)


phyloseq_pos_seq

otumphlan <- t(otu_table(mphlan_dataset_physeq))
otumphlan <- as.data.frame(otumphlan)

taxmphlan <- (tax_table(mphlan_dataset_physeq))
taxmphlan <- as.data.frame(taxmphlan)

fun_sample_tab <- sample_data(fungi_selected_samples)
fun_otu_tab <- t(otu_table(fungi_selected_samples))
fun_o <- as.data.frame(otu_tab)
fun_tax_tab <- tax_table(fungi_selected_samples)



batch_sample_tab <- sample_data(phyloseq_mphlan_batch)
batch_otu_tab <- t(otu_table(phyloseq_mphlan_batch))
batch_tax_tab <- tax_table(phyloseq_mphlan_batch)




pre_dataset <- microtable$new(sample_table = sample_tab, otu_table = otu_tab, tax_table = tax_tab)
pre_dataset

fungi_pre_dataset <- microtable$new(sample_table = fun_sample_tab, otu_table = fun_otu_tab, tax_table = fun_tax_tab)
fungi_pre_dataset

batch_pre_dataset <- microtable$new(sample_table = batch_sample_tab, otu_table = batch_otu_tab, tax_table = batch_tax_tab)
batch_pre_dataset

bacteria_250_dataset <- microtable$new(sample_table = sample_tab_250, otu_table = otu_tab_250, tax_table = tax_tab_250)
bacteria_250_dataset

pos_dataset <- microtable$new(sample_table = sample_tab_pos, otu_table = otu_tab_pos, tax_table = tax_tab_pos)
pos_dataset


merge_new_set_dataset <- microtable$new(sample_table = sample_tab_set, otu_table = otu_tab_set, tax_table = tax_tab_set)
merge_new_set_dataset 


merge_new_set_dataset_its <- microtable$new(sample_table = sample_tab_set_its, otu_table = otu_tab_set_its, tax_table = tax_tab_set_its)
merge_new_set_dataset_its


new_set_dataset <- microtable$new(sample_table = sample_tab_mod, otu_table = otu_tab_mod, tax_table = tax_tab_mod)
new_set_dataset 

bac_sel_dataset <- microtable$new(sample_table = sample_tab_sel, otu_table = otu_tab_sel, tax_table = tax_tab_sel)
bac_sel_dataset 

fun_sel_dataset <- microtable$new(sample_table = sample_tab_sel_its, otu_table = otu_tab_sel_its, tax_table = tax_tab_sel_its)
fun_sel_dataset



fungi_250_dataset <- microtable$new(sample_table = sample_tab_250_its, otu_table = otu_tab_250_its, tax_table = tax_tab_250_its)
fungi_250_dataset



pre_dataset$cal_alphadiv(PD = FALSE)
pre_dataset$cal_betadiv(unifrac = FALSE)
pre_dataset$cal_abund()

t_alpha <- trans_alpha$new(dataset = pre_dataset, by_group = "Condition_disease_variant" )

trans_alpha$new


saveRDS(pre_dataset, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/pre_dataset.rds")



t_alpha$data_alpha

t_alpha$cal_diff(method = "t.test")
# return t1$res_diff
head(t1$res_diff)


pre_data_abund<- pre_dataset$cal_abund()

tax_count <- pre_dataset$taxa_abund

pre_dataset$cal_betadiv(unifrac = FALSE)

pre_dataset$alpha_diversity

sample_tab_fun <- sample_data(fungi_selected_samples)
otu_tab_fun <- t(otu_table(fungi_selected_samples))
tax_tab_fun <- tax_table(fungi_selected_samples)

fungi_dataset <- microtable$new(sample_table = sample_tab_fun, otu_table = otu_tab_fun, tax_table = tax_tab_fun)



#Taxonomic_abundance
save(pre_dataset, file = "/Users/cimi_bioinformatics/Desktop/pre_dataset.RData")
zip("/Users/cimi_bioinformatics/Desktop/pre_dataset.zip", "/Users/cimi_bioinformatics/Desktop/pre_dataset.RData")

pre_dataset$sample_table$Condition_disease_variant %<>% factor(., levels = c("CTRL", "CIS", "RR", "SP"))
fungi_pre_dataset$sample_table$DiseaseEvolution %<>% factor(., levels = c("CTRL", "CIS", "RR", "SP"))

bacteria_250_eco$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS", "RR", "SP"))
fungi_250_eco$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS", "RR", "SP"))
batch_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))
batch_dataset$sample_table$Type %<>% factor(., levels = c("Pos", "Neg"))
batch_cis_dataset$sample_table$Type %<>% factor(., levels = c("Pos", "Neg"))
batch_control_dataset$sample_table$Type %<>% factor(., levels = c("Pos", "Neg"))
posseq_dataset$sample_table$Notes.1 %<>% factor(., levels = c("IgApos_IgGpos", "IgAneg_IgGneg"))
new_dataset$sample_table$Study %<>% factor(., levels = c("IgA_pos", "IgA_neg"))
merge_new_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS", "MS"))
merge_new_dataset_its$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS", "MS"))
mod_control_dataset$sample_table$Study %<>% factor(., levels = c("IgA_pos", "IgA_neg"))
mod_cis_dataset$sample_table$Study %<>% factor(., levels = c("IgA_pos", "IgA_neg"))
mod_control_cis_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))

mod_pos_pruned$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))
mod_pos_edited$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))
mod_pos_edited$sample_table$Study %<>% factor(., levels = c("IgA_pos", "IgA_neg"))
mod_negsamples_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))
mphlan_pos_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))
mphlan_neg_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS"))

new_dataset_posneg$sample_table$Study %<>% factor(., levels = c("IgA_pos", "IgA_neg"))

bacteria_selected_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS", "MS"))
fungi_selected_dataset$sample_table$Condition %<>% factor(., levels = c("CTRL", "CIS", "MS"))

saveRDS(bacteria_selected_dataset, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/bacteria_selected_dataset.rds")

saveRDS(fungi_selected_dataset, "/Users/cimi_bioinformatics/Desktop/Phyloseq_Objects/fungi_selected_dataset.rds")

mphlan_control_dataset$sample_table$Type %<>% factor(., levels = c("Pos", "Neg"))
mphlan_cis_dataset$sample_table$Type %<>% factor(., levels = c("Pos", "Neg"))

str(pre_dataset$sample_table)

pre_abundance <- trans_abund$new(dataset = pre_dataset, taxrank = "Species", ntaxa = 20, )
fungi_pre_abundance <- trans_abund$new(dataset = fungi_pre_dataset, taxrank = "Species", ntaxa = 20, )
data_abund <- pre_abundance$data_abund


batch_abundance <- trans_abund$new(dataset = batch_dataset, taxrank = "Species", ntaxa = 20, )
batch_data_abund <- batch_abundance$data_abund



pre_abundance$plot_bar(others_color = "grey70", facet = "Condition_disease_variant", xtext_keep = FALSE, legend_text_italic = FALSE)

pre_abundance$plot_heatmap(facet = "Condition_disease_variant", xtext_keep = FALSE, withmargin = FALSE)
batch_abundance$plot_heatmap(facet = "Type", xtext_keep = FALSE, withmargin = FALSE, ytext_size = 16)


fungi_pre_abundance$plot_heatmap(facet = "DiseaseEvolution", xtext_keep = FALSE, withmargin = FALSE, plot_text_size = 20)


t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)


#Beta diversity

beta_pre_div <- trans_beta$new(dataset = bacteria_250_eco, measure = "bray", group = "Condition")
beta_pre_div$cal_ordination(ordination = "PCoA")
beta_pre_div$plot_ordination(plot_color = "Condition", plot_shape = "Condition", plot_type = "point") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2)


beta_fun_div <- trans_beta$new(dataset = fungi_250_eco, measure = "bray", group = "Condition")
beta_fun_div$cal_ordination(ordination = "PCoA")
beta_fun_div$plot_ordination(plot_color = "Condition", plot_shape = "Condition", plot_type = c("point", "ellipse")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


beta_batch_div <- trans_beta$new(dataset = mphlan_cis_dataset, measure = "bray", group = "Type")
beta_batch_div$cal_ordination(ordination = "PCoA")
beta_batch_div$plot_ordination(plot_color = "Type", plot_shape = "Type", color_values = RColorBrewer::brewer.pal(n = 6, name = "Set1"), centroid_segment_linetype = 3, plot_type = c("point", "ellipse")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))  + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) +
  theme(axis.text.y=element_text(size=15)) + theme(axis.text.x=element_text(size=15)) + 
  theme(axis.title.y=element_text(size=17)) + theme(axis.title.x=element_text(size=17)) 



beta_pre_div$cal_group_distance(within_group = TRUE)
beta_pre_div$cal_group_distance_diff(method = "anova")
beta_pre_div$plot_group_distance(boxplot_add = "mean")
beta_pre_div$plot_group_distance(plot_group_order = c("CTRL", "CIS", "RR", "SP"), boxplot_add = "mean")


beta_fun_div$cal_group_distance(within_group = TRUE)
beta_fun_div$cal_group_distance(within_group = FALSE)
beta_fun_div$cal_group_distance_diff(method = "wilcox")
beta_fun_div$plot_group_distance(plot_group_order = c("CTRL", "CIS", "RR", "SP"), boxplot_add = "mean")


beta_batch_div$cal_group_distance(within_group = TRUE)
beta_batch_div$cal_group_distance_diff(method = "wilcox")
beta_pre_div$plot_group_distance(boxplot_add = "mean")
beta_batch_div$plot_group_distance(plot_group_order = c("CTRL", "CIS"), boxplot_add = "mean")





beta_pre_div$plot_clustering(group = "Condition")
beta_fun_div$plot_clustering(group = "Condition")


pre_dataset$sample_table

#Alpha_Diversity
alpha_pre_div <- trans_alpha$new(dataset = pre_dataset, group = "Condition_disease_variant", taxa_level )
alpha_fun_div <- trans_alpha$new(dataset = fungi_pre_dataset, group = "DiseaseEvolution")
alpha_merge_div <- trans_alpha$new(dataset = merge_new_dataset_its, group = "Condition")
alpha_mod_div <- trans_alpha$new(dataset = mod_negsamples_dataset, group = "Condition")


alpha_bacsel_div <- trans_alpha$new(dataset = bacteria_selected_dataset, group = "Condition")
alpha_funsel_div <- trans_alpha$new(dataset = fungi_selected_dataset, group = "Condition")

alpha_pos_div <- trans_alpha$new(dataset = mphlan_pos_dataset, group = "Condition")
alpha_neg_div <- trans_alpha$new(dataset = mod_negsamples_dataset, group = "Condition")

alpha_ctrl_div <- trans_alpha$new(dataset = mphlan_control_dataset, group = "Type")
alpha_cis_div <- trans_alpha$new(dataset = mod_cis_dataset, group = "Study")


alpha_bac_250 <- trans_alpha$new(dataset = bacteria_250_eco, group = "Condition")
alpha_bac_250level <- trans_alpha$new(dataset = bacteria_250_eco, group = "Condition", taxa_level = "Species")

alpha_fun_250 <- trans_alpha$new(dataset = fungi_250_eco, group = "Condition")

head(alpha_pre_div$data_stat)

alpha_pre_div$data_stat$
alpha_fun_div$data_stat
alpha_pos_div$data_stat


alpha_bac_250$data_stat
alpha_fun_250$data_stat

alpha_pre_div$cal_diff(method = "")
alpha_fun_div$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE, measure = "Chao1")

alpha_mod_div$cal_diff(method = "wilcox")
alpha_merge_div$cal_diff(method = "wilcox")

alpha_ctrl_div$cal_diff(method = "wilcox")
alpha_cis_div$cal_diff(method = "wilcox")

alpha_bacsel_div$cal_diff(method = "wilcox")
alpha_funsel_div$cal_diff(method = "wilcox", measure = "Shannon")

alpha_pos_div$cal_diff(method = "wilcox")
alpha_pos_div$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)

alpha_bac_250$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
alpha_bac_250$cal_diff(method = "wilcox")

alpha_fun_250$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE, measure ="Simpson")
alpha_fun_250$cal_diff(method = "wilcox", measure = "Shannon")


significance_fungi_amplicon <- as.data.frame(alpha_fun_250$res_diff)
write.csv(significance_fungi_amplicon, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_fungi_amplicon.csv")


significance_merge_bacteria_amplicon <- as.data.frame(alpha_merge_div$res_diff)
write.csv(significance_merge_bacteria_amplicon, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_merge_bacteria_amplicon.csv")


significance_merge_fungi_amplicon <- as.data.frame(alpha_merge_div$res_diff)
write.csv(significance_merge_fungi_amplicon, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_merge_fungi_amplicon.csv")

alpha_pre_div$plot_alpha(measure = "Chao1", shape = "Condition_disease", add_sig = TRUE, order_x = c("CTRL", "CIS", "MS") ) 
  #alpha_pre_div$plot_alpha(measure = "Chao1", shape = "Condition_disease
  

significance_bacteria_amplicon <- as.data.frame(alpha_bac_250$res_diff)
write.csv(significance_bacteria_amplicon, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_bacteria_amplicon.csv")

significance_neg_controlcis <- as.data.frame(alpha_mod_div$res_diff)
write.csv(significance_neg_controlcis, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_neg_controlcis.csv")

significance_posneg_cis <- as.data.frame(alpha_mod_div$res_diff)
write.csv(significance_posneg_cis, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_posneg_cis.csv")

significance_controlcis <- as.data.frame(alpha_mod_div$res_diff)
write.csv(significance_controlcis, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_controlcis.csv")

significance_controlcisms <- as.data.frame(alpha_bacsel_div$res_diff)
write.csv(significance_controlcisms, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_controlcisms_bacsel.csv")

significance_controlcisms_fungi <- as.data.frame(alpha_funsel_div$res_diff)
write.csv(significance_controlcisms_fungi, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_controlcisms_funsel.csv")

significance_controlposneg <- as.data.frame(alpha_ctrl_div$res_diff)
write.csv(significance_controlposneg, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_controlposneg.csv")

significance_cisposneg <- as.data.frame(alpha_cis_div$res_diff)
write.csv(significance_cisposneg, "/Users/cimi_bioinformatics/Desktop/Bioinf/significance_cisposneg.csv")



plot_bac_shannon <- alpha_bacsel_div$plot_alpha(measure = "Shannon", shape = "Study", add_sig_text_size = 10, xlab = "Sample",
                                 color_values = RColorBrewer::brewer.pal(n = 6, name = "Set1")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.text.x = element_text(size = 25)) 



plot_bac_shannon



plot_bac_shannon <- alpha_funsel_div$plot_alpha(measure = "Shannon", shape = "Condition", add_sig_text_size = 10, xlab = "Sample",
                                                color_values = RColorBrewer::brewer.pal(n = 7, name = "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.text.x = element_text(size = 25)) 



plot_bac_shannon







plot_bac_Observed <- alpha_mod_div$plot_alpha(measure = "Observed", shape = "Condition", xlab = "Sample",
                                         color_values = RColorBrewer::brewer.pal(7, "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.text.x = element_text(size = 25))

plot_bac_Observed


plot_bac_Chao <- alpha_bacsel_div$plot_alpha(measure = "Chao1", shape = "Condition", xlab = "Sample",
                                              color_values = RColorBrewer::brewer.pal(7, "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.text.x = element_text(size = 25))


plot_bac_Chao


plot_bac_Simpson <- alpha_bac_250$plot_alpha(measure = "Simpson", shape = "Condition", xlab = "Sample",
                                          color_values = RColorBrewer::brewer.pal(7, "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.title.x=element_text(size=25))


plot_bac_Simpson


#Fungi_Alpha_Diversity_Plots


plot_fun_shannon <- alpha_fun_250$plot_alpha(measure = "Shannon", shape = "Condition", xlab = "Sample",
                                         color_values = RColorBrewer::brewer.pal(7, "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.title.x=element_text(size=25))

plot_fun_shannon


plot_fun_observed <- alpha_fun_250$plot_alpha(measure = "Observed", shape = "Condition", xlab = "Sample",
                                         color_values = RColorBrewer::brewer.pal(7, "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.title.x=element_text(size=25))

plot_fun_observed

plot_fun_simpson <- alpha_fun_250$plot_alpha(measure = "Simpson", shape = "Condition", xlab = "Sample",
                                              color_values = RColorBrewer::brewer.pal(7, "Dark2")) +
  theme(axis.title.y=element_text(size=25,face="bold")) + theme(axis.title.x=element_text(size=25))

plot_fun_simpson




plot_chao1 <- alpha_pre_div$plot_alpha(measure = "Chao1", shape = "Condition_disease_variant", xlab = "Sample",
                                         color_values = RColorBrewer::brewer.pal(7, "Dark2"))

plot_chao1


#Differential abundance
fungi_pre_diff <- trans_diff$new(dataset = fungi_pre_dataset, method = "lefse", group = "DiseaseEvolution", alpha = 0.01,
                           p_adjust_method = "none")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
fungi_pre_diff$plot_diff_bar(threshold = 0.5)
# we show 20 taxa with the highest LDA (log10)
fungi_pre_diff$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("CTRL", "CIS", "RR", "SP"))


# clade_label_level 5 represent phylum level in this analysis
# require ggtree package
fungi_pre_diff$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 50, clade_label_level = 5, group_order = c("CTRL", "RR", "CIS", "SP"))
                             


# Differential Abundance
t1_bac_250 <- trans_diff$new(dataset = bacteria_250_eco, method = "DESeq2", group = "Condition", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Species")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1_bac_250$plot_diff_bar(threshold = 4)
# we show 20 taxa with the highest LDA (log10)
t1_bac_250$plot_diff_bar(use_number = 1:20, width = 0.8, group_order = c("CTRL", "CIS", "RR", "SP"))

t1_bac_250$res_diff %<>% subset(Significance %in% "***")
t1_bac_250$plot_diff_abund(use_number = 1:30, add_sig = T, add_sig_label = "Significance", select_taxa = t1_bac_250$plot_diff_bar_taxa)

t1_bac_250$plot_diff_abund(group_order = c("CTRL", "CIS", "RR", "SP"), select_taxa = t1_bac_250$plot_diff_bar_taxa,
                           add_sig = T, add_sig_label = "Significance")


t1_bac_250$plot_diff_abund(add_sig = T, add_sig_label = "Significance")


t1_bac_250$plot_diff_abund(use_number = 1:15,add_sig = T, add_sig_label = "Significance", group_order = c("CTRL", "CIS", "RR", "SP")) +
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))


t1_bac_250$plot_diff_bar(width = 0.4, select_group = "CTRL - RR") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))



t1 <- trans_diff$new(dataset = bacteria_250_eco, method = "wilcox", group = "Condition", taxa_level = "Species", filter_thres = 0.001)
# filter something not needed to show
t1$res_diff %<>% subset(Significance %in% "***")
t1$plot_diff_abund(use_number = 1:20, add_sig = T, add_sig_label = "Significance")


#Bacteria


t12f_pos <- trans_diff$new(dataset = fungi_selected_dataset, method = "ancombc2", group = "Condition", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Species", p_adjust_method = "BH")

theme(axis.title.y=element_text(size=25,face="bold")) 

# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t12_pos$plot_diff_bar(select_group = "IgAneg_IgGneg - IgApos_IgGpos")
# we show 20 taxa with the highest LDA (log10)
t1_pos$plot_diff_bar(width = 0.4, select_group = "IgAneg_IgGneg - IgApos_IgGpos") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))

t12_pos$res_diff %<>% subset(Significance %in% c("*","**","***","****"))
t12_pos$plot_diff_abund(add_sig = T, add_sig_label = "Significance") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))


t12_pos$plot_diff_abund(add_sig = T, add_sig_label = "Significance", coord_flip = TRUE, 
                        barwidth= 0.7, color_values = RColorBrewer::brewer.pal(n = 6, name = "Set1")) + 
  theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size = 20) +
                                                         theme(axis.title.y=element_text(size=25,face="italic"))) 
                                                       

t12f_pos$plot_diff_abund(add_sig = T, add_sig_label = "Significance", ) + 
  theme(axis.text.y = element_text(size = 20, face = "italic")) + theme(axis.text.x = element_text(size = 15)) 



###For funghi in dystopian story 


#Fungi

t1_fun_250 <- trans_diff$new(dataset = fungi_250_eco, method = "DESeq", group = "Condition", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Species")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1_fun_250$plot_diff_bar(select_group = "CTRL - RR")
# we show 20 taxa with the highest LDA (log10)
t1_fun_250$plot_diff_bar(width = 0.4, select_group = "CTRL - SP") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))

t1_fun_250$res_diff %<>% subset(Significance %in% "***")
t1_fun_250$plot_diff_abund(add_sig = T, add_sig_label = "Significance") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))


t1_fun_250$plot_diff_abund(add_sig = T, add_sig_label = "Significance")


#Bacteria Shotgun


t1_batch <- trans_diff$new(dataset = bacteria_selected_dataset, method = "anc", group = "Type", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Species")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1_batch$plot_diff_bar()
# we show 20 taxa with the highest LDA (log10)
t1_batch$plot_diff_bar(width = 0.4, select_group = "Neg - Pos") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))

t1_batch$res_diff %<>% subset(Significance %in% "*")
t1_batch$plot_diff_abund(add_sig = T, add_sig_label = "Significance") + 
  theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15))




#Bacteria_Alpha_Diversity_Amplicon

#CTRL - 6.794294e+00
#CIS - 6.817368e+00
#MS. - 6.586964e+00


#Fungi_Alpha_Diversity_Amplicon

#CTRL - 2.252920e+00
#CIS - 2.026261e+00
#MS - 1.948724e+00


##PosNeg_ControlCIS
#Control - 7.627966e+00
#CIS - 7.688442e+00

#PosNeg_Control
#Pos - 7.631946e+00
#Neg - 7.623701e+00


#PosNeg_CIS
#Pos - 7.805989e+00
#Neg - 7.551304e+00















