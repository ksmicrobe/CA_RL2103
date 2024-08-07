library(microbiomeMarker)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggplot2)
library(microViz)
library(ape)
library(mia)
library(ggtree)
library(pheatmap)
library(RColorBrewer)
library(rcartocolor)
library(rcompanion)
library(vegan)
library(wesanderson)
library(readxl)
library(otu2ot)
library(bios2mds)
library(ANCOMBC)
library(speedyseq)
library(msa)

###SECTION ONE: Prepping data & PERMANOVAs###
marine_counts_taxa <- data.frame(readRDS("total_taxa_counts_marine")) %>% 
  dplyr::rename(Species = sequence) #turn sequence into "Species" for phyloseq incorporation
marine_counts <- data.frame(marine_counts_taxa[,8:85])
marine_taxa <- data.frame(marine_counts_taxa[,1:7])
marine_counts <- as.matrix(marine_counts)
marine_taxa <- as.matrix(marine_taxa)

CA_phylo <- phyloseq(otu_table(marine_counts, taxa_are_rows = TRUE), 
                     tax_table(marine_taxa))
marine_phylo <- prune_taxa(taxa_sums(CA_phylo) > 0, CA_phylo) #Trim zero-count ASVs (zeros in all samples)

marine_tax <- data.frame(tax_table(marine_phylo))
rownames(marine_tax) <- 1:length(rownames(marine_tax))

marine_asv <- data.frame(otu_table(marine_phylo))
rownames(marine_asv) <- rownames(marine_tax)

CA_marine_meta <- data.frame(read_xlsx("~/OneDrive/Documents/CA_marine_meta_dTRC.xlsx"))
rownames(CA_marine_meta) <- CA_marine_meta$sample
CA_marine_meta <- CA_marine_meta %>% mutate(station = as.character(station))

ASVs <- lessR::to("ASV", nrow(marine_asv))
rownames(marine_asv) <- ASVs
rownames(marine_tax) <- ASVs

saveRDS(marine_asv, "raw_counts_RL2103")
saveRDS(marine_tax, "silva_phytoref_taxa_RL2103")

marine_asv <- readRDS("raw_counts_RL2103")
marine_tax <- readRDS("silva_phytoref_taxa_RL2103")

CA_marine_ps <- phyloseq(otu_table(as.matrix(marine_asv), taxa_are_rows = T), 
                         tax_table(as.matrix(marine_tax)), 
                         sample_data(CA_marine_meta))

CA_marine_ps <- microbiomeMarker::normalize(CA_marine_ps, "TSS") #TSS normalization

CA_marine_ps <- subset_taxa(CA_marine_ps, Kingdom != "NA")
CA_marine_ps <- microViz::tax_fix(CA_marine_ps, unknowns = "NA") #All unkowns auto-assigned to highest tax. level
marine_genus <- microbiomeMarker::aggregate_taxa(CA_marine_ps, "Genus")
subset_taxa(marine_genus, Genus != "Unknown") #No unknown taxa left, woohoo!

dist <- phyloseq::distance(marine_genus, method = "bray")

adonis2(dist ~ line, data = data.frame(sample_data(marine_genus))) #0.001
adonis2(dist ~ upwell_strength, data = data.frame(sample_data(marine_genus))) #0.001
adonis2(dist ~ depth_cat_1, data = data.frame(sample_data(marine_genus))) #0.001
adonis2(dist ~ depth_cat_2, data = data.frame(sample_data(marine_genus))) #0.001
###


###SECTION TWO: Station-by-station heatmap###

marine_counts <- data.frame(otu_table(marine_phylo)) #need raw counts, not TSS-relative abundance
rownames(CA_marine_meta) <- colnames(marine_counts) 

marine_taxa <- data.frame(tax_table(marine_phylo))
rownames(marine_taxa) <- rownames(marine_counts)

marine_counts <- as.matrix(marine_counts)
marine_taxa <- as.matrix(marine_taxa)

marine_phylo_reads <- phyloseq(otu_table(marine_counts, taxa_are_rows = T), 
                               tax_table(marine_taxa), 
                               sample_data(CA_marine_meta))


marine_phylo_reads_station <- merge_samples(marine_phylo_reads,"station") #samples merged by station and now we can TSS transform
otu_table(marine_phylo_reads_station) <- t(otu_table(marine_phylo_reads_station)) #if you don't do this then your otu_table will be incorrectly oriented!
marine_phylo_reads_station_tss <- microbiomeMarker::normalize(marine_phylo_reads_station, "TSS") #now they are normalized for unequal read depth by station

marine_phylo_reads_station_tss <- microViz::tax_fix(marine_phylo_reads_station_tss, unknowns = "NA")
station_genus <- microbiomeMarker::aggregate_taxa(marine_phylo_reads_station_tss, "Genus")
top_genus_station <- names(sort(taxa_sums(station_genus), 
                                decreasing = TRUE))[1:38] #~70%-93% of total relative abundance per station
top_genus_marine <- prune_taxa(top_genus_station, station_genus)
marine_genus_treesum <- makeTreeSummarizedExperimentFromPhyloseq(top_genus_marine)
mat <- assay(marine_genus_treesum)

taxa_hclust <- hclust(vegdist(mat, "bray"), method = "complete") #perform hierarchical clustering for top genera based on dissimilarity matrix
taxa_tree <- as.phylo(taxa_hclust)
taxa_tree <- ggtree(taxa_tree) + theme(plot.margin=margin(0,0,0,0))
taxa_ordered <- get_taxa_name(taxa_tree)

station_hclust <- hclust(vegdist(t(mat), "bray"), method = "complete") #perform hierarchical clustering for stations
station_tree <- as.phylo(station_hclust)
station_tree <- ggtree(station_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins
stations_ordered <- rev(get_taxa_name(station_tree))

# Creates clusters
station_clusters <- factor(cutree(tree = station_hclust, k = 6)) #clusters are essentially equal to monophyletic groups- picked from dendogram

# Converts into data frame
station_data <- data.frame(clusters = station_clusters)

# Order data so that it's same as in phylo tree
station_data <- station_data[stations_ordered, , drop = FALSE]
station_data$"line" <- c("FAR", "FAR", "SUR", "MON", "MON", 
                         'NAV', "MON", "NAV", "FAR", "FAR", "FAR", 
                         "NAV", "NAV", "NAV", "SUR", "SUR", "MON", 
                         "MON", "MON", "SUR")
station_data$line <- factor(station_data$line, levels = c("NAV", "FAR", "MON", "SUR"))

station_data <- station_data %>% dplyr::rename(Line = line, Cluster = clusters)
station_data$Cluster <- factor(station_data$Cluster, levels = c(2,3,5,1,4,6)) #colors follow clusters (cluster ordering not important)

annoCol<-list(Line=c(NAV="#E41A1C", FAR="#377EB8", MON="#4DAF4A", SUR="#984EA3"))

(stations_heatmap <- pheatmap(100*mat, annotation_col = station_data, annotation_colors = annoCol, 
                              legend_breaks = c(1, 5, 10, 15, 20), cluster_rows = taxa_hclust, 
                              cluster_cols = station_hclust, 
                              color = rev(brewer.pal(10, 'Spectral')), 
                              legend_labels = c('1%', '5%', '10%', '15%', '20%')))
  
ggsave(stations_heatmap, filename = "stations_dendogram_heatmap.png", height = 10, width = 8)


trcs <- CA_marine_meta %>% select(sample, station, B1, HMP, AmMP, cHET, HET)
trcs_station <- trcs %>% group_by(station) %>%
  select(-sample) %>% 
  summarise_each(funs(mean=mean(., na.rm = TRUE)))
trcs_station <- data.frame(trcs_station)
rownames(trcs_station) <- trcs_station$station
trcs_station <- trcs_station %>% select(-station)
trcs_station <- as.matrix(trcs_station)
trc_hclust <- hclust(vegdist(trcs_station, "bray"), method = "complete")
trc_tree <- as.phylo(trc_hclust)
trc_tree <- ggtree(trc_tree) + theme(plot.margin=margin(0,0,0,0))

trc_station_clusters <- factor(cutree(tree = trc_hclust, k = 6)) #clusters are essentially equal to monophyletic groups- picked from dendogram
trc_station_data <- data.frame(clusters = trc_station_clusters) %>% 
  arrange(clusters)

trc_station_data$"line" <- c("FAR", "FAR", "SUR", "FAR", "FAR", 
                         'FAR', "NAV", "NAV", "NAV", "MON", "MON", 
                         "MON", "MON", "MON", "MON", "NAV", "SUR", 
                         "SUR", "NAV", "SUR")
trc_station_data$line <- factor(trc_station_data$line, levels = c("NAV", "FAR", "MON", "SUR"))

trc_station_data <- trc_station_data %>% dplyr::rename(Line = line, Cluster = clusters)


stations_hclust <- hclust(vegdist(t(trcs_station), "bray"), method = "complete") #perform hierarchical clustering for stations
stations_tree <- as.phylo(stations_hclust)
stations_tree <- ggtree(stations_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0))
print(stations_tree)

annoCol<-list(Line=c(NAV="#E41A1C", FAR="#377EB8", MON="#4DAF4A", SUR="#984EA3"))

trcs_station <- t(trcs_station)

pheatmap(trcs_station, annotation_col = trc_station_data, annotation_colors = annoCol, 
         cluster_rows = stations_hclust, cluster_cols = trc_hclust)



###SECTION THREE: Stacked barplots###

bar_pal <- rcartocolor::carto_pal(12, "Safe")

station_order <- c("SUR", "MON", "FAR", "NAV")
depth_order <- c("surface", "mixed", "dcm", "deep")

marine_phylo_reads <- tax_fix(marine_phylo_reads)

marine_phylo_reads_line <- merge_samples(marine_phylo_reads,"line")
otu_table(marine_phylo_reads_line) <- t(otu_table(marine_phylo_reads_line))
marine_phylo_reads_line_tss <- microbiomeMarker::normalize(marine_phylo_reads_line, "TSS") 

marine_phylo_reads_depth <- merge_samples(marine_phylo_reads,"depth_cat_2")
otu_table(marine_phylo_reads_depth) <- t(otu_table(marine_phylo_reads_depth))
marine_phylo_reads_depth_tss <- microbiomeMarker::normalize(marine_phylo_reads_depth, "TSS") 


##Separate phyloseq objects for proks and euks -> line and depth 
CA_marine_ps_alg_line <- subset_taxa(marine_phylo_reads_line_tss, 
                                     Kingdom == "Eukaryota")
CA_marine_ps_proks_line <- subset_taxa(marine_phylo_reads_line_tss, 
                                       Kingdom != "Eukaryota")

CA_marine_ps_alg_depth <- subset_taxa(marine_phylo_reads_depth_tss, 
                                      Kingdom == "Eukaryota")
CA_marine_ps_proks_depth <- subset_taxa(marine_phylo_reads_depth_tss, 
                                        Kingdom != "Eukaryota")
##

CA_marine_ps_proks_line %>% 
  comp_barplot(tax_level = "Order", n_taxa = 14, sample_order = station_order, other_name = "Other Taxa") +
  labs(x = NULL, y = NULL) + theme(legend.position = "bottom", 
                                   axis.text.y = element_text(face = "bold", size = 12)) + 
  coord_flip() +
  scale_color_manual(values = bar_pal)

CA_marine_ps_alg_line %>% 
  comp_barplot(tax_level = "Order", n_taxa = 14, sample_order = station_order) +
  labs(x = NULL, y = NULL) + theme(legend.position = "bottom", 
                                   axis.text.y = element_text(face = "bold", size = 12)) + 
  coord_flip() +
  scale_color_manual(values = bar_pal)

CA_marine_ps_proks_depth %>% 
  comp_barplot(tax_level = "Order", n_taxa = 14, sample_order = rev(depth_order)) +
  labs(x = NULL, y = NULL) + theme(legend.position = "bottom", 
                                   axis.text.y = element_text(face = "bold", size = 12)) + 
  coord_flip() +
  scale_color_manual(values = bar_pal)

CA_marine_ps_alg_depth %>% 
  comp_barplot(tax_level = "Order", n_taxa = 14, sample_order = rev(depth_order)) +
  labs(x = NULL, y = NULL) + theme(legend.position = "bottom", 
                                   axis.text.y = element_text(face = "bold", size = 12)) + 
  coord_flip() +
  scale_color_manual(values = bar_pal)
###


###SECTION FOUR: Beta diversity plots (microbial communities)

##Vectors plots
CA_div <- estimate_richness(CA_marine_ps, measures = c("Shannon", "Simpson", "InvSimpson"))
CA_meta <- data.frame(sample_data(CA_marine_ps))

CA_marine_ps_alg <- phyloseq::subset_taxa(marine_phylo_reads, Kingdom == "Eukaryota")
algal_richness <- colSums(otu_table(CA_marine_ps_alg))
CA_meta$"algae" <- algal_richness

vectors_data <- data.frame(shan_diversity = CA_div$Shannon, diversity = CA_div$InvSimpson, 
                           thiamine = CA_meta$B1, HMP = CA_meta$HMP, 
                           AmMP = CA_meta$AmMP, cHET = CA_meta$cHET, HET = CA_meta$HET, 
                           temp = CA_meta$temp, den = CA_meta$density, 
                           sal = CA_meta$sal, chl = CA_meta$chl, oxy = CA_meta$oxygen, 
                           trans = CA_meta$trans, irrad = CA_meta$irrad, 
                           algal_rich = CA_meta$algae)

rownames(vectors_data) <- rownames(CA_div)

vectors_data <- vectors_data %>% mutate(across(shan_diversity:algal_rich, as.numeric)) %>% 
  mutate(across(shan_diversity:algal_rich, transformTukey)) #Data now transformed by Tukey's Ladder of Powers

marine_genus_ord <- phyloseq::ordinate(marine_genus, "NMDS", "bray") %>% print()

set.seed(1577)
envfit_gen <- envfit(marine_genus_ord, 
                     vectors_data, permutations = 9999)
envfit_gen

en <- as.list(envfit_gen$vectors)
pvals.en <- as.data.frame(en$pvals)
arrows.en <- as.data.frame(en$arrows*sqrt(en$r))
C.en <- cbind(arrows.en, pvals.en)
Cred.en <- subset(C.en, pvals.en<0.05) #Only use significant vectors
Cred.en$"Variable" <- rownames(Cred.en)

data.scores.16S <- vegan::scores(marine_genus_ord)
data.scores.16S <- as.data.frame(data.scores.16S$sites)
data.scores.16S$"Line" <- CA_meta$line
data.scores.16S$Line <- factor(data.scores.16S$Line, 
                               levels = c("NAV", "FAR", "MON", "SUR"))
data.scores.16S$"Depth" <- CA_meta$depth_cat_2
data.scores.16S$Depth <- factor(data.scores.16S$Depth, levels = c("surface", "mixed", "deep", "dcm"))
data.scores.16S$"Upwelling" <- CA_meta$upwell_strength
data.scores.16S$Upwelling <- factor(data.scores.16S$Upwelling, levels = c("weak", "intermediate", "strong"))

rownames(data.scores.16S) <- CA_meta$sample

set1 <- RColorBrewer::brewer.pal(4, "Set1")

ggplot(data = data.scores.16S, aes(x = NMDS1, y = NMDS2))  +
  theme_classic() +
  geom_point(data = data.scores.16S, aes(shape = Upwelling, color = Line), size = 3) +
  scale_color_manual(values = set1) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = Cred.en, size =1, colour = "black", 
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text_repel(data = Cred.en, aes(x = NMDS1, y = NMDS2), 
                  colour = "black", size = 3.5, label = Cred.en$Variable, 
                  direction = "both")
##

##CAP ordination and Alpha Diversity

marine_genus_ord_cap <- ordinate(marine_genus, distance = "bray", 
                                 method = "CAP", 
                                 formula = ~depth_cat_2+line+upwell_strength)

depthpal <- c("#FD8D3C", "#E31A1C", "#800026", "#238B45")

plot_ordination(marine_genus, marine_genus_ord_cap, 
                type = "samples", color = "depth_cat_2") +  theme_classic() +
  geom_point(size = 4) +   
  scale_color_manual(values = depthpal, name = "Depth", limits = c("surface", "mixed", "deep", "dcm"), #without shape
                                              labels = c("Surface", "Intermediate", "Deep", "DCM")) +
  stat_ellipse(aes(group = line), type = "norm") 



plot_ordination(marine_genus, marine_genus_ord_cap, 
                type = "samples", , color = "depth_cat_2") +  theme_classic() +
  geom_point(size = 3.5, aes(shape = line)) +   
  scale_color_manual(values = depthpal, name = "Depth", limits = c("surface", "mixed", "deep", "dcm"), #with shape
                     labels = c("Surface", "Intermediate", "Deep", "DCM")) +
  scale_shape_manual(values = c(10,19,9,18)) +
  stat_ellipse(aes(group = line), type = "norm") 

##

###SECTION FOUR: dTRCs correlations and PCA

vectors_data %>% select('var_1', 'var_2' %>% pairs(lower.panel = NULL) #Used this to check for linear assocations 

#thiamine
cor.test(vectors_data$thiamine, vectors_data$den)
cor.test(vectors_data$thiamine, vectors_data$temp)
cor.test(vectors_data$thiamine, vectors_data$oxy)
cor.test(vectors_data$thiamine, vectors_data$chl)
cor.test(vectors_data$thiamine, vectors_data$sal)
cor.test(vectors_data$thiamine, vectors_data$diversity) #0.057
cor.test(vectors_data$thiamine, vectors_data$shan_diversity)
cor.test(vectors_data$thiamine, vectors_data$irrad)
cor.test(vectors_data$thiamine, vectors_data$trans)
cor.test(vectors_data$thiamine, vectors_data$algal_rich)
#

#HMP
cor.test(vectors_data$HMP, vectors_data$den)
cor.test(vectors_data$HMP, vectors_data$temp)
cor.test(vectors_data$HMP, vectors_data$oxy)
cor.test(vectors_data$HMP, vectors_data$chl)
cor.test(vectors_data$HMP, vectors_data$sal)
cor.test(vectors_data$HMP, vectors_data$diversity) 
cor.test(vectors_data$HMP, vectors_data$shan_diversity)
cor.test(vectors_data$HMP, vectors_data$irrad) #0.039
cor.test(vectors_data$HMP, vectors_data$trans)
cor.test(vectors_data$HMP, vectors_data$algal_rich)
#

#cHET
cor.test(vectors_data$cHET, vectors_data$den)
cor.test(vectors_data$cHET, vectors_data$temp)
cor.test(vectors_data$cHET, vectors_data$oxy)
cor.test(vectors_data$cHET, vectors_data$chl) #0.041
cor.test(vectors_data$cHET, vectors_data$sal)
cor.test(vectors_data$cHET, vectors_data$diversity) 
cor.test(vectors_data$cHET, vectors_data$shan_diversity)
cor.test(vectors_data$cHET, vectors_data$irrad)
cor.test(vectors_data$cHET, vectors_data$trans) #0.0011
cor.test(vectors_data$cHET, vectors_data$algal_rich) #0.00015
#

#AmMP
cor.test(vectors_data$AmMP, vectors_data$den)
cor.test(vectors_data$AmMP, vectors_data$temp)
cor.test(vectors_data$AmMP, vectors_data$oxy)
cor.test(vectors_data$AmMP, vectors_data$chl) #0.021
cor.test(vectors_data$AmMP, vectors_data$sal)
cor.test(vectors_data$AmMP, vectors_data$diversity) #0.0027
cor.test(vectors_data$AmMP, vectors_data$shan_diversity) 
cor.test(vectors_data$AmMP, vectors_data$irrad)
cor.test(vectors_data$AmMP, vectors_data$trans) #0.039
cor.test(vectors_data$AmMP, vectors_data$algal_rich)
#

#HET
cor.test(vectors_data$HET, vectors_data$den)
cor.test(vectors_data$HET, vectors_data$temp)
cor.test(vectors_data$HET, vectors_data$oxy) #0.059
cor.test(vectors_data$HET, vectors_data$chl) #0.031
cor.test(vectors_data$HET, vectors_data$sal)
cor.test(vectors_data$HET, vectors_data$diversity) #0.00065
cor.test(vectors_data$HET, vectors_data$shan_diversity)
cor.test(vectors_data$HET, vectors_data$irrad)
cor.test(vectors_data$HET, vectors_data$trans) #0.082
cor.test(vectors_data$HET, vectors_data$algal_rich)
#

#dTRC co-correlations
cor.test(vectors_data$thiamine, vectors_data$HMP)
cor.test(vectors_data$thiamine, vectors_data$AmMP) #7.24e-8
cor.test(vectors_data$thiamine, vectors_data$cHET)
cor.test(vectors_data$thiamine, vectors_data$HET)#3.06e-8
cor.test(vectors_data$HMP, vectors_data$AmMP)#0.0051
cor.test(vectors_data$HMP, vectors_data$cHET)
cor.test(vectors_data$HMP, vectors_data$HET)#0.047
cor.test(vectors_data$cHET, vectors_data$HET)#1.15e-5
cor.test(vectors_data$AmMP, vectors_data$cHET)#2.31e-7
cor.test(vectors_data$AmMP, vectors_data$HET)#6.831e-8
#

#NAV465_19 and FAR138_3 not included since there's no microbial community data for them

##Vectors plot
marine_trc_dist <- vegdist(as.matrix(marine_trc), "bray")
adonis2(marine_trc_dist ~ line, CA_meta, permutations = 9999) #.0004
adonis2(marine_trc_dist ~ depth_cat_1, CA_meta, permutations = 9999) 
adonis2(marine_trc_dist ~ depth_cat_2, CA_meta, permutations = 9999) 
adonis2(marine_trc_dist ~ upwell_strength, CA_meta, permutations = 9999) #.0004

marine_trc <- vectors_data %>% select(thiamine:HET)
marine_trc_pca <- vegan::rda(marine_trc)

set.seed(1579)
envfit_trc <- envfit(marine_trc_pca, 
                     vectors_data, permutations = 9999)
envfit_trc

en <- as.list(envfit_trc$vectors)
pvals.en <- as.data.frame(en$pvals)
arrows.en <- as.data.frame(en$arrows*sqrt(en$r))
C.en <- cbind(arrows.en, pvals.en)
Cred.en <- subset(C.en, pvals.en<0.05)
Cred.en$"Variable" <- rownames(Cred.en)

data.scores.16S <- vegan::scores(marine_trc_pca)
data.scores.16S <- as.data.frame(data.scores.16S$sites)
data.scores.16S$"Line" <- CA_meta$line
data.scores.16S$Line <- factor(data.scores.16S$Line, 
                               levels = c("NAV", "FAR", "MON", "SUR"))
data.scores.16S$"Upwelling" <- CA_meta$upwell_strength
data.scores.16S$Upwelling <- factor(data.scores.16S$Upwelling, 
                                    levels = c("weak", "intermediate", "strong"))
rownames(data.scores.16S) <- CA_meta$sample

en_coord_cont.16S = as.data.frame(vegan::scores(envfit_trc, "vectors")) * ordiArrowMul(envfit_trc) 

ggplot(data = data.scores.16S, aes(x = PC1, y = PC2))  +
  theme_classic() +
  geom_point(data = data.scores.16S, aes(shape = Upwelling, color = Line), size = 3) +
  scale_color_manual(values = set1) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               data = Cred.en, size =1, colour = "black", 
               arrow = arrow(length = unit(0.2,"cm")))
##
###

###SECTION FIVE: Differential abundance analysis and MaAsLin2###
marine_asv <- data.frame(readRDS("raw_counts_RL2103"))
marine_tax <- data.frame(readRDS("silva_phytoref_taxa_RL2103"))

marine_tax$Kingdom <- rownames(marine_tax)
marine_tax <- marine_tax[,-7]
marine_tax <- as.matrix(marine_tax)
marine_asv <- as.matrix(marine_asv)
rownames(marine_tax) <- rownames(marine_asv)

marine_phylo_reads <- phyloseq(otu_table(marine_asv, taxa_are_rows = T), 
                               tax_table(marine_tax), 
                               sample_data(CA_marine_meta))

sample_data(marine_phylo_reads)$depth_cat_2 <- factor(data.frame(sample_data(marine_phylo_reads))$depth_cat_2, 
                                       c("surface", "mixed", "dcm", "deep"))

reads_per_asv <- taxa_sums(marine_phylo_reads)
tot <- print(sum(reads_per_asv))
rare <- sum(reads_per_asv[reads_per_asv < 16]) #10.0% of reads

marine_phylo_reads_abund <- subset_taxa(marine_phylo_reads, taxa_sums(marine_phylo_reads) > 16)
marine_phylo_reads_abund <- tax_fix(marine_phylo_reads_abund, unknowns = "NA")

rare/tot #1% of ASV relative abundance

set.seed(1032)
ancom_res <- ancombc2(marine_phylo_reads_abund, tax_level = "Kingdom", 
                     fix_formula = "B1+cHET+HET+HMP+AmMP+depth_cat_2", 
                     verbose = TRUE, pseudo_sens = TRUE)
saveRDS(ancom_res, "ancombc2_RL2103_results")
ancom_res <- readRDS("ancombc2_RL2103_results")


ancom_res$res %>% filter(diff_B1 == 'TRUE') %>% dim() #0
ancom_res$res %>% filter(diff_HMP == 'TRUE') %>% dim() #2
ancom_res$res %>% filter(diff_cHET == 'TRUE') %>% dim() #48
ancom_res$res %>% filter(diff_HET == 'TRUE') %>% dim() #4
ancom_res$res %>% filter(diff_AmMP == 'TRUE') %>% dim() #0

diffabund_hmp <- data.frame(ancom_res$res %>% 
                              filter(diff_HMP == 'TRUE')) 
diffabund_chet <- data.frame(ancom_res$res %>% 
                               filter(diff_cHET == 'TRUE')) 
diffabund_het <- data.frame(ancom_res$res %>% 
                              filter(diff_HET == 'TRUE')) 
diffabund_depth <- data.frame(ancom_res$res %>% 
                                filter(diff_depth_cat_2mixed == 'TRUE' | diff_depth_cat_2dcm == 'TRUE' | diff_depth_cat_2deep == 'TRUE'))


diffabund_trc <- rbind(diffabund_hmp, diffabund_chet)
diffabund_trc <- rbind(diffabund_trc, diffabund_het)
diffabund_trc <- diffabund_trc %>% select(diff_cHET:diff_HMP, diff_depth_cat_2mixed:diff_depth_cat_2deep)

diffabund_trc$"Kingdom" <- rownames(diffabund_trc)
diffabund_asvs <- diffabund_trc$Kingdom

marine_taxa_diffabund <- data.frame(tax_table(marine_phylo_reads_abund)) %>% 
  filter(Kingdom %in% diffabund_asvs) %>% left_join(diffabund_trc, 'Kingdom')

marine_taxa_diffabund_cHET <- marine_taxa_diffabund %>% filter(diff_cHET == "TRUE")
marine_taxa_diffabund_HET <- marine_taxa_diffabund %>% filter(diff_HET == "TRUE")
marine_taxa_diffabund_HMP <- marine_taxa_diffabund %>% filter(diff_HMP == "TRUE")

top_genera <- data.frame(tax_table(top_genus_marine))

unique(marine_taxa_diffabund_cHET$Genus) %in% top_genera$Genus
unique(marine_taxa_diffabund_HET$Genus) %in% top_genera$Genus
unique(marine_taxa_diffabund_HMP$Genus) %in% top_genera$Genus

top_genera_trc <- top_genera %>% mutate(chet_diffabund = unique(top_genera$Genus) %in% marine_taxa_diffabund_cHET$Genus)

sum(top_genera_trc$chet_diffabund == 'TRUE') #25
sum(top_genera_trc$chet_diffabund == 'FALSE') #13

#25/38 top genera are differentially abundant with cHET
#Only 1/38 are differentially abundant with HET (and cHET, Roseibacillus), and none are with HMP

top_genera_chet <- top_genera_trc %>% filter(chet_diffabund == 'TRUE')
top_genera_het <- top_genera_trc %>% filter(het_diffabund == 'TRUE')
top_genera_chet_not <- top_genera_trc %>% filter(chet_diffabund == 'FALSE')


phylum <- microbiome::aggregate_taxa(CA_marine_ps, "Phylum")
sort(colSums(otu_table(phylum)["Desulfobacterota",])*100)

#Get total relative abundances of differentially abundant ASVs

marine_asv_diffabund <- data.frame(otu_table(marine_phylo_reads)) %>% 
  mutate(Kingdom = rownames(marine_phylo_reads@otu_table)) %>% 
  filter(Kingdom %in% diffabund_asvs) %>% select(-Kingdom)

tot_reads <- sum(sample_sums(marine_phylo_reads))
diffabund_ASV <- data.frame(t(marine_asv_diffabund))
tot_reads_diffabund <- sum(colSums(diffabund_ASV))
(tot_reads_diffabund/tot_reads)*100 #differentially abundant ASVs have a cumulative relative abundance of 31.6% across all samples. 

#MaAsLin2 input
asvs <- rownames(CA_marine_ps@otu_table)                        
CA_marine_ps_sub <- mutate_tax_table(CA_marine_ps, Species = asvs)
high_abund_diffabund_asvs <- colnames(diffabund_ASV)                        
marine_top_diffabund <- phyloseq::subset_taxa(CA_marine_ps_sub, 
                                              Species %in% high_abund_diffabund_asvs)
                        
diffabund_asv_tab <- data.frame(t(marine_top_diffabund@otu_table))
diffabund_tax_tab <- data.frame(marine_top_diffabund@tax_table)
diffabund_tax_tab$combined <- paste(diffabund_tax_tab$Species, diffabund_tax_tab$Genus, 
                                    sep = "; ") 
colnames(diffabund_asv_tab) <- diffabund_tax_tab$combined
diffabund_samp_tab <- data.frame(marine_top_diffabund@sam_data)
                        
library(Maaslin2)

diffabund_samp_tab <- diffabund_samp_tab %>% 
  mutate(B1_binned = as.factor(quantcut(B1, q = 4, 
                                           labels = c("low", "normal", "normal", "high")))) %>% 
  mutate(HMP_binned = as.factor(quantcut(HMP, q = 4, 
                                           labels = c("low", "normal", "normal", "high")))) %>%
  mutate(cHET_binned = as.factor(quantcut(cHET, q = 4, 
                                           labels = c("low", "normal", "normal", "high")))) %>% 
  select(B1, B1_binned, HMP, HMP_binned, cHET, cHET_binned, date, depth_cat_1, upwell_strength)
  
diffabund_samp_tab$B1_binned <- factor(diffabund_samp_tab$B1_binned, 
                                       levels = c("normal", "low", "high"))
diffabund_samp_tab$HMP_binned <- factor(diffabund_samp_tab$HMP_binned, 
                                       levels = c("normal", "low", "high"))
diffabund_samp_tab$cHET_binned <- factor(diffabund_samp_tab$cHET_binned, 
                                       levels = c("normal", "low", "high"))
diffabund_samp_tab$upwell_strength <- factor(diffabund_samp_tab$upwell_strength, 
                                       levels = c("intermediate", "weak", "strong"))

diffabund_samp_tab$depth_cat_1 <- factor(diffabund_samp_tab$depth_cat_1, 
                                             levels = c("surface", "intermediate", "deep"))

Maaslin2(input_metadata = diffabund_samp_tab, input_data = diffabund_asv_tab, 
                             min_prevalence = 0, normalization = "NONE", transform = "NONE", min_abundance = 0.001,
         output = "RL2103_Maaslin2_results_hmp_chet", 
         random_effects = c("date"), 
         fixed_effects = c("HMP_binned", "cHET_binned", "depth_cat_1"), 
         reference = c("HMP_binned,normal", "cHET_binned,normal", "depth_cat_1,surface"))

Maaslin2(input_metadata = diffabund_samp_tab, input_data = diffabund_asv_tab, 
         min_prevalence = 0, normalization = "NONE", transform = "NONE", min_abundance = 0.001,
         output = "RL2103_Maaslin2_results_b1_het_ammp", 
         random_effects = c("date"), 
         fixed_effects = c("B1_binned", "HET_binned", "AmMP_binned", "depth_cat_1"), 
         reference = c("B1_binned,normal", "HET_binned,normal", 
                       "AmMP_binned,normal", "depth_cat_1,surface"))

Maaslin2(input_metadata = diffabund_samp_tab, input_data = diffabund_asv_tab, 
         min_prevalence = 0, normalization = "NONE", transform = "NONE", min_abundance = 0.001,
         output = "RL2103_Maaslin2_results_upwelling", 
         random_effects = c("date"), 
         fixed_effects = "upwell_strength", 
         reference = "upwell_strength,intermediate")


###Section Six: Oligotyping
phylo_all <- readRDS("phylo_for_oligotyping")
phylo_sar11 <- subset_taxa(phylo_all, Order == "SAR11 clade")
phylo_sar11 <- merge_samples(phylo_sar11,"station")
otu_table(phylo_sar11) <- t(otu_table(phylo_sar11))
phylo_sar11 <- microbiomeMarker::normalize(phylo_sar11, "TSS") 
sar11_taxa <- data.frame(tax_table(phylo_sar11))
sar11_cladeIa <- sar11_taxa %>% filter(Genus == "Clade Ia")
cladeIa_seqs <- as.character(sar11_clade1a$Species)
sar11_asvs <- rownames(data.frame(otu_table(phylo_sar11)))

seqs <- DNAStringSet(cladeIa_seqs, use.names = TRUE)
names(seqs) <- rownames(sar11_clade1a)

sar11_aln <- msa::msa(seqs, method = "Muscle") #Need to export this as an alignment fasta file then import
alignment2Fasta(sar11_aln, filename = "path_to_/sar11_aln.fa")
OT.seq.concat1 <- MED('~/OneDrive/RStudio/sar11_aln.fa', Plot = TRUE, minseq = 21, entropymin = 0.6)
ENV <- GetEnvironmentDatafromFileR("~/OneDrive/RStudio/sar11_aln.fa",Start=2,Stop=9,test=FALSE)
Table0 <- SampleXOT_Table(OT.seq.concat = OT.seq.cconcat1, ENV = ENV, mosaicPlot = TRUE, filterByMinAbund = 0)   

#Shows entropy positions at bp, 67 and 105 (previously used by Bolanos et al., 2022 to define 1a.1 and 1a.3)
#Entropy position at the end of the reads likely from MiSeq error rate increasing

sar11_Ia <- subset_taxa(phylo_sar11_noseqs, Genus == "Clade Ia")

sar11_Ia_taxa <- c("Clade Ia.3", "Clade Ia.1", "Clade Ia.3", "Clade Ia.3", "Clade Ia.1", 
                   "Clade Ia.1", "Clade Ia.1", "Clade Ia.3", "Clade Ia.1", "Clade Ia.1", 
                   "Clade Ia.1", "Clade Ia.1", "Clade Ia.1", "Clade Ia.1", "Clade Ia.1", 
                   "Clade Ia.1", "Clade Ia.3", "Clade Ia.1", "Clade Ia.3")

sar11_Ia <- sar11_Ia %>% mutate_tax_table(Genus = sar11_Ia_taxa)
sar11_nonIa <- subset_taxa(phylo_sar11, Genus != "Clade Ia")
phylo_sar11 <- merge_phyloseq(sar11_Ia, sar11_nonIa) #Now they have clade Ia ecotype names
phylo_sar11 <- phylo_sar11 %>% mutate_tax_table(Species = sar11_asvs)

station_order <- c("NAV_466", "NAV_465", "NAV_464", "NAV_463", "NAV_461", 
                   "FAR_154", "FAR_152", "FAR_138", "FAR_139", "FAR_237", 
                   "MON_212", "MON_211", "MON_110", "MON_113", "MON_116", "MON_114", 
                   "SUR_105", "SUR_104", "SUR_103", "SUR_101")

phylo_sar11 %>% 
  comp_barplot(tax_level = "Species", n_taxa = 19, sample_order = rev(station_order), 
               other_name = "Other Taxa") +
  labs(x = NULL, y = NULL) + theme(legend.position = "right", 
                                   axis.text.y = element_text(face = "bold", size = 12)) + 
  coord_flip()

###Section 7: Alpha Diversity###

CA_div <- estimate_richness(CA_marine_ps, measures = c("Shannon", "Simpson", "InvSimpson"))

CA_div$"depth" <- CA_meta$depth_cat_2
CA_div$"line" <- CA_meta$line
CA_div$depth <- factor(CA_div$depth, levels = c("surface", "mixed", "deep", "dcm"))
CA_div$line <- factor(CA_div$line, levels = c("NAV", "FAR", "MON", "SUR"))

x.trans.norm <- function(x.trans) {
  (x.trans-min(x.trans))/(max(x.trans)-min(x.trans))
}

CA_div <- CA_div %>% mutate(Shannon = x.trans.norm(Shannon))

ggplot(CA_div, aes(x=depth, y=Shannon,fill=depth)) +
  geom_boxplot()+
  geom_point(size = 1)+
  facet_wrap(~line,scales="free_y",nrow=1)+
  ggh4x::facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(3.2,5)), 
    scale_y_continuous(limits = c(3.2,5)),
    scale_y_continuous(limits = c(3.2,5)),
    scale_y_continuous(limits = c(3.2,5)))) +
theme_bw() +
  scale_fill_manual(values=depthpal, name = "Depth", 
                    labels = c("Surface", "Intermediate", "Deep", "DCM")) +
  labs(x = "", y = "Shannon Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = .9)) +
  ggpubr::stat_compare_means(Shannon~line, 
                             data = CA_div, 
                             method = "kruskal.test",
                             mapping = aes(label =  ..p.signif..)) +
  scale_x_discrete(labels = c("Surface", "Intermediate", "Deep", "DCM"))

ggplot(CA_div, aes(x=depth, y=Shannon,fill=depth)) +
  geom_boxplot()+
  geom_point(size = 2) +
  theme_bw() +
  scale_fill_manual(values=depthpal, name = "Depth", 
                    labels = c("Surface", "Intermediate", "Deep", "DCM")) +
  labs(x = "", y = "Shannon Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = .9)) +
  ggpubr::stat_compare_means(Shannon~depth, 
                             data = CA_div, 
                             method = "kruskal.test",
                             mapping = aes(label =  ..p.signif..)) +
  scale_x_discrete(labels = c("Surface", "Intermediate", "Deep", "DCM"))

ggplot(CA_div, aes(x=line, y=Shannon,fill=line)) +
  geom_boxplot()+
  geom_point(size = 2)+
  theme_bw() +
  scale_fill_manual(values=set1, name = "Line") +
  labs(x = "", y = "Shannon Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = .9)) +
  ggpubr::stat_compare_means(Shannon~line, 
                             data = CA_div, 
                             method = "kruskal.test",
                             mapping = aes(label =  ..p.signif..))

###Section 8: Mantel Tests and Stepwise AIC Regression###
CA_marine_ps_mixed <- subset_samples(marine_genus, depth_cat_1 == "surface" | depth_cat_1 == "intermediate") #sig
CA_marine_ps_mixed <- prune_taxa(taxa_sums(CA_marine_ps_mixed) > 0, CA_marine_ps_mixed) 

CA_marine_meta <- data.frame(sample_data(CA_marine_ps))
mixed_meta <- data.frame(sample_data(CA_marine_ps_mixed))

abund_mixed <- data.frame(t(otu_table(CA_marine_ps_mixed)))
abund <- data.frame(t(otu_table(marine_genus)))

abund_mixed_dist <- vegdist(abund_mixed, method = "bray")
vegdist_den <- vegdist(mixed_meta$density, method = "euclidean")
vegdist_sal <- vegdist(mixed_meta$sal, method = "euclidean")
vegdist_oxy <- vegdist(mixed_meta$oxygen, method = "euclidean")

set.seed(1032)
mantel(abund_mixed_dist, vegdist_den, method = "spearman", permutations = 9999, na.rm = TRUE) #Community dissimilarity in the mixed layer increases with density, i.e. under upwelling 
mantel(abund_mixed_dist, vegdist_sal, method = "spearman", permutations = 9999, na.rm = TRUE) #Same for salinity
mantel(abund_mixed_dist, vegdist_oxy, method = "spearman", permutations = 9999, na.rm = TRUE) #Same for oxygen

abund_dist <- vegdist(abund, method = "bray")

vegdist_depth <- vegdist(CA_marine_meta$actual.depth...30, method = "euclidean")
mantel(abund_dist, vegdist_depth, method = "spearman", permutations = 9999, na.rm = TRUE) #Community dissimilarity increases with depth

#Stepwise AIC model 

stepwise_data <- data.frame(CA_marine_ps@sam_data) %>% dplyr::select(CUTI, BEUTI, B1:HET, oxygen, density, chl)

stepwise_data <- stepwise_data %>% mutate(across(B1:chl, x.trans.norm)) %>% 
  mutate(across(B1:density, transformTukey))
stepwise_data$upwelling <- diffabund_samp_tab$upwell_strength

trc_mod <- lm(CUTI ~ B1+HMP+cHET+HET+AmMP++oxygen+density+B1:AmMP+B1:HET+
                B1:cHET+B1:HMP, 
   data = stepwise_data)

set.seed(1111)
trc_stepmod <- stepAIC(trc_mod, direction = "both")
summary(trc_stepmod) #chlorophyll not deemed an important variable

trc_mod_beut <- lm(BEUTI ~ B1+HMP+cHET+HET+AmMP+oxygen+density+B1:AmMP+B1:HET+HMP:AmMP+cHET:HET+
                B1:cHET+B1:HMP+HMP:cHET+HMP:HET+AmMP:HET+AmMP:cHET, 
              data = stepwise_data)
set.seed(1112)
trc_stepmod_beut <- stepAIC(trc_mod_beut, direction = "both")
summary(trc_stepmod_beut)











