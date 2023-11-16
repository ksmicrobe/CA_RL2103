library(dada2)
library(phyloseq)
library(tidyverse)
library(decontam)
library(lessR)

path <- "path_to_fastq_files"

trimgalore_to_seqtab <- function(path){
  fnFs <- sort(list.files(path, pattern="_R1_001_trimmed.fq.gz", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="_R2_001_trimmed.fq.gz", full.names = TRUE))
  print("files read in")
  sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  filtFs <- file.path(path, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample_names
  print(sample_names)
  names(filtFs) <- sample_names
  names(filtRs) <- sample_names
  out <<- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        maxN=0, truncQ=0, rm.phix=TRUE,
                        compress=FALSE, multithread=TRUE, matchIDs = TRUE, 
                        maxEE = 2, trimLeft = c(8,120)) #parameters picked based on multiqc results
  print("sequences trimmed")
  derepF <- derepFastq(filtFs)
  derepR <- derepFastq(filtRs)
  print("reads dereplicated")
  errF <<- learnErrors(derepF, multithread=TRUE)
  errR <<- learnErrors(derepR, multithread=TRUE)
  print("error rates learned")
  dadaFs <<- dada(derepF, err=errF, multithread=TRUE, pool = "pseudo")
  dadaRs <<- dada(derepR, err=errR, multithread=TRUE, pool = "pseudo")
  print("ASVs assigned")
  mergers <<- mergePairs(dadaFs, derepF, dadaRs, derepR, verbose=TRUE)
  print("reads merged")
  seqtab <<- makeSequenceTable(mergers)
  seqtab_nochim <<- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  print(table(nchar(getSequences(seqtab_nochim))))
  getN <- function(x) sum(getUniques(x))
  track <<- cbind(out, sapply(dadaFs, getN), 
                  sapply(dadaRs, getN), 
                  sapply(mergers, getN), 
                  rowSums(seqtab_nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample_names
  View(track)
}

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)

mont_seqtab_nocontam_nochim <- readRDS("MONT_WC_seqtab_nocontam_nochim") #Previously-processed data from Monterey Bay samples (DADA2 default parameters, w/ pseudo-pooling)
fresh_WC_nocontam_nochim <- as.data.frame(t(readRDS("SA_WC_nocontam_nochim"))) #Fresh water data from Sacramento River

silva <- "path_to_silva_db"

neg_ctrl_seqtab <- seqtab_nochim_marine[c(54,55),]

nc_taxa <- assignTaxonomy(neg_ctrl_seqtab, refFasta = silva)
sub_marine_asvs <- lessR::to("ASV", nrow(seqtab_nochim_marine))

dim(seqtab_nochim_marine)
dim(mont_seqtab_nocontam_nochim)
dim(fresh_WC_nocontam_nochim)

seqtab1 <- as.matrix(seqtab_nochim)
seqtab2 <- as.matrix(mont_seqtab_nocontam_nochim)
seqtab3 <- as.matrix(fresh_WC_nocontam_nochim)

fresh_marine_seqtab_nochim <- mergeSequenceTables(table1 = seqtab1, 
                                                  table2 = seqtab2, 
                                                  table3 = seqtab3)

total_sample_names <- c("NAV463_3", "NAV463_15", "NAV463_32", "NAV463_100", "NAV461_3", 
                        "NAV461_10", "NAV461_30", "NAV461_65", "NAV466_3", "NAV466_15", 
                        "NAV466_40", "NAV466_100", "NAV465_3", "NAV465_10", "NAV465_100", 
                        "NAV464_3", "NAV464_10", "NAV464_20", "NAV464_100", "FAR154_3", 
                        "FAR154_11", "FAR154_40", "FAR154_100", "FAR152_3", "FAR152_12", 
                        "FAR152_30", "FAR152_95", "FAR138_10", "FAR138_25", "FAR138_40", 
                        "SVR105_3", "SVR105_12", "SVR105_50", "SVR105_100", "SVR104_3", 
                        "SVR104_10", "SVR104_25", "SVR104_100", "SVR103_3", "SVR103_10", 
                        "SVR103_25", "SVR103_85", "SVR101_3", "SVR101_20", "SVR101_40", 
                        "SVR101_65", "FAR237_3", "FAR237_8", "FAR237_15", "FAR237_50", 
                        "FAR139_3", "FAR139_6", "FAR139_12", "FAR139_40", "KCS_NC_water", 
                        "KCS_NC_sterivex", "MON110_100", "MON110_40", "MON116_15", "MON110_3", 
                        "MON113_100", "MON113_3", "MON116_3", "MON113_40", "MON212_40", 
                        "MON110_20", "MON113_20", "MON114_12", "MON116_110", "MON114_65", 
                        "MON211_40", "MON211_25", "MON114_3", "MON212_20", "MON211_100", 
                        "MON116_30", "MON211_3", "MON114_20", "MON212_100", "MON212_3", 
                        "BCC_3", "BCC_1", "BCC_4", "BCC_2", "CCH_3", 
                        "CCH_1", "CCH_4", "CCH_2", "CCM_3", "CCM_1", 
                        "CCM_4", "CCM_2", "FRD_3", "KLA_3", "KLA_2", 
                        "NFB_3", "NFB_1", "NFC_4", "NFC_2", "SAC_3", 
                        "SAC_1", "SAC_4", "SAC_2", "SFB_3", "SFB_1", 
                        "SFB_4", "SFB_2", "SSC_3", "SSC_1", "SSC_4", 
                        "SSC_2", "FRU_1", "FRD_1", "FRU_2", "FRD_2", 
                        "FRU_3", "FR3_3", "FR2_3", "FRU_4", "FRD_4")

tot_sample_type <- c("marine", "marine", "marine", "marine", "marine", 
                     "marine", "marine", "marine", "marine", "marine", 
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "marine", "marine", "marine", "marine", "marine",
                     "river", "river", "river", "river", "river", 
                     "river", "river", "river", "river", "river",
                     "river", "river", "river", "river", "river",
                     "river", "river", "river", "river", "river",
                     "river", "river", "river", "river", "river",
                     "river", "river", "river", "river", "river",
                     "river", "river", "river", "river", "river",
                     "river", "river", "river", "river", "river")

tot_meta <- data.frame(sample = total_sample_names, type = tot_sample_type)
rownames(tot_meta) <- tot_meta$sample
rownames(fresh_marine_seqtab_nochim) <- tot_meta$sample

fresh_marine_taxa <- assignTaxonomy(fresh_marine_seqtab_nochim, refFasta = silva, 
                                    verbose = TRUE)

phytoref <- "path_to_phytoref_db"


fresh_marine_taxa$"sequence" <- rownames(fresh_marine_taxa)
fresh_marine_asvs <- lessR::to("ASV", nrow(fresh_marine_taxa))
rownames(fresh_marine_taxa) <- fresh_marine_asvs

algal_taxa <- fresh_marine_taxa %>% 
  filter(Order == "Chloroplast") %>% 
  select(sequence)

algal_taxa_phytoref <- assignTaxonomy(algal_taxa$sequence, refFasta = phytoref, 
                                      taxLevels = algal_levels)

algal_levels <- c("Domain", "Sub-Kingdom", "Phylum", "Class", "Sub-Class", 
                  "Order", "Sub-Order", "Family", "Genus", "Species")

phytoref_taxa <- as.data.frame(algal_taxa_phytoref)
phytoref_taxa$"sequence" <- rownames(algal_taxa_phytoref)

raw_counts <- data.frame(t(fresh_marine_seqtab_nochim))
raw_counts$"sequence" <- rownames(raw_counts)
rownames(raw_counts) <- fresh_marine_asvs
algal_taxa_counts <- fresh_marine_taxa %>% 
  left_join(raw_counts, by = "sequence") %>% 
  filter(Order == "Chloroplast") %>% 
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>% 
  left_join(phytoref_taxa, by = "sequence") %>% 
  select(-'Sub-Kingdom', -'Sub-Class', -'Sub-Order', -Species) %>% 
  dplyr::rename(Kingdom = Domain)

tot_taxa_counts <- fresh_marine_taxa %>% 
  left_join(raw_counts, by = "sequence") %>% 
  dplyr::filter(!Order %in% "Chloroplast", 
                !Family %in% "Mitochondria", 
                !Kingdom %in% "Eukaryota")

tot_taxa_counts <- rbind(tot_taxa_counts, algal_taxa_counts)

dim(tot_taxa_counts)
dim(algal_taxa_counts) #Sanity check

fresh_marine_asvs <- lessR::to("ASV", nrow(tot_taxa_counts))
rownames(tot_taxa_counts) <- fresh_marine_asvs

###phytoref-assigned ASVs finalized and all non-algal Eukaryotes are removed
#now make into phyloseq object and run through decontam

fresh_marine_counts <- tot_taxa_counts %>% 
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -sequence)
fresh_marine_tax <- tot_taxa_counts %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus)

ps_contam <- phyloseq(otu_table(as.matrix(fresh_marine_counts), taxa_are_rows = TRUE), 
                      tax_table(as.matrix(fresh_marine_tax)), 
                      sample_data(tot_meta))


sample_or_control <- c("sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "control", 
                       "control", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample", 
                       "sample", "sample", "sample", "sample", "sample")

sample_data(ps_contam)$"sample_or_control" <- sample_or_control
sample_data(ps_contam)$is.neg <- sample_data(ps_contam)$sample_or_control == "control"
contamdf_prev <- isContaminant(ps_contam, method = "prevalence", neg = "is.neg", threshold=0.5)
table(contamdf_prev$contaminant)

bad_asvs <- c("ASV01906", "ASV02322", "ASV05692", "ASV08104", "ASV10949", 
              "ASV11590", "ASV17980", "ASV17990", "18024")

ps.pa <- transform_sample_counts(ps_contam, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "sample", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    +                     contaminant=contamdf.prev$contaminant)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
                      +     xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
                    
                    ##based on the plot, it appears "contaminant" asvs are found both in positive and negative samples, 
                    #they are all marine-associated taxa, and the negative controls for this project displayed no gel bands. 
                    #Therefore I will not remove any of the decontam-identified ASVs because they are likely true ASVs that were cross contaminated into negative control samples. 
                    
                    fresh_marine_counts <- fresh_marine_counts[,-c(55,56)]
                    fresh_marine_counts <- fresh_marine_counts[,-119]
                    
                    fresh_marine_counts <- fresh_marine_counts %>% 
                      dplyr::rename(MON113_3 = MON133_3)
                    
                    ###
                    ps_marine_fresh <- phyloseq(otu_table(as.matrix(fresh_marine_counts), taxa_are_rows = TRUE), 
                                                tax_table(as.matrix(fresh_marine_tax)), 
                                                sample_data(tot_meta))
                    
                    ps_marine_fresh <- prune_taxa(taxa_sums(ps_marine_fresh) > 0, ps_marine_fresh)
                    ###
                    
                    
                    
                    saveRDS(ps_marine_fresh, file = "ps_marine_fresh_unrarefied")
                    
                    
                    
                    
                    
                    
                    