library(data.table)
library(tidyverse)
library(vegan)
library(RcppCNPy)
library(windowscanr)
library(ggman)

setwd("/Volumes/Accelsior 4M2/Projects/Lumpfish/")
Lump_WGS_meta <- fread("Lump_WGS_HD_meta.txt")
OSC_2016 <- str_subset(Lump_WGS_meta$BAM, "OSC_2016")
OSC_2017 <- str_subset(Lump_WGS_meta$BAM, "OSC_2017")

Table1_Pub<- Lump_WGS_meta %>%   count(Pop, Region, Lat, Lon)  
write.table(Table1_Pub, "Table1_Pub.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
Lump_WGS_meta %>%  count(Region)
Lump_WGS_meta_map <- Lump_WGS_meta %>%  select(Pop, Lat, Lon, Region) %>%  distinct() %>%  filter(!Region %in% c("Europe", "OSC"))

#2 filter for PCANGSD plotting
Lump_WGS_meta$Filter <- 0
Lump_WGS_meta$Filter[!(Lump_WGS_meta$Region %in% c("Europe", "OSC"))] <- 1 
Lump_WGS_meta_NA  <-Lump_WGS_meta %>% filter(Filter == 1)

Lump_WGS_meta$Filter <- 0
Lump_WGS_meta$Filter[(Lump_WGS_meta$Region %in% c("Newfoundland"))] <- 1 
Lump_WGS_meta_NL  <-Lump_WGS_meta %>% filter(Filter == 1)

Lump_WGS_meta$Filter <- 0
Lump_WGS_meta$Filter[(Lump_WGS_meta$Region %in% c("Newfoundland", "Gulf of St. Lawrence"))] <- 1 
Lump_WGS_meta_NL  <-Lump_WGS_meta %>% filter(Filter == 1)

#1 map
library(marmap)
bathydata <- getNOAA.bathy(-71,-50,41,55, res=3,keep=T)
plot(bathydata)
map=autoplot(bathydata, geom=c("r", "c"), colour="grey", size=0.1) +
  scale_fill_etopo(guide = FALSE) +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "none")
map + geom_point(data =Lump_WGS_meta_map , aes(x = Lon, y = Lat, col = Region), inherit.aes = F) + theme(legend.position="right") 

#pop structure
covmat_NA <-  read.table("PCANGSD_NA/Cyclopterus_lumpus_sizechecked_NA.cov", quote="\"",
                         comment.char="", stringsAsFactors=FALSE)
PCA_NA <- eigen(covmat_NA)

#var explained per axis
eigenvalues <- PCA_NA$values
(eigenvalues[1]/sum(eigenvalues))*100
(eigenvalues[2]/sum(eigenvalues))*100
(eigenvalues[3]/sum(eigenvalues))*100

PCAscores_NA <- data.frame(cbind(Lump_WGS_meta_NA,
                                 PCA_NA$vectors[,1:10]), stringsAsFactors = F)

ggplot() + geom_point(data = PCAscores_NA, aes(x = V1, y = V2, colour=Region)) + 
  theme_classic() + 
  geom_hline(yintercept = -0.03) +
  geom_vline(xintercept = -0.05)


JCLUSTER_NA <- PCAscores_NA %>%  filter(V1 > -0.05, V2 < -0.035)


ggplot() + geom_point(data = PCAscores_NA, aes(x = V1, y = V2, colour=Pop)) + theme_classic()


#Admix NA
NA_admix <- fread("PCANGSD_NA/Cyclopterus_lumpus_sizechecked_NA.admix.3.Q")
Lump_admix_NA  <- data.frame(cbind(Lump_WGS_meta_NA, NA_admix))

Admix_table_NA <- Lump_admix_NA
rownames(Admix_table_NA ) <- Lump_admix_NA$BAM

plot_data_NA <-  Admix_table_NA  %>% 
  gather('pop', 'prob', V1:V3) %>% 
  group_by(Region)

ggplot(plot_data_NA, aes(BAM, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Region, scales = 'free', space = 'free')


#now all NL
##PCA NL
covmat_NL <-  read.table("PCANGSD_NL/Cyclopterus_lumpus_sizechecked_NL.cov", quote="\"",
                         comment.char="", stringsAsFactors=FALSE)
PCA_NL <- eigen(covmat_NL)

#var explained per axis
eigenvalues <- PCA_NL$values
(eigenvalues[1]/sum(eigenvalues))*100
(eigenvalues[2]/sum(eigenvalues))*100
(eigenvalues[3]/sum(eigenvalues))*100

PCAscores_NL <- data.frame(cbind(Lump_WGS_meta_NL,
                                 PCA_NL$vectors[,1:10]), stringsAsFactors = F)

ggplot() + geom_point(data = PCAscores_NL, aes(x = V1, y = V2, colour=Region)) + theme_classic()
ggplot() + geom_point(data = PCAscores_NL, aes(x = V1, y = V2, colour=Pop)) + theme_classic() + geom_vline(xintercept = -0.03)


#Admix NL
NL_admix <- fread("PCANGSD_NL/Cyclopterus_lumpus_sizechecked_NL.admix.2.Q")
Lump_admix_NL  <- data.frame(cbind(Lump_WGS_meta_NL, NL_admix))

Admix_table_NL <- Lump_admix_NL
rownames(Admix_table_NL ) <- Lump_admix_NL$BAM

plot_data_NL <-  Admix_table_NL  %>% 
  gather('pop', 'prob', V1:V2) %>% 
  group_by(Pop)

ggplot(plot_data_NL, aes(BAM, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Pop, scales = 'free', space = 'free')

ggplot(plot_data_NL, aes(BAM, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Region, scales = 'free', space = 'free')

#for North South comparisons
NORTH_samps_FST <- PCAscores_NA %>%  filter(Region %in% c("Newfoundland", "Gulf of St. Lawrence")) %>% select(BAM) %>%  sample_n(64)
write.table(NORTH_samps_FST, "NORTH_samps_FST.tsv", col.names = F, row.names = F, sep = "\t", quote = F)


SOUTH_samps_FST <- PCAscores_NA %>%  filter(Region %in% c("Maritimes", "Gulf of Maine")) %>%  select(BAM)
write.table(SOUTH_samps_FST , "SOUTH_samps_FST.tsv", col.names = F, row.names = F, sep = "\t", quote = F)


#for NL Cluster 1 Cluster 2 comparisons

JCLUSTER_FST <- JCLUSTER_NA %>%  slice_min(order_by = V2, n = 50) %>%  
  select(BAM) 
write.table(JCLUSTER_FST, "JCLUSTER_samps_FST_50.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

Newfoundland_samps <- PCAscores_NA %>%  filter(Region %in% "Newfoundland", !BAM %in% JCLUSTER_NA$BAM) %>% select(BAM)
Newfoundland_samps_FST <- Newfoundland_samps %>%  sample_n(50)
write.table(Newfoundland_samps_FST, "Newfoundland_samps_FST_50.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

##FST comparisons
system("vcftools --gzvcf Cyclopterus_lumpus_sizechecked_phased_chromchange.vcf.gz --maf 0.05 --weir-fst-pop NORTH_samps_FST.tsv --weir-fst-pop SOUTH_samps_FST.tsv --out NORTH_SOUTH  ")

system("vcftools --gzvcf  Cyclopterus_lumpus_sizechecked_phased_chromchange.vcf.gz --maf 0.05 --weir-fst-pop Newfoundland_samps_FST_50.tsv --weir-fst-pop JCLUSTER_samps_FST_50.tsv --out JCluster_Newfoundland_50  ")

#NA_PCLoadings
#import pcangsd fastPCA loadings per site
NA_PCLoadings <- data.frame(npyLoad("PCANGSD_NA/Cyclopterus_lumpus_sizechecked_NA.selection.npy")) 
#import pcangsd sites output
NA_sites <- fread("PCANGSD_NA/Cyclopterus_lumpus_sizechecked_NA.sites", data.table = F,
                  stringsAsFactors = F, header = F)

Loc_pos <- fread("Marker_order_lumpfishbeag")
setnames(Loc_pos, new = c("CHROM", "POS"))

NA_Loc_pos <- bind_cols(Loc_pos, NA_sites)
NA_Loc_pos <- NA_Loc_pos %>% filter(V1 %in% "1")

NA_PC_pos <- bind_cols(NA_Loc_pos , NA_PCLoadings)

colnames(NA_PC_pos)[3:5] <- c("KEEP", "PC1", "PC2")
NA_PC_pos$pvals_PC1 <- 1- pchisq(q=NA_PC_pos$PC1, df = 1)

NA_PC_pos <- NA_PC_pos %>%  mutate(SNP = paste0(CHROM, "_", POS))

ChromConversion <- fread("ChromosomeConversion.txt") %>% mutate(CHROM= SequenceName) %>%  select(-SequenceName)

NA_PC_pos <- inner_join(ChromConversion, NA_PC_pos)


#this can be made nicer...
PCA_10KB_WINDO <- function(frame){frame_maxes <-frame  %>% 
  group_by(CHROM) %>%
  filter(POS == max(POS))
frame_maxes <- frame_maxes %>%  mutate(POS = POS + 2000)
frame_maxed <- rbind(frame,frame_maxes)
PC_10kwin <- winScan(x = frame_maxed, groups  = "AssignedMoleculeChromosome", 
                    position = "POS", 
                    win_size = 10000, 
                    win_step = 10000,  
                    values = "PC1", 
                    funs = "mean",
                    cores = 28)
PC_10kwin  <- PC_10kwin %>%  filter(PC1_n > 0) %>% drop_na()
return(PC_10kwin) }

#10KB window
PC1_NA_10K <- PCA_10KB_WINDO(NA_PC_pos)
PC1_NA_10K <- PC1_NA_10K %>% mutate(POS_ID = paste0(AssignedMoleculeChromosome, "_", win_mid))
PC1_NA_10K_99_OL<- PC1_NA_10K  %>%  filter(PC1_mean > quantile(PC1_mean, 0.99))
min(PC1_NA_10K_99_OL$PC1_mean)
ggman(gwas = PC1_NA_10K , snp = "POS_ID", bp = "win_start", chrom = "AssignedMoleculeChromosome", pvalue = "PC1_mean", logTransform = F, pointSize = 1, ymax = 150, ymin = 0, sigLine = 9.69, ylabel = "PC1 NA") + theme_classic()


#NL_PCLoadings
#import pcangsd fastPCA loadings per site
NL_PCLoadings <- data.frame(npyLoad("PCANGSD_NL/Cyclopterus_lumpus_sizechecked_NL.selection.npy")) 
#import pcangsd sites output
NL_sites <- fread("PCANGSD_NL/Cyclopterus_lumpus_sizechecked_NL.sites", data.table = F,
                  stringsAsFactors = F, header = F)

Loc_pos <- fread("Marker_order_lumpfishbeag")
setnames(Loc_pos, new = c("CHROM", "POS"))

NL_Loc_pos <- bind_cols(Loc_pos, NL_sites)
NL_Loc_pos <- NL_Loc_pos %>% filter(V1 %in% "1")

NL_PC_pos <- bind_cols(NL_Loc_pos , NL_PCLoadings)

colnames(NL_PC_pos)[3:4] <- c("KEEP", "PC1")
NL_PC_pos$pvals_PC1 <- 1- pchisq(q=NL_PC_pos$PC1, df = 1)

NL_PC_pos <- NL_PC_pos %>%  mutate(SNP = paste0(CHROM, "_", POS))

NL_PC_pos <- inner_join(ChromConversion, NL_PC_pos)

PC1_NL_10K <- PCA_10KB_WINDO(NL_PC_pos)
PC1_NL_10K <- PC1_NL_10K %>% mutate(POS_ID = paste0(AssignedMoleculeChromosome, "_", win_mid))
PC1_NL_10K_99_OL<- PC1_NL_10K  %>%  filter(PC1_mean > quantile(PC1_mean, 0.99))
min(PC1_NL_10K_99_OL$PC1_mean)
ggman(gwas = PC1_NL_10K , snp = "POS_ID", bp = "win_start", chrom = "AssignedMoleculeChromosome", pvalue = "PC1_mean", logTransform = F, pointSize = 1, ymax =80, ymin = 0, sigLine =7.65, ylabel = "PC1 NL") + theme_classic()





#and compare with FSTs
#FSTs
FST_10KB_WINDO <- function(frame){frame_maxes <-frame  %>% 
  group_by(CHROM) %>%
  filter(POS == max(POS))
frame_maxes <- frame_maxes %>%  mutate(POS = POS + 2000)
frame_maxed <- rbind(frame,frame_maxes)
FST_10kwin <- winScan(x = frame_maxed, groups  = "AssignedMoleculeChromosome", 
                     position = "POS", 
                     win_size = 10000, 
                     win_step = 10000,  
                     values = "WEIR_AND_COCKERHAM_FST", 
                     funs = "mean",
                     cores = 28)
FST_10kwin  <- FST_10kwin %>%  filter(WEIR_AND_COCKERHAM_FST_n > 0) %>% drop_na()
return(FST_10kwin) }

NORTH_SOUTH_50.fst <- fread("NORTH_SOUTH.weir.fst") %>% 
  mutate(CHROM = paste0("S", CHROM)) %>%  
  mutate(SNP = paste0(CHROM, "_", POS)) %>%  
  inner_join(ChromConversion)
##Weir and Cockerham mean Fst estimate: 0.013997
###Weir and Cockerham weighted Fst estimate: 0.017387

NS_10KB_FST <- FST_10KB_WINDO(NORTH_SOUTH_50.fst)
NS_10KB_FST_99_OL <- NS_10KB_FST %>% filter(WEIR_AND_COCKERHAM_FST_mean > quantile(WEIR_AND_COCKERHAM_FST_mean, 0.99))


#get VCF PC1 overlap North South
NORTH_SOUTH_PC_FST_OL <- inner_join(NS_10KB_FST_99_OL, PC1_NA_10K_99_OL)
NORTH_SOUTH_PC_FST_OL_bed <- inner_join(NORTH_SOUTH_PC_FST_OL, ChromConversion) %>% 
  select(CHROM, win_start, win_end, PC1_mean, WEIR_AND_COCKERHAM_FST_mean ) %>%  mutate(win_end = win_end+1)  ##add 1 to avoid sci notations that make bedtools explode
write.table(NORTH_SOUTH_PC_FST_OL_bed,
            "NORTH_SOUTH_PC1_FST_99OL.bed",
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)
system("bedtools intersect -a Lump_genes_entrez.bed -b NORTH_SOUTH_PC1_FST_99OL.bed | cut -f4 | sort | uniq > NORTH_SOUTH_PC1_FST_99OL_entrezIDS.tsv")
system("bedtools intersect -a Lump_genes_entrez.bed -b NORTH_SOUTH_PC1_FST_99OL.bed -wb > NORTH_SOUTH_PC1_FST_99OL_entrez_ALL.tsv ")
NORTH_SOUTH_PC1_FST_OL_entrez <- fread("NORTH_SOUTH_PC1_FST_99OL_entrez_ALL.tsv") %>%  mutate(CHROM = V1,
                                                                                       win_start = V2 -1, 
                                                                                       win_end = V3, 
                                                                                       PC1_mean = V8, 
                                                                                       entrezID = V4, 
                                                                                       FST_mean = V9) %>% 
  select(CHROM, win_start, win_end,  PC1_mean, FST_mean, entrezID) %>% 
  inner_join(ChromConversion)

NORTH_SOUTH_PC1_FST_OL_entrez<- NORTH_SOUTH_PC1_FST_OL_entrez %>% arrange(desc(PC1_mean))
write.table(NORTH_SOUTH_PC1_FST_OL_entrez, "NORTH_SOUTH_PC1_FST_OL_entrez.tsv", col.names = T, 
            row.names = F, 
            sep =  "\t", 
            quote = F)


## Newfoundland comparison now
JCluster_Newfoundland_50.fst <- fread("JCluster_Newfoundland_50.weir.fst")  %>% 
  mutate(CHROM = paste0("S", CHROM)) %>%  
  mutate(SNP = paste0(CHROM, "_", POS)) %>%  
  inner_join(ChromConversion)
###Weir and Cockerham mean Fst estimate: 0.00927
###Weir and Cockerham weighted Fst estimate: 0.011612


NL12_10KB_FST <- FST_10KB_WINDO(JCluster_Newfoundland_50.fst)
NL12_10KB_FST_99_OL <- NL12_10KB_FST %>% filter(WEIR_AND_COCKERHAM_FST_mean > quantile(WEIR_AND_COCKERHAM_FST_mean, 0.99))

#get VCF PC1 overlap North South
NL12_PC_FST_OL <- inner_join(NL12_10KB_FST_99_OL, PC1_NL_10K_99_OL)
NL12_PC_FST_OL_bed <- inner_join(NL12_PC_FST_OL, ChromConversion) %>% 
  select(CHROM, win_start, win_end, PC1_mean, WEIR_AND_COCKERHAM_FST_mean ) %>%  mutate(win_end = win_end+1)  ##add 1 to avoid sci notations that make bedtools explode
write.table(NL12_PC_FST_OL_bed,
            "NL12_PC1_FST_99OL.bed",
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)
system("bedtools intersect -a Lump_genes_entrez.bed -b NL12_PC1_FST_99OL.bed | cut -f4 | sort | uniq > NL12_PC1_FST_99OL_entrezIDS.tsv")
system("bedtools intersect -a Lump_genes_entrez.bed -b NL12_PC1_FST_99OL.bed -wb > NL12_PC1_FST_99OL_entrez_ALL.tsv ")
NL12_PC1_FST_OL_entrez <- fread("NL12_PC1_FST_99OL_entrez_ALL.tsv") %>%  mutate(CHROM = V1,
                                                                                              win_start = V2 -1, 
                                                                                              win_end = V3, 
                                                                                              PC1_mean = V8, 
                                                                                              entrezID = V4, 
                                                                                              FST_mean = V9) %>% 
  select(CHROM, win_start, win_end,  PC1_mean, FST_mean, entrezID) %>% 
  inner_join(ChromConversion)

NL12_PC1_FST_OL_entrez<- NL12_PC1_FST_OL_entrez %>% arrange(desc(PC1_mean))
write.table(NL12_PC1_FST_OL_entrez, "NL12_PC1_FST_OL_entrez.tsv", col.names = T, 
            row.names = F, 
            sep =  "\t", 
            quote = F)




inner_join(NL12_PC1_FST_OL_entrez, NORTH_SOUTH_PC1_FST_OL_entrez, by = "entrezID")
write.table(inner_join(NL12_PC1_FST_OL_entrez, NORTH_SOUTH_PC1_FST_OL_entrez, by = "entrezID"), "NS_NL12_FSTPCOL_entrezoverlap.tsv", col.names = T, 
            row.names = F, quote = F, sep = "\t")

