## -------------------------------------------------------------------------
##
## Script name: 20240308_R_script_blast_results_heatmap_PCA.R
##
## Purpose: short-BLASTN searches and visualization of results as homology-based heatmap for SARS-CoV-2 and SARS-CoV-1 genomic segments 
##
## Authors: Amber R. Paulson, Vincent Montoya
##
## Date Created: 2024-03-08
##
##
## Copyright (c) Amber Paulson, 2025 - governed by CC-BY-NC 4.0 - Attribution-Non-Commercial 4.0 International 
## see https://creativecommons.org/licenses/by-nc/4.0/
##
## Email: Amber[dot]Rose[dot]Paulson[at]gmail[dot]com
## Github: https://github.com/damselflywingz/EpiVIRTEX/
##
##
## ------------------------------------------------------------------------
## Notes:
##   
## ### 2025-08-09: updates made to Palm Civet genome for the correct species 
##
##
##
## ------------------------------------------------------------------------

# !Important Note on implementation: 
# The R Project should be opened with RStudio first, and then next proceed to open and run the script.


##-------------------------------------------------------------------------
##              Virus and host reference genomes                          -
##-------------------------------------------------------------------------

# Reference genomes downloaded and short-BLASTN performed in Linux environment
# example of working directory on server is: ~/working_dir/Local_blast/20240301_SCV2_vs_R.ferr/batch_blast_test 
# SARS-CoV-2 Wuhan-Hu-1 (NC_045512.2) and SARS-CoV-1 Tor2 (AY274119.3) were accessed from NCBI and sequenced split into sequence of 25 nt length

# October 27, 2023 running some local blasts - examples only shown

# downloaded genomes from NCBI Genbank, some from the Bloom paper regarding animals at the wet market
# noted Rhizomys pruinosus, Erinaceau amurensis, and Marmota himalayana were not of good quality or not available as a genomes
# Masked palm civet gene from CNGBdb

# Used Rhinophus ferrumequinum (mRhiFer1_v1), the RatG13 host (horseshoe bat) has a poor genome assembly
# DSBC_Hbra_1.0 = Hystrix brachyura
# mManPen7.hap1 = Manis pentadactyla
# NYPRO = Nyctereutes procyonoides
# Plar_chr = Paguma larvata

# control hosts
# bGalGal1.mat.broiler.GRCg7b = Gallus gallus (Chicken)
# rDiaPun1.0.p = Diadophis punctatus similis (San Diego ring-necked snake)

# blast is on the server version 2.6.0+

# bat genome
# header has spaces and ,
cat GCF_004115265.2_mRhiFer1_v1.p_genomic.fna | tr " " _ > GCF_004115265.2_mRhiFer1_v1.p_genomic_.fna
cat GCF_004115265.2_mRhiFer1_v1.p_genomic_.fna | tr "," _ > GCF_004115265.2_mRhiFer1_v1.p_genomic_.fna

makeblastdb -dbtype nucl -parse_seqids -in GCF_004115265.2_mRhiFer1_v1.p_genomic__.fna -out R_ferrumequinum_bat -title R_ferrumequinum_genome

# August 9, 2024 - updating for the masked palm civet from Chinese Genome 
makeblastdb -dbtype nucl -parse_seqids -in Plar_chr.fa -out Masked_palm_civet -title Masked_palm_civet


# the script "../Scripts/bulk_bash.sh" used for execution of short-BLASTN with 100 % identity and e-value 0.8 thresholds, these parameters were identified through web-based blast searches
# this is done on the server 

#to run the script on the server
./bulk_bash.sh /path/to/input_directory/ /path/to/output_directory

#example for masked palm civet and scv2
./bulk_bash.sh ~/working_dir/20240301_SCV2_genome_chop/NC_045512.2_split.fa.split/ ~/working_dir/Local_blast/20240301_SCV2_vs_R.ferr/batch_blast_test/Plar_chr_out/
./bulk_bash.sh ~/working_dir/20240601_SCV1_genome_chop/AY274119.3_split.fa.split/ ~/working_dir/Local_blast/20240301_SCV2_vs_R.ferr/batch_blast_test/Plar_chr_out_SCV1/


# !once completed there will be too many files in the folder to cat so need to run a loop  
for i in *.blastn_out; do cat $i >> Plar_chr_out_all.csv; done # the csv is next imported into R in RStudio (see README)


##-------------------------------------------------------------------------
##                  processing of BLASTN results                          -
##-------------------------------------------------------------------------

# June 02, 2024 - added script from V. Montoya (blastRepair.txt), slightly adapted below:
gs <- read.csv(here('data/genomeSizes2.txt'),sep = '\t')

##SCV2

#bat
bat <- read.csv(here("data/bat_out_all.csv"),header = F)
colnames(bat) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
#library(tidyr)
#library(dplyr)
batc <- bat %>% group_by(qseqid) %>% summarise(n = n())
bat <- left_join(bat,batc, by = "qseqid")

#library(stringr)

bat[c('prefix','frag')] <- str_split_fixed(bat$qseqid, ':', 2)
bat[c('start','end')] <- str_split_fixed(bat$frag, '-', 2)
bat <- bat %>% mutate(center = as.integer(start) + 12)

bat <- bat %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

bat <- left_join(bat,gs,by ="Host")
bat <- bat %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

batu <- bat %>% distinct(center,.keep_all = T)
fields <- head(batu,1)

iterations = 29903
variables = 1
gdf <- matrix(ncol=variables, nrow=iterations)
for(i in 1:iterations){
  gdf[i,] <- i
  
}
gdf <- data.frame(gdf)
gdf <- rename(gdf,center = gdf)
gdf$n <- 0


batu <- full_join(batu,gdf,by = "center") 
batu <- batu %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

batu <- batu %>% mutate(Host = fields$Host)
batu <- batu %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
batu <- batu %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
batu <- batu %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
batout <- batu %>% select(center,countNorm)
batout[is.na(batout)] <- 0

#rename column
#dplyr

batout = batout %>% rename(Bat = countNorm)
write.table(batout, here('data/batout.csv'), row.names = F, col.name=T, sep = ";")

# removed from bat only positions 1764-1771 (outlier, magnitudes order greater Blast hits per Gbp)

batout_noout = read.table(here('data/batout_noout.txt'), header = T)

#sort

batout_sort = batout_noout %>%
  arrange((center))

# create a column in the dataframe that includes the running mean using the caTools package to calculate running averages per window size
# library(caTools)

batout_avg = runmean(batout_sort$Bat,20)
batout_sort['V3'] = batout_avg
colnames(batout_sort)[3] = "Bat_avg"
batout_sort$Bat_avg = as.numeric(as.character(batout_sort$Bat_avg))
batout_sort$Bat_avg = as.numeric(batout_sort$Bat_avg)

head(batout_sort, n=50)


#Malaysian porcupine
M_porcu <- read.csv(here("data/DSBC_Hbra_out_all.csv"),header = F)
colnames(M_porcu) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
M_porcuc <- M_porcu %>% group_by(qseqid) %>% summarise(n = n())
M_porcu <- left_join(M_porcu,M_porcuc, by = "qseqid")

M_porcu[c('prefix','frag')] <- str_split_fixed(M_porcu$qseqid, ':', 2)
M_porcu[c('start','end')] <- str_split_fixed(M_porcu$frag, '-', 2)
M_porcu <- M_porcu %>% mutate(center = as.integer(start) + 12)

M_porcu <- M_porcu %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

M_porcu <- left_join(M_porcu,gs,by ="Host")
M_porcu <- M_porcu %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

M_porcuu <- M_porcu %>% distinct(center,.keep_all = T)
fields <- head(M_porcuu,1)

M_porcuu <- full_join(M_porcuu,gdf,by = "center") 
M_porcuu <- M_porcuu %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

M_porcuu <- M_porcuu %>% mutate(Host = fields$Host)
M_porcuu <- M_porcuu %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
M_porcuu <- M_porcuu %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
M_porcuu <- M_porcuu %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
M_porcuout <- M_porcuu %>% select(center,countNorm)
M_porcuout[is.na(M_porcuout)] <- 0


#rename column
#dplyr

M_porcuout = M_porcuout %>% rename(Mporcu = countNorm)

#sort

M_porcuout_sort = M_porcuout %>%
  arrange((center))

#running average
M_porcuout_avg = runmean(M_porcuout_sort$Mporcu,20)
M_porcuout_sort['V3'] = M_porcuout_avg
colnames(M_porcuout_sort)[3] = "Mporcu_avg"
M_porcuout_sort$Mporcu_avg = as.numeric(as.character(M_porcuout_sort$Mporcu_avg))
head(M_porcuout_sort)


#Human

Human <- read.csv(here("data/Human_out_all.csv"),header = F)
colnames(Human) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
Humanc <- Human %>% group_by(qseqid) %>% summarise(n = n())
Human <- left_join(Human,Humanc, by = "qseqid")

Human[c('prefix','frag')] <- str_split_fixed(Human$qseqid, ':', 2)
Human[c('start','end')] <- str_split_fixed(Human$frag, '-', 2)
Human <- Human %>% mutate(center = as.integer(start) + 12)

Human <- Human %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

Human <- left_join(Human,gs,by ="Host")
Human <- Human %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

Human <- Human %>% distinct(center,.keep_all = T)
fields <- head(Human,1)

Human <- full_join(Human,gdf,by = "center") 
Human <- Human %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

Human <- Human %>% mutate(Host = fields$Host)
Human <- Human %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
Human <- Human %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
Human <- Human %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
Humanout <- Human %>% select(center,countNorm)
Humanout[is.na(Humanout)] <- 0

#rename column
#dplyr

Humanout = Humanout %>% rename(Human = countNorm)

#sort

Humanout_sort = Humanout %>%
  arrange((center))

#running average
Humanoutt_avg = runmean(Humanout_sort$Human,20)
Humanout_sort['V3'] = Humanoutt_avg
colnames(Humanout_sort)[3] = "Human_avg"
Humanout_sort$Human_avg = as.numeric(as.character(Humanout_sort$Human_avg))
head(Humanout_sort)


#Racoon Dog
RDog <- read.csv(here("data/RD_all_out.csv"),header = F)
colnames(RDog) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
RDogc <- RDog %>% group_by(qseqid) %>% summarise(n = n())
RDog <- left_join(RDog, RDogc, by = "qseqid")

RDog[c('prefix','frag')] <- str_split_fixed(RDog$qseqid, ':', 2)
RDog[c('start','end')] <- str_split_fixed(RDog$frag, '-', 2)
RDog <- RDog %>% mutate(center = as.integer(start) + 12)

RDog <- RDog %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

RDog <- left_join(RDog,gs,by ="Host")
RDog <- RDog %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

RDogu <- RDog %>% distinct(center,.keep_all = T)
fields <- head(RDogu,1)

RDogu <- full_join(RDogu,gdf,by = "center") 
RDogu <- RDogu %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

RDogu <- RDogu %>% mutate(Host = fields$Host)
RDogu <- RDogu %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
RDogu <- RDogu %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
RDogu <- RDogu %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
RDogout <- RDogu %>% select(center,countNorm)
RDogout[is.na(RDogout)] <- 0

#rename column
#dplyr

RDogout = RDogout %>% rename(RD = countNorm)

#sort

RDogout_sort = RDogout %>%
  arrange((center))

#running average
RDogout_avg = runmean(RDogout_sort$RD,20)
RDogout_sort['V3'] = RDogout_avg
colnames(RDogout_sort)[3] = "RDog_avg"
RDogout_sort$RDog_avg = as.numeric(as.character(RDogout_sort$RDog_avg))
head(RDogout_sort)



#Pangolin
Pang <- read.csv(here("data/Pang_out_all.csv"),header = F)
colnames(Pang) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
Pangc <- Pang %>% group_by(qseqid) %>% summarise(n = n())
Pang <- left_join(Pang,Pangc, by = "qseqid")

Pang[c('prefix','frag')] <- str_split_fixed(Pang$qseqid, ':', 2)
Pang[c('start','end')] <- str_split_fixed(Pang$frag, '-', 2)
Pang <- Pang %>% mutate(center = as.integer(start) + 12)

Pang <- Pang %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

Pang <- left_join(Pang,gs,by ="Host")
Pang <- Pang %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

Pangu <- Pang %>% distinct(center,.keep_all = T)
fields <- head(Pangu,1)

Pangu <- full_join(Pangu,gdf,by = "center") 
Pangu <- Pangu %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

Pangu <- Pangu %>% mutate(Host = fields$Host)
Pangu <- Pangu %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
Pangu <- Pangu %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
Pangu <- Pangu %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
Pangout <- Pangu %>% select(center,countNorm)
Pangout[is.na(Pangout)] <- 0

#rename column
#dplyr

Pangout = Pangout %>% rename(Pang = countNorm)

#sort

Pangout_sort = Pangout %>%
  arrange((center))

#running average
Pangout_avg = runmean(Pangout_sort$Pang,20)
Pangout_sort['V3'] = Pangout_avg
colnames(Pangout_sort)[3] = "Pang_avg"
Pangout_sort$Pang_avg = as.numeric(as.character(Pangout_sort$Pang_avg))
head(Pangout_sort)

write.table(Pangout, here('data/Pangout_new.csv'), row.names = F, col.name=T, sep = ";")



# July 6, 2025 - updated for Paguma larvata genome (masked palm civet), the initial pipeline was developed with the Asian palm civet genome (Paradoxurus hermaphroditus)
# Palm Civet
PCivet <- read.csv(here("data/Plar_chr_out_all_try.csv"),header = F)
colnames(PCivet) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
PCivetc <- PCivet %>% group_by(qseqid) %>% summarise(n = n())
PCivet <- left_join(PCivet,PCivetc, by = "qseqid")

PCivet[c('prefix','frag')] <- str_split_fixed(PCivet$qseqid, ':', 2)
PCivet[c('start','end')] <- str_split_fixed(PCivet$frag, '-', 2)
PCivet <- PCivet %>% mutate(center = as.integer(start) + 12)

PCivet <- PCivet %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

PCivet <- left_join(PCivet,gs,by ="Host")
PCivet <- PCivet %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

PCivetu <- PCivet %>% distinct(center,.keep_all = T)
fields <- head(PCivetu,1)

PCivetu <- full_join(PCivetu,gdf,by = "center") 
PCivetu <- PCivetu %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

PCivetu <- PCivetu %>% mutate(Host = fields$Host)
PCivetu <- PCivetu %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
PCivetu <- PCivetu %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
PCivetu <- PCivetu %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
PCivetout <- PCivetu %>% select(center,countNorm)
PCivetout[is.na(PCivetout)] <- 0

#rename column
#dplyr

PCivetout = PCivetout %>% rename(Pcivet = countNorm)

#sort

PCivetout_sort = PCivetout %>%
  arrange((center))

#running average
PCivetout_avg = runmean(PCivetout_sort$Pcivet,20)
PCivetout_sort['V3'] = PCivetout_avg
colnames(PCivetout_sort)[3] = "Pcivet_avg"
PCivetout_sort$Pcivet_avg = as.numeric(as.character(PCivetout_sort$Pcivet_avg))

#write.table(PCivetout, here('data/PCivetout_new.csv'), row.names = F, col.name=T, sep = ";")



#snake
snake <- read.csv(here("data/snake_out_all.csv"),header = F)
colnames(snake) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
snakec <- snake %>% group_by(qseqid) %>% summarise(n = n())
snake <- left_join(snake,snakec, by = "qseqid")

snake[c('prefix','frag')] <- str_split_fixed(snake$qseqid, ':', 2)
snake[c('start','end')] <- str_split_fixed(snake$frag, '-', 2)
snake <- snake %>% mutate(center = as.integer(start) + 12)

snake <- snake %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

snake <- left_join(snake,gs,by ="Host")
snake <- snake %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

snakeu <- snake %>% distinct(center,.keep_all = T)
fields <- head(snakeu,1)

snakeu <- full_join(snakeu,gdf,by = "center") 
snakeu <- snakeu %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

snakeu <- snakeu %>% mutate(Host = fields$Host)
snakeu <- snakeu %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
snakeu <- snakeu %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
snakeu <- snakeu %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
snakeout <- snakeu %>% select(center,countNorm)
snakeout[is.na(snakeout)] <- 0

#rename column
#dplyr

snakeout = snakeout %>% rename(Snake = countNorm)


#sort

snakeout_sort = snakeout %>%
  arrange((center))

#running average
snakeout_avg = runmean(snakeout_sort$Snake,20)
snakeout_sort['V3'] = snakeout_avg
colnames(snakeout_sort)[3] = "snake_avg"
snakeout_sort$snake_avg = as.numeric(as.character(snakeout_sort$snake_avg))



#chicken
chick <- read.csv(here("data/chick_out_all.csv"),header = F)
colnames(chick) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
chickc <- chick %>% group_by(qseqid) %>% summarise(n = n())
chick <- left_join(chick,chickc, by = "qseqid")

chick[c('prefix','frag')] <- str_split_fixed(chick$qseqid, ':', 2)
chick[c('start','end')] <- str_split_fixed(chick$frag, '-', 2)
chick <- chick %>% mutate(center = as.integer(start) + 12)

chick <- chick %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferremequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferremequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

chick <- left_join(chick,gs,by ="Host")
chick <- chick %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

chicku <- chick %>% distinct(center,.keep_all = T)
fields <- head(chicku,1)

chicku <- full_join(chicku,gdf,by = "center") 
chicku <- chicku %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

chicku <- chicku %>% mutate(Host = fields$Host)
chicku <- chicku %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
chicku <- chicku %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
chicku <- chicku %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
chickout <- chicku %>% select(center,countNorm)
chickout[is.na(chickout)] <- 0

#rename column
#dplyr

chickout = chickout %>% rename(Chick = countNorm)

#sort

chickout_sort = chickout %>%
  arrange((center))

#running average
chickout_avg = runmean(chickout_sort$Chick,20)
chickout_sort['V3'] = chickout_avg
colnames(chickout_sort)[3] = "chick_avg"
chickout_sort$chick_avg = as.numeric(as.character(chickout_sort$chick_avg))


# join SCV2 datasets together to replace the input

SCV2_result_ = merge(batout_sort,M_porcuout_sort, by = "center")
SCV2_result_ = merge(SCV2_result_,Humanout_sort, by = "center")
SCV2_result_ = merge(SCV2_result_,RDogout_sort, by = "center")
SCV2_result_ = merge(SCV2_result_,Pangout_sort, by = "center")
SCV2_result_ = merge(SCV2_result_,PCivetout_sort, by = "center")
SCV2_result_ = merge(SCV2_result_,snakeout_sort, by = "center")
SCV2_result_ = merge(SCV2_result_,chickout_sort, by = "center")
SCV2_result_ = SCV2_result_ %>% rename(Location = center)



#SCV1

bat_scv1 <- read.csv(here("data/bat_out_scv1_all.csv"),header = F)
colnames(bat_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
batc_scv1 <- bat_scv1 %>% group_by(qseqid) %>% summarise(n = n())
bat_scv1 <- left_join(bat_scv1,batc_scv1, by = "qseqid")

bat_scv1[c('prefix','frag')] <- str_split_fixed(bat_scv1$qseqid, ':', 2)
bat_scv1[c('start','end')] <- str_split_fixed(bat_scv1$frag, '-', 2)
bat_scv1 <- bat_scv1 %>% mutate(center = as.integer(start) + 12)

bat_scv1 <- bat_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

bat_scv1 <- left_join(bat_scv1,gs,by ="Host")
bat_scv1 <- bat_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

bat_scv1u <- bat_scv1 %>% distinct(center,.keep_all = T)
fields <- head(bat_scv1u,1)

iterations_scv1 = 29751
variables = 1
gdf_scv1 <- matrix(ncol=variables, nrow=iterations_scv1)
for(i in 1:iterations_scv1){
  gdf_scv1[i,] <- i
  
}
gdf_scv1 <- data.frame(gdf_scv1)
gdf_scv1 <- rename(gdf_scv1,center = gdf_scv1)
gdf_scv1$n <- 0


bat_scv1u <- full_join(bat_scv1u,gdf_scv1,by = "center") 
bat_scv1u <- bat_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

bat_scv1u <- bat_scv1u %>% mutate(Host = fields$Host)
bat_scv1u <- bat_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
bat_scv1u <- bat_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
bat_scv1u <- bat_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
bat_scv1out <- bat_scv1u %>% select(center,countNorm)
bat_scv1out[is.na(bat_scv1out)] <- 0

#rename column
#dplyr

bat_scv1out = bat_scv1out %>% rename(Bat = countNorm)

#sort

bat_scv1out_sort = bat_scv1out %>%
  arrange((center))

# running average

bat_scv1out_avg = runmean(bat_scv1out_sort$Bat,20)
bat_scv1out_sort['V3'] = bat_scv1out_avg
colnames(bat_scv1out_sort)[3] = "Bat_avg"
bat_scv1out_sort$Bat_avg = as.numeric(as.character(bat_scv1out_sort$Bat_avg))

head(bat_scv1out_sort, n=50)



#Malaysian porcupine
M_porcu_scv1 <- read.csv(here("data/DSBC_Hbra_out_scv1_all.csv"),header = F)
colnames(M_porcu_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
M_porcuc_scv1 <- M_porcu_scv1 %>% group_by(qseqid) %>% summarise(n = n())
M_porcu_scv1 <- left_join(M_porcu_scv1,M_porcuc_scv1, by = "qseqid")

M_porcu_scv1[c('prefix','frag')] <- str_split_fixed(M_porcu_scv1$qseqid, ':', 2)
M_porcu_scv1[c('start','end')] <- str_split_fixed(M_porcu_scv1$frag, '-', 2)
M_porcu_scv1 <- M_porcu_scv1 %>% mutate(center = as.integer(start) + 12)

M_porcu_scv1 <- M_porcu_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

M_porcu_scv1 <- left_join(M_porcu_scv1,gs,by ="Host")
M_porcu_scv1 <- M_porcu_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

M_porcu_scv1u <- M_porcu_scv1 %>% distinct(center,.keep_all = T)
fields <- head(M_porcu_scv1u,1)

M_porcu_scv1u <- full_join(M_porcu_scv1u,gdf_scv1,by = "center") 
M_porcu_scv1u <- M_porcu_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

M_porcu_scv1u <- M_porcu_scv1u %>% mutate(Host = fields$Host)
M_porcu_scv1u <- M_porcu_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
M_porcu_scv1u <- M_porcu_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
M_porcu_scv1u <- M_porcu_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
M_porcu_scv1out <- M_porcu_scv1u %>% select(center,countNorm)
M_porcu_scv1out[is.na(M_porcu_scv1out)] <- 0

#rename column
#dplyr

M_porcu_scv1out = M_porcu_scv1out %>% rename(Mporcu = countNorm)


#sort

M_porcu_scv1out_sort = M_porcu_scv1out %>%
  arrange((center))

#running average
M_porcu_scv1out_avg = runmean(M_porcu_scv1out_sort$Mporcu,20)
M_porcu_scv1out_sort['V3'] = M_porcu_scv1out_avg
colnames(M_porcu_scv1out_sort)[3] = "Mporcu_avg"
M_porcu_scv1out_sort$Mporcu_avg = as.numeric(as.character(M_porcu_scv1out_sort$Mporcu_avg))
head(M_porcu_scv1out_sort)



##Pangolin
Pang_scv1 <- read.csv(here("data/Pan_out_scv1_all.csv"),header = F)
colnames(Pang_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
Pangc_scv1 <- Pang_scv1 %>% group_by(qseqid) %>% summarise(n = n())
Pang_scv1 <- left_join(Pang_scv1,Pangc_scv1, by = "qseqid")

Pang_scv1[c('prefix','frag')] <- str_split_fixed(Pang_scv1$qseqid, ':', 2)
Pang_scv1[c('start','end')] <- str_split_fixed(Pang_scv1$frag, '-', 2)
Pang_scv1 <- Pang_scv1 %>% mutate(center = as.integer(start) + 12)

Pang_scv1 <- Pang_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

Pang_scv1 <- left_join(Pang_scv1,gs,by ="Host")
Pang_scv1 <- Pang_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

Pang_scv1u <- Pang_scv1 %>% distinct(center,.keep_all = T)
fields <- head(Pang_scv1u,1)

Pang_scv1u <- full_join(Pang_scv1u,gdf_scv1,by = "center") 
Pang_scv1u <- Pang_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

Pang_scv1u <- Pang_scv1u %>% mutate(Host = fields$Host)
Pang_scv1u <- Pang_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
Pang_scv1u <- Pang_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
Pang_scv1u <- Pang_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
Pang_scv1out <- Pang_scv1u %>% select(center,countNorm)
Pang_scv1out[is.na(Pang_scv1out)] <- 0

#rename column
#dplyr

Pang_scv1out = Pang_scv1out %>% rename(Pang = countNorm)

#sort

Pang_scv1out_sort = Pang_scv1out %>%
  arrange((center))

write.table(Pang_scv1out_sort, here('data/Pang_scv1_out_sort.csv'), row.names = F, col.name=T, sep = ";")

# removed from Pangolin only positions 22840-22847 by replacing with zero

pang_scv1out_noout = read.table(here('data/Pang_scv1_out_sort_nooutliers_.txt'), header = T)


#running average
Pangout_avg = runmean(pang_scv1out_noout$Pang,20)
pang_scv1out_noout['V3'] = Pangout_avg
colnames(pang_scv1out_noout)[3] = "Pang_avg"
pang_scv1out_noout$Pang_avg = as.numeric(as.character(pang_scv1out_noout$Pang_avg))
head(pang_scv1out_noout)



###Raccoon Dog

RDog_scv1 <- read.csv(here("data/NYPRO_out_scv1_all.csv"),header = F)
colnames(RDog_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
RDogc_scv1 <- RDog_scv1%>% group_by(qseqid) %>% summarise(n = n())
RDog_scv1 <- left_join(RDog_scv1, RDogc_scv1, by = "qseqid")

RDog_scv1[c('prefix','frag')] <- str_split_fixed(RDog_scv1$qseqid, ':', 2)
RDog_scv1[c('start','end')] <- str_split_fixed(RDog_scv1$frag, '-', 2)
RDog_scv1 <- RDog_scv1 %>% mutate(center = as.integer(start) + 12)

RDog_scv1 <- RDog_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

RDog_scv1 <- left_join(RDog_scv1,gs,by ="Host")
RDog_scv1 <- RDog_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

RDog_scv1u <- RDog_scv1 %>% distinct(center,.keep_all = T)
fields <- head(RDog_scv1u,1)

RDog_scv1u <- full_join(RDog_scv1u,gdf_scv1,by = "center") 
RDog_scv1u <- RDog_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

RDog_scv1u <- RDog_scv1u %>% mutate(Host = fields$Host)
RDog_scv1u <- RDog_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
RDog_scv1u <- RDog_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
RDog_scv1u <- RDog_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
RDog_scv1out <- RDog_scv1u %>% select(center,countNorm)
RDog_scv1out[is.na(RDog_scv1out)] <- 0

#rename column
#dplyr

RDog_scv1out = RDog_scv1out %>% rename(RD = countNorm)


#sort

RDog_scv1out_sort = RDog_scv1out %>%
  arrange((center))

#running average
RDog_scv1out_avg = runmean(RDog_scv1out_sort$RD,20)
RDog_scv1out_sort['V3'] = RDog_scv1out_avg
colnames(RDog_scv1out_sort)[3] = "RDog_avg"
RDog_scv1out_sort$RDog_avg = as.numeric(as.character(RDog_scv1out_sort$RDog_avg))
head(RDog_scv1out_sort)



###Human

Human_scv1 <- read.csv(here("data/Human_out_scv1_all.csv"),header = F)
colnames(Human_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
Huamnc_scv1 <- Human_scv1 %>% group_by(qseqid) %>% summarise(n = n())
Human_scv1 <- left_join(Human_scv1,Humanc_scv1, by = "qseqid")

Human_scv1[c('prefix','frag')] <- str_split_fixed(Human_scv1$qseqid, ':', 2)
Human_scv1[c('start','end')] <- str_split_fixed(Human_scv1$frag, '-', 2)
Human_scv1 <- Human_scv1 %>% mutate(center = as.integer(start) + 12)

Human_scv1 <- Human_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

Human_scv1 <- left_join(Human_scv1,gs,by ="Host")
Human_scv1 <- Human_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

Human_scv1u <- Human_scv1 %>% distinct(center,.keep_all = T)
fields <- head(Human_scv1u,1)

Human_scv1u <- full_join(Human_scv1u,gdf_scv1,by = "center") 
Human_scv1u <- Human_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

Human_scv1u <- Human_scv1u %>% mutate(Host = fields$Host)
Human_scv1u <- Human_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
Human_scv1u <- Human_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
Human_scv1u <- Human_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
Human_scv1out <- Human_scv1u %>% select(center,countNorm)
Human_scv1out[is.na(Human_scv1out)] <- 0

#rename column
#dplyr

Human_scv1out = Human_scv1out %>% rename(Human = countNorm)
#sort

Human_scv1out_sort = Human_scv1out %>%
  arrange((center))

write.table(Human_scv1out_sort, here('data/Humanout_scv1_sort.csv'), row.names = F, col.name=T, sep = ";")

# removed from human only positions 15261-15270 by replacing with zero

human_scv1out_noout = read.table(here('data/Humanout_scv1_sort_noout_.txt'), header = T)

#running average
Human_scv1out_avg = runmean(human_scv1out_noout$Human,20)
human_scv1out_noout['V3'] = Human_scv1out_avg
colnames(human_scv1out_noout)[3] = "Human_avg"
human_scv1out_noout$Human_avg = as.numeric(as.character(human_scv1out_noout$Human_avg))
head(human_scv1out_noout)



### Palm Civet

PCivet_scv1 <- read.csv(here("data/Plar_chr_out_all_SCV1_try.csv"),header = F)
colnames(PCivet_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
PCivetc_scv1 <- PCivet_scv1 %>% group_by(qseqid) %>% summarise(n = n())
PCivet_scv1 <- left_join(PCivet_scv1,PCivetc_scv1, by = "qseqid")

PCivet_scv1[c('prefix','frag')] <- str_split_fixed(PCivet_scv1$qseqid, ':', 2)
PCivet_scv1[c('start','end')] <- str_split_fixed(PCivet_scv1$frag, '-', 2)
PCivet_scv1 <- PCivet_scv1 %>% mutate(center = as.integer(start) + 12)

PCivet_scv1 <- PCivet_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

PCivet_scv1 <- left_join(PCivet_scv1,gs,by ="Host")
PCivet_scv1 <- PCivet_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

PCivet_scv1u <- PCivet_scv1 %>% distinct(center,.keep_all = T)
fields <- head(PCivet_scv1u,1)

PCivet_scv1u <- full_join(PCivet_scv1u,gdf_scv1,by = "center") 
PCivet_scv1u <- PCivet_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

PCivet_scv1u <- PCivet_scv1u %>% mutate(Host = fields$Host)
PCivet_scv1u <- PCivet_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
PCivet_scv1u <- PCivet_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
PCivet_scv1u <- PCivet_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
PCivet_scv1out <- PCivet_scv1u %>% select(center,countNorm)
PCivet_scv1out[is.na(PCivet_scv1out)] <- 0

#rename column
#dplyr

PCivet_scv1out = PCivet_scv1out %>% rename(Pcivet = countNorm)


#sort

PCivet_scv1out_sort = PCivet_scv1out %>%
  arrange((center))

#running average
PCivet_scv1out_avg = runmean(PCivet_scv1out_sort$Pcivet,20)
PCivet_scv1out_sort['V3'] = PCivet_scv1out_avg
colnames(PCivet_scv1out_sort)[3] = "Pcivet_avg"
PCivet_scv1out_sort$Pcivet_avg = as.numeric(as.character(PCivet_scv1out_sort$Pcivet_avg))



##snake

snake_scv1 <- read.csv(here("data/snake_out_scv1_all.csv"),header = F)
colnames(snake_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
snakec_scv1 <- snake_scv1 %>% group_by(qseqid) %>% summarise(n = n())
snake_scv1 <- left_join(snake_scv1,snakec_scv1, by = "qseqid")

snake_scv1[c('prefix','frag')] <- str_split_fixed(snake_scv1$qseqid, ':', 2)
snake_scv1[c('start','end')] <- str_split_fixed(snake_scv1$frag, '-', 2)
snake_scv1 <- snake_scv1 %>% mutate(center = as.integer(start) + 12)

snake_scv1 <- snake_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

snake_scv1 <- left_join(snake_scv1,gs,by ="Host")
snake_scv1 <- snake_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

snake_scv1u <- snake_scv1 %>% distinct(center,.keep_all = T)
fields <- head(snake_scv1u,1)

snake_scv1u <- full_join(snake_scv1u,gdf_scv1,by = "center") 
snake_scv1u <- snake_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

snake_scv1u <- snake_scv1u %>% mutate(Host = fields$Host)
snake_scv1u <- snake_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
snake_scv1u <- snake_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
snake_scv1u <- snake_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
snake_scv1out <- snake_scv1u %>% select(center,countNorm)
snake_scv1out[is.na(snake_scv1out)] <- 0

#rename column
#dplyr

snake_scv1out = snake_scv1out %>% rename(Snake = countNorm)


#sort

snake_scv1out_sort = snake_scv1out %>%
  arrange((center))

write.table(snake_scv1out_sort, here('data/Snake_scv1_out_sort.csv'), row.names = F, col.name=T, sep = ";")

# removed from Pangolin only positions 3189-3198 by replacing with zero

snake_scv1out_noout = read.table(here('data/Snake_scv1_out_sort_nooutliers_.txt'), header = T)

#running average
snake_scv1out_avg = runmean(snake_scv1out_noout$Snake,20)
snake_scv1out_noout['V3'] = snake_scv1out_avg
colnames(snake_scv1out_noout)[3] = "snake_avg"
snake_scv1out_noout$snake_avg = as.numeric(as.character(snake_scv1out_noout$snake_avg))



##chick

chick_scv1 <- read.csv(here("data/chick_out_scv1_all.csv"),header = F)
colnames(chick_scv1) <- c('qseqid','sseqid','length','qlen','slen','qstart','qend','sstart','send','pident','mismatch','gapopen','evalue','bitscore')
chickc_scv1 <- chick_scv1 %>% group_by(qseqid) %>% summarise(n = n())
chick_scv1 <- left_join(chick_scv1,chickc_scv1, by = "qseqid")


chick_scv1[c('prefix','frag')] <- str_split_fixed(chick_scv1$qseqid, ':', 2)
chick_scv1[c('start','end')] <- str_split_fixed(chick_scv1$frag, '-', 2)
chick_scv1 <- chick_scv1 %>% mutate(center = as.integer(start) + 12)

chick_scv1 <- chick_scv1 %>% mutate(Host = case_when(
  grepl('Rhinolophus_ferrumequinum', sseqid, fixed = TRUE) ~ "Rhinolophus_ferrumequinum",
  grepl('Nyctereutes_procyonoides', sseqid, fixed = TRUE) ~ "Nyctereutes_procyonoides",
  grepl('Manis_pentadactyla', sseqid, fixed = TRUE) ~ "Manis_pentadactyla",
  grepl('Homo_sapien', sseqid, fixed = TRUE) ~ "Homo_sapien",
  grepl('Paguma_larvata', sseqid, fixed = TRUE) ~ "Paguma_larvata",
  grepl('Hystrix_brachyura', sseqid, fixed = TRUE) ~ "Hystrix_brachyura",
  grepl('Gallus_gallus', sseqid, fixed = TRUE) ~ "Gallus_gallus",
  grepl('Diadophis_punctatus_similis', sseqid, fixed = TRUE) ~ "Diadophis_punctatus_similis",
  
))

chick_scv1 <- left_join(chick_scv1,gs,by ="Host")
chick_scv1 <- chick_scv1 %>% mutate(countNorm = round((n/GenomeSize)*1e9,3))

chick_scv1u <- chick_scv1 %>% distinct(center,.keep_all = T)
fields <- head(chick_scv1u,1)

chick_scv1u <- full_join(chick_scv1u,gdf_scv1,by = "center") 
chick_scv1u <- chick_scv1u %>% mutate(countNorm = ifelse(is.na(countNorm),0,countNorm))

chick_scv1u <- chick_scv1u %>% mutate(Host = fields$Host)
chick_scv1u <- chick_scv1u %>% mutate(GenomeSize = ifelse(is.na(GenomeSize),fields$GenomeSize,GenomeSize))
chick_scv1u <- chick_scv1u %>% mutate(start = ifelse(is.na(start),ifelse((center - 12) < 0,0,center - 12),start))
chick_scv1u <- chick_scv1u %>% mutate(end = ifelse(is.na(end),center + 12,end))

##select the columns
chick_scv1out <- chick_scv1u %>% select(center,countNorm)
chick_scv1out[is.na(chick_scv1out)] <- 0

#rename column
#dplyr

chick_scv1out = chick_scv1out %>% rename(Chick = countNorm)

#sort

chick_scv1out_sort = chick_scv1out %>%
  arrange((center))

#running average
chick_scv1out_avg = runmean(chick_scv1out_sort$Chick,20)
chick_scv1out_sort['V3'] = chick_scv1out_avg
colnames(chick_scv1out_sort)[3] = "chick_avg"
chick_scv1out_sort$chick_avg = as.numeric(as.character(chick_scv1out_sort$chick_avg))

# join SCV1 datasets together to replace the input

SCV1_result_ = merge(bat_scv1out_sort,M_porcu_scv1out_sort, by = "center")
SCV1_result_ = merge(SCV1_result_,human_scv1out_noout, by = "center")
SCV1_result_ = merge(SCV1_result_,RDog_scv1out_sort, by = "center")
SCV1_result_ = merge(SCV1_result_,pang_scv1out_noout, by = "center")
SCV1_result_ = merge(SCV1_result_,PCivet_scv1out_sort, by = "center")
SCV1_result_ = merge(SCV1_result_,snake_scv1out_noout, by = "center")
SCV1_result_ = merge(SCV1_result_,chick_scv1out_sort, by = "center")
SCV1_result_ = SCV1_result_ %>% rename(Location = center)

#29751


##-------------------------------------------------------------------------
##                          PCA for outliers                              -
##-------------------------------------------------------------------------

## PCA to explore the data for outliers - this was not presented in the manuscript - outliers removed are provided above

# below is older script, see above - use SCV1_result_ as object - June 9, 2024, and the final heatmap plots of normalized rolling averages are provided at the end of this script


#library(ggplot2)
#library(viridis)


heat = read.table(here('data/To_r_4_heatmap.txt'), header=T)
head(heat)

library(tidyr)
heat2 = heat %>% gather(sample, value, -Location)
#hard to see try log
heat3 = heat2
heat3$value = heat2$value+1
head(heat3)


heat3[,3:3] = log2(heat3[3:3])
head(heat3)
tail(heat3)
ggplot(heat3, aes(x = sample, y = Location, fill = value)) + geom_tile() + scale_fill_viridis(option="plasma",limits = c(0, 350), direction = -1) + ylab("SCV2 genome") + xlab("host")
ggplot(heat2, aes(x = sample, y = Location, fill = value)) + geom_tile() + scale_fill_viridis(option="turbo", direction=-1) + ylab("SCV2 genome") + xlab("host")


ggplot(heat2, aes(x = Location, y = sample, fill = value)) + geom_tile() + scale_fill_gradient(low="gray100", high ="midnightblue", limit = c(0, 165))+ ylab("host") + xlab("SCV2 genome")+
  theme_bw()

##the bat is a lot different than the other hosts, use pca to identify outliers
##https://www.geeksforgeeks.org/how-to-make-pca-plot-with-r/
heat.pca = prcomp(heat2[,c(1,3)],
                  center = TRUE,
                  scale. = TRUE)

summary(heat.pca)
str(heat.pca)

library(ggfortify)

heat.pca.plot = autoplot(heat.pca,
                         data=heat2, color = 'sample')


# outliers in bat at locations: 1764 - 1771 ranging from 894-897 total blast hits
# replaced these with zeros

heat_nout = read.table(here('data/To_r_4_heatmap_no_outlier.txt'), header=T)
heat2_nout = heat_nout %>% gather(sample, value, -Location)


heat_ = read.table(here('data/20240524_SCV2_blast_chick_snake_noout.txt'),header=T) # this is normalized 
heat_2 = heat_ %>% gather(sample, value, -Location)


# June 2 - SCV2_result object to replace heat_ 
heat_3 = SCV2_result %>% gather(sample, value, -Location) #works
write.table(heat_3,file = here('data/20240603_norm_counts_SCV2_heatmap.txt'), row.names = F, col.names=T)

heat_3_scv1 = SCV1_result %>% gather(sample, value, -Location) #works
write.table(heat_3_scv1,file = here('data/20240603_norm_counts_SCV1_heatmap.txt'), row.names = F, col.names=T)

# additional PCA exploratory plots below

heatnout.pca = prcomp(heat2_nout[,c(1,3)],
                      center = TRUE,
                      scale. = TRUE)
summary(heatnout.pca)
str(heatnout.pca)

heatnout.pca.plot = autoplot(heatnout.pca,
                             label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                             data=heat2_nout, color = 'sample')+ theme_bw()


heatnout2.pca = prcomp(heat_2[,c(1,3)],
                      center = TRUE,
                      scale. = TRUE)
summary(heatnout2.pca)
str(heatnout2.pca)


heatnout2.pca.plot = autoplot(heatnout2.pca,
                             label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                             data=heat_2, color = 'sample')+ theme_bw()


heatnout3.pca = prcomp(heat_3[,c(1,3)],
                       center = TRUE,
                       scale. = TRUE)
summary(heatnout3.pca)
str(heatnout3.pca)


heatnout3.pca.plot = autoplot(heatnout3.pca,
                              label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                              data=heat_2, color = 'sample')+ theme_bw()


heatnout3_scv1.pca = prcomp(heat_3_scv1[,c(1,3)],
                       center = TRUE,
                       scale. = TRUE)
summary(heatnout3_scv1.pca)
str(heatnout3_scv1.pca)


heatnout3_scv1.pca.plot = autoplot(heatnout3_scv1.pca,
                              label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                              data=heat_3_scv1, color = 'sample')+ theme_bw()


##-------------------------------------------------------------------------
##                normalization by genome length                          -
##-------------------------------------------------------------------------

# normalizing by "hits per billion" - with the billion being # of bps in the reference genome for each host

heat_nout_norm = read.table(here('data/To_r_4_heatmap_no_outlier_norm.txt'), header=T)
heat2_nout_norm = heat_nout_norm %>% gather(sample, value, -Location)

heatnout_norm.pca = prcomp(heat2_nout_norm[,c(1,3)],
                           center = TRUE,
                           scale. = TRUE)

heatnout_norm.pca.plot = autoplot(heatnout_norm.pca,
                                  label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                                  data=heat2_nout_norm, color = 'sample')+ theme_bw()

# this can help pick out clusters that are most specific to a given host compared to the others

# test with the 'zeros' removed

heat_nout_norm_nozero = read.table(here('data/hits_norm_nozero_4PCA.txt'), header=T)
heat2_nout_norm_nozero = heat_nout_norm_nozero %>% gather(sample, value, -Location)

heatnout_norm_nozero.pca = prcomp(heat2_nout_norm_nozero[,c(1,3)],
                                  center = TRUE,
                                  scale. = TRUE)

heatnout_norm_nozero.pca.plot = autoplot(heatnout_norm_nozero.pca,
                                         label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                                         data=heat2_nout_norm_nozero, color = 'sample')+ theme_bw()



# re-plot the normalize without the outlier for the heatmap
# this helps to avoid some of the overlapping points


out_plot = ggplot(heat2_nout_norm, aes(x = Location, y = sample, fill = value)) + geom_tile() + scale_fill_gradient2(low="gray100", mid="mediumblue", high ="midnightblue", midpoint=30)+ ylab("host") + xlab("SCV2 genome")+
  theme_bw()

out_plot2 = out_plot + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

pdf(here("20240308_Heatmap_blast_SCV2.pdf"), width=8.5, height=1) # whole-genome not presented in the manuscript
out_plot2
dev.off()

heat_2$sample = factor(heat_2$sample, levels = c("Snake", "Chick", "Mporcu", "Pcivet", "RD", "Pang","Bat", "Pang"))


out_plot3 = ggplot(heat_2, aes(x = Location, y = sample, fill = value)) + geom_tile() + scale_fill_gradient2(low="gray100", mid="mediumblue", high ="midnightblue", midpoint=30)+ ylab("host") + xlab("SCV2 genome")+
  theme_bw()

out_plot4 = out_plot3 + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

pdf(here("20240524_Heatmap_blast_SCV2_chicksnake_noout.pdf"), width=8.5, height=2)
out_plot4
dev.off()

#June 3 - removed from bat only positions 1764-1771 (outlier, magnitudes order greater Blast hits per Gbp)

heat3_nout_norm = read.table(here("data/20240603_norm_counts_SCV2_heatmap_nooutlier.txt"), header=T)
# heat3_nout_norm = heat3_nout_norm %>% gather(sample, value, -Location)

heat3nout_norm.pca = prcomp(heat3_nout_norm[,c(1,3)],
                           center = TRUE,
                           scale. = TRUE)

heat3nout_norm.pca.plot = autoplot(heat3nout_norm.pca,
                                  label = TRUE, label.label = "Location", label.size=1, label.vjust=1, label.hjust =2,
                                  data=heat3_nout_norm, color = 'sample')+ theme_bw()

pdf(here("data/20240603_Heatmap_blast_SCV2_chicksnake_noout.pdf"), width=8.5, height=2)
heat3nout_norm.pca.plot
dev.off()

heat3_nout_norm_scv1$sample = factor(heat3_nout_norm_scv1$sample, levels = c("Snake", "Chick", "Mporcu", "Pcivet", "RD", "Pang","Bat", "Human"))


out_plot7 = ggplot(heat3_nout_norm, aes(x = Location, y = sample, fill = value)) + geom_tile() + scale_fill_gradient2(low="gray100", mid="mediumblue", high ="midnightblue", midpoint=30)+ ylab("host") + xlab("SCV2 genome")+
  theme_bw()

out_plot8 = out_plot7 + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))


# June 9, 2024 - re-try, but log2 + 1 the values
heat3_nout_norm_log = heat3_nout_norm
heat3_nout_norm_log$transformed_value = log(heat3_nout_norm_log$value + 1) 


out_plot7_log = ggplot(heat3_nout_norm_log, aes(x = Location, y = sample, fill = transformed_value)) + geom_tile() + scale_fill_gradient2(low="gray100", mid="mediumblue", high ="midnightblue", midpoint=3)+ ylab("host") + xlab("SCV2 genome")+
  theme_bw()

out_plot8_log = out_plot7_log + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

# the pattern is not visible with log2 or log10

##-------------------------------------------------------------------------
##         update main manuscript figures with rolling average            -
##-------------------------------------------------------------------------

# June 9, 2024 - re-gather the results including the running averages

SCV2_result_avg = SCV2_result_[c("Location", "Bat_avg", "Mporcu_avg","Human_avg","RDog_avg", "Pang_avg", "Pcivet_avg", "snake_avg", "chick_avg")]
heat_4_avg = SCV2_result_avg %>% gather(sample, value, -Location) # works

write.table(heat_4_avg,here("data/20250705_heat_4_avg.csv"), row.names=F, col.names=T, sep = ";")

heat_4_avg$sample = factor(heat_4_avg$sample, levels = c("snake_avg", "chick_avg", "Mporcu_avg", "Pcivet_avg", "RDog_avg", "Pang_avg","Bat_avg", "Human_avg"))

#library(ggplot2)
out_plot_scv2avg = ggplot(heat_4_avg, aes(x = Location, y = sample, fill = value)) + geom_tile() + scale_fill_gradient2(low="white", mid="mediumblue", high ="midnightblue", midpoint=20)+ ylab("host") + xlab("SCV2 genome")+
  theme_bw()

out_plot_scv2avg_ = out_plot_scv2avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

out_plot_scv2avg_b = out_plot_scv2avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(28274, 30000)+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

# N ORF + 3 UTR
pdf(here("data/20250706_Heatmap_blast_SCV2_chicksnake_noout_N3utr_avg_win20.pdf"), width=7.5, height=2)
out_plot_scv2avg_b
dev.off()

# S ORF
out_scv2_plotS_avg = out_plot_scv2avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(21563, 25384)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20250707_Heatmap_blast_SCV2_chicksnake_noout_S_ang.pdf"), width=7.5, height=2)
out_scv2_plotS_avg
dev.off()

# July 7, 2025

# nsp2
out_scv2_plot_nsp2_avg = out_plot_scv2avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(806, 2719)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240720_Heatmap_blast_SCV2_chicksnake_noout_nsp2_avg.pdf"), width=7.5, height=2)
out_scv2_plot_nsp2_avg
dev.off()

# nsp3
out_scv2_plot_nsp3_avg = out_plot_scv2avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(2720, 8554)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240720_Heatmap_blast_SCV2_chicksnake_noout_nsp3_avg.pdf"), width=7.5, height=2)
out_scv2_plot_nsp3_avg
dev.off()

# E
out_scv2_plot_E_avg = out_plot_scv2avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(26245, 26472)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240720_Heatmap_blast_SCV2_chicksnake_noout_E_avg.pdf"), width=4.25, height=2)
out_scv2_plot_E_avg
dev.off()


# scv1_avg
SCV1_result_avg = SCV1_result_[c("Location", "Bat_avg", "Mporcu_avg","Human_avg","RDog_avg", "Pang_avg", "Pcivet_avg", "snake_avg", "chick_avg")]
heat_5_avg = SCV1_result_avg %>% gather(sample, value, -Location) #works

write.table(heat_5_avg,here("data/20250706_heat_5_avg.csv"), row.names=F, col.names=T, sep = ";")

heat_5_avg$sample = factor(heat_5_avg$sample, levels = c("snake_avg", "chick_avg", "Mporcu_avg", "Pcivet_avg", "RDog_avg", "Pang_avg","Bat_avg", "Human_avg"))

out_plot_scv1avg = ggplot(heat_5_avg, aes(x = Location, y = sample, fill = value)) + geom_tile() + scale_fill_gradient2(low="white", mid="mediumblue", high ="midnightblue", midpoint=45)+ ylab("host") + xlab("SCV1 genome")+
  theme_bw()

out_plot_scv1avg_ = out_plot_scv1avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

out_plot_scv1avg_b = out_plot_scv1avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(28120, 29750)+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4),legend.key.size = unit(.2,'cm'))

pdf(here("data/20250706_Heatmap_blast_SCV1_chicksnake_noout_N3utr_avg_win20.pdf"), width= 7.5, height=2)
out_plot_scv1avg_b
dev.off()


# S ORF
out_plotS_scv1_avg = out_plot_scv1avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(21492, 25259)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240719_Heatmap_blast_SCV1_chicksnake_noout_S_avg.pdf"), width=7.5, height=2)
out_plotS_scv1_avg
dev.off()

# nsp2
out_scv1_plot_nsp2_avg = out_plot_scv1avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(805, 2718)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240720_Heatmap_blast_SCV1_chicksnake_noout_nsp2_avg.pdf"), width=7.5, height=2)
out_scv1_plot_nsp2_avg
dev.off()

# nsp3
out_scv1_plot_nsp3_avg = out_plot_scv1avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(2719, 8484)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240720_Heatmap_blast_SCV1_chicksnake_noout_nsp3_avg.pdf"), width=7.5, height=2)
out_scv1_plot_nsp3_avg
dev.off()

# E ORF
out_scv1_plot_E_avg = out_plot_scv1avg + guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  xlim(26117, 26347)+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.size = unit(.5,'cm'))

pdf(here("data/20240720_Heatmap_blast_SCV1_chicksnake_noout_E_avg.pdf"), width=4.25, height=2)
out_scv1_plot_E_avg
dev.off()

#### End of script

