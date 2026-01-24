## ------------------------------------------------------------------------
##
## Script name: 20251101_DEG_Venn_Enrichment_PaulsonAR_2025.R
##
## Purpose of script: Functional enrichment of known SARS-CoV-2 responsive Differential Expressed Genes (DEG) in Calu-3 cells, human lung, and patients with severe COVID-19 
## 
##
##
## Author: Amber R. Paulson, PhD.
##
## Date Created: 2025-11-01
##
## Copyright (c) Amber Paulson, 2025 - governed by CC-BY-NC 4.0 - Attribution-Non-Commercial 4.0 International 
## see https://creativecommons.org/licenses/by-nc/4.0/
##
## Email: Amber[dot]Rose[dot]Paulson[at]gmail[dot]com
## Github: https://github.com/damselflywingz/EpiVIRTEX/
##
## ------------------------------------------------------------------------
##
## Notes: Analysis for enrichment of  gene targets of svRNA-TFOs and DEGs 
## from three transcriptomes studies based on available unique Gene symbols
## (see Materials and Methods in manuscript)
##        
##        
## ------------------------------------------------------------------------


##-------------------------------------------------------------------------
##                                Venn diagrams                           -
##-------------------------------------------------------------------------

# adapted from https://www.r-graph-gallery.com/14-venn-diagramm.html

# install.packages("VennDiagram")
library(VennDiagram)
library(here)
# need to get a vectored list 

E_in= read.table(here("data/E_targets.txt"), header =T)
N_in = read.table(here("data/N_targets.txt"), header =T)
nsp2_in = read.table(here("data/nsp2_targets.txt"), header =T)
nsp3_in = read.table(here("data/nsp3_targets.txt"), header =T)
S.1_in = read.table(here("data/S_targets.txt"), header =T)
S.2_in = read.table(here("data/S.2_targets.txt"), header=T)
All_in_ = read.table(here("data/All_targets_six_svRNA_TFOs.txt"), header = T)

# Aug 2, 2024 trying more venns / ORF region

# E
# need to get a vectored list 

calu3_DE = E_in$Calu3_DE
remove = "-"

calu3_DE = setdiff(calu3_DE, remove)

lung_DE = E_in$Lung_DE
lun_DE = setdiff(lung_DE, remove)

COVID19_up = E_in$COVID19_up
COVID19_up = setdiff(COVID19_up, remove)


## Chart
venn.diagram(
  x = list(calu3_DE, lung_DE, COVID19_up),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = here("data/20240803_venn_genes_E_only.png"),
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-E"
  
)


E.venn = venn.diagram(
  x = list(calu3_DE, lung_DE, COVID19_up),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = NULL,
  output=TRUE,
                      
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
                      
  # Circles
  lwd = 2,
  lty = 'blank',
  cex = 0.75,
  fill = myCol3, 
                      
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-E"
                  
)

grid.newpage()
grid.draw(E.venn)


lapply(E.venn, names)
lapply(E.venn,function(i) i$label)

x = list(calu3_DE, lung_DE, COVID19_up)
overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

E.venn[[7]]$label  <- paste(overlaps$a1, collapse="\n")  
# 
E.venn[[8]]$label  <- paste(overlaps$a2, collapse="\n")  
# 
E.venn[[9]]$label  <- paste(overlaps$a3, collapse="\n")  
# 
E.venn[[10]]$label  <- paste(overlaps$a4, collapse="\n")  
# 
E.venn[[11]]$label  <- paste(overlaps$a6, collapse="\n")  
# 
E.venn[[12]]$label  <- paste(overlaps$a7, collapse="\n")  
# 
 
grid.newpage()
grid.draw(E.venn) # save this as svg image


## N
## need to get a vectored list 
calu3_DE_N = N_in$Calu3_DE
remove = "-"

calu3_DE_N = setdiff(calu3_DE_N, remove)

lung_DE_N = N_in$Lung_DE
lung_DE_N = setdiff(lung_DE_N, remove)

COVID19_up_N = N_in$COVID19_up
COVID19_up_N = setdiff(COVID19_up_N, remove)

## Chart
N.venn = venn.diagram(
  x = list(calu3_DE_N, lung_DE_N, COVID19_up_N),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = NULL,
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  cex = 0.55,
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-N"
  
)

grid.newpage()
grid.draw(N.venn)


lapply(N.venn, names)
lapply(N.venn,function(i) i$label)

x = list(calu3_DE_N, lung_DE_N, COVID19_up_N)
overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

N.venn[[7]]$label  <- paste(overlaps$a1, collapse="\n")  
# 
N.venn[[8]]$label  <- paste(overlaps$a2, collapse="\n")  
# 
N.venn[[9]]$label  <- paste(overlaps$a3, collapse="\n")  
# 
N.venn[[10]]$label  <- paste(overlaps$a4, collapse="\n")  
# 
N.venn[[11]]$label  <- paste(overlaps$a6, collapse="\n")  
# 
N.venn[[12]]$label  <- paste(overlaps$a7, collapse="\n")  
# 

grid.newpage()
grid.draw(N.venn) # save this as svg image

## nsp2
## need to get a vectored list 

calu3_DE_nsp2 = nsp2_in$Calu3_DE
calu3_DE_nsp2 = setdiff(calu3_DE_nsp2, remove)

lung_DE_nsp2 = nsp2_in$Lung_DE
lung_DE_nsp2 = setdiff(lung_DE_nsp2, remove)

COVID19_up_nsp2 = nsp2_in$COVID19_up
COVID19_up_nsp2 = setdiff(COVID19_up_nsp2, remove)

## Chart
venn.diagram(
  x = list(calu3_DE_nsp2, lung_DE_nsp2, COVID19_up_nsp2),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = here("data/20240802_venn_genes_nsp2_only.png"),
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-nsp2"
  
)

# try to get labels for the manuscript - https://stackoverflow.com/questions/25019794/venn-diagram-with-item-labels

x = list(calu3_DE_nsp2, lung_DE_nsp2, COVID19_up_nsp2)

nsp2.venn = venn.diagram(x,
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = NULL,
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  cex = 0.5,
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-nsp2"
  
)

grid.newpage()
grid.draw(nsp2.venn)


lapply(nsp2.venn, names)
lapply(nsp2.venn,function(i) i$label)

overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

# 
#nsp2.venn[[7]]$label  <- paste(overlaps$a1) 
# 
nsp2.venn[[7]]$label  <- paste(overlaps$a1, collapse="\n")  
# 
nsp2.venn[[8]]$label  <- paste(overlaps$a2, collapse="\n")  
# 
nsp2.venn[[9]]$label  <- paste(overlaps$a3, collapse="\n")  
# 
nsp2.venn[[10]]$label <- paste(overlaps$a4, collapse="\n") 
#
nsp2.venn[[11]]$label <- paste(overlaps$a6, collapse="\n")  
# 
nsp2.venn[[12]]$label <- paste(overlaps$a7, collapse="\n") 

grid.newpage()
grid.draw(nsp2.venn)

## nsp3
## need to get a vectored list 

calu3_DE_nsp3 = nsp3_in$Calu3_DE
calu3_DE_nsp3 = setdiff(calu3_DE_nsp3, remove)

lung_DE_nsp3 = nsp3_in$Lung_DE
lung_DE_nsp3 = setdiff(lung_DE_nsp3, remove)

COVID19_up_nsp3 = nsp3_in$COVID19_up
COVID19_up_nsp3 = setdiff(COVID19_up_nsp3, remove)

## Chart
venn.diagram(
  x = list(calu3_DE_nsp3, lung_DE_nsp3, COVID19_up_nsp3),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = here("data/20240803_venn_genes_nsp3_only.png"),
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-nsp2"
  
)

x = list(calu3_DE_nsp3, lung_DE_nsp3, COVID19_up_nsp3)

nsp3.venn = venn.diagram(x,
                         category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
                         filename = NULL,
                         output=TRUE,
                         
                         #Output features
                         #height = 480,
                         #width = 480 , 
                         #resolution = 300,
                         #compression = "lzw",
                         
                         # Circles
                         lwd = 2,
                         lty = 'blank',
                         cex = 0.75,
                         fill = myCol3, 
                         
                         # Numbers
                         fontface = "bold",
                         fontfamily = "sans",
                         main = "svRNA-TFO-nsp3"
                         
)

grid.newpage()
grid.draw(nsp3.venn) # no intersections to report


lapply(nsp3.venn, names)
lapply(nsp3.venn,function(i) i$label)

overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

nsp3.venn[[7]]$label  <- paste(overlaps$a1, collapse="\n")  
 
grid.newpage()
grid.draw(nsp3.venn)

## S.1
## need to get a vectored list 
calu3_DE_S = S_in$Calu3_DE
remove = "-"

calu3_DE_S = setdiff(calu3_DE_S, remove)

lung_DE_S = S_in$Lung_DE
lung_DE_S = setdiff(lung_DE_S, remove)

COVID19_up_S = S_in$COVID19_up
COVID19_up_S = setdiff(COVID19_up_S, remove)

## Chart
S.venn = venn.diagram(
  x = list(calu3_DE_S, lung_DE_S, COVID19_up_S),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = NULL,
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  cex = 0.75,
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-S"
  
)

grid.newpage()
grid.draw(S.venn)

lapply(S.venn, names)
lapply(S.venn,function(i) i$label)

x = list(calu3_DE_S, lung_DE_S, COVID19_up_S)
overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

S.venn[[7]]$label  <- paste(overlaps$a3, collapse="\n")  
# 
S.venn[[8]]$label  <- paste(overlaps$a6, collapse="\n")  
# 
S.venn[[9]]$label  <- paste(overlaps$a7, collapse="\n")  
# 
S.venn[[10]]$label  <- paste(overlaps$a4, collapse="\n")  
# 
S.venn[[11]]$label  <- paste(overlaps$a1, collapse="\n")  
# 

grid.newpage()
grid.draw(S.venn) # save this as svg image

## S.2
## need to get a vectored list 
calu3_DE_S.2 = S.2_in$Calu3_DE
remove = "-"

calu3_DE_S.2 = setdiff(calu3_DE_S.2, remove)

lung_DE_S.2 = S.2_in$Lung_DE
lung_DE_S.2 = setdiff(lung_DE_S.2, remove)

COVID19_up_S.2 = S.2_in$COVID19_up
COVID19_up_S.2 = setdiff(COVID19_up_S.2, remove)

## Chart
S.2.venn = venn.diagram(
  x = list(calu3_DE_S.2, lung_DE_S.2, COVID19_up_S.2),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = NULL,
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  cex = 0.75,
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "svRNA-TFO-S.2"
  
)

grid.newpage()
grid.draw(S.2.venn)

lapply(S.2.venn, names)
lapply(S.2.venn,function(i) i$label)

x = list(calu3_DE_S.2, lung_DE_S.2, COVID19_up_S.2)
overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

S.2.venn[[7]]$label  <- paste(overlaps$a3, collapse="\n")  
# 
S.2.venn[[8]]$label  <- paste(overlaps$a6, collapse="\n")  
# 
S.2.venn[[9]]$label  <- paste(overlaps$a7, collapse="\n")  
# 
S.2.venn[[10]]$label  <- paste(overlaps$a4, collapse="\n")  
# 
S.2.venn[[11]]$label  <- paste(overlaps$a1, collapse="\n")  
 

grid.newpage()
grid.draw(S.2.venn) # save this as svg image


# Oct 5, 2025 note - opted not to include the all diagram in the manuscript
## All datasets combined
## need to get a vectored list 

calu3_DE_all_ = All_in_$Calu3_DE
remove = "-"

calu3_DE_all_ = setdiff(calu3_DE_all_, remove)

lung_DE_all_ = All_in_$Lung_DE
lung_DE_all_ = setdiff(lung_DE_all_, remove)

COVID19_up_all_ = All_in_$COVID19_up
COVID19_up_all_ = setdiff(COVID19_up_all_, remove)

# Chart
All.venn_ = venn.diagram(
  x = list(calu3_DE_all_, lung_DE_all_, COVID19_up_all_),
  category.names = c("calu3 DE", "lung DE", "COVID-19 up"),
  filename = NULL,
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  cex = 0.65,
  fill = myCol3, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main = "Targets of all svRNA-TFOs combined "
  
)

grid.newpage()
grid.draw(All.venn_)

lapply(All.venn_, names)
lapply(All.venn_,function(i) i$label)

x = list(calu3_DE_all_, lung_DE_all_, COVID19_up_all_)
overlaps = calculate.overlap(x)

# Over-write labels (specific labels are chosen by manual check of labels)

All.venn_[[7]]$label  <- paste(overlaps$a1, collapse="\n")  
# 
All.venn_[[8]]$label  <- paste(overlaps$a2, collapse="\n")  
# 
All.venn_[[9]]$label  <- paste(overlaps$a3, collapse="\n")  
# 
All.venn_[[10]]$label  <- paste(overlaps$a4, collapse="\n")  
# 
All.venn_[[11]]$label  <- paste(overlaps$a6, collapse="\n")  
# 
All.venn_[[12]]$label  <- paste(overlaps$a7, collapse="\n")  

grid.newpage()
grid.draw(All.venn_) # save this as svg image



##-------------------------------------------------------------------------
##                          ENrichment Analysis                           -
##-------------------------------------------------------------------------

# Oct 5, 2025,
# downloaded the raw count .tsv for each of the transcriptome studies - GSE148729 (Wyler et al. 2021; Calu3 12 and 24 hpi host libraries only used for transcriptome analysis top 500 up and down each); GSE147507 (Blanco-melo et al. 2020; Calu3 and lung up down top 500 each); and GSE171110 (Lévy et al. 2021; COVID19 up only 1730 non-redundant genes used)
# note Dec 28, 2025 - updates to remove all numeric style gene accession style gene annotations these were not included in the analysis (see supplemental Table S3)
# note Dec 28, 2025 - removed the "-AS1, -AS2, -AS3, -AS4, and -AS5" and "-IT1, IT2, IT3" containing annotations which are for non-coding RNAs ommited from this analysis

# use enrichment analysis to see if matching of 120 predicted target genes (out of 350 Uniprot GeneIDs that could be assigned to the 375 Ensemble ENSG IDs) for the six svRNA-TFOs identified in the study
# hypergeometric distribution used for sampling without replacement

N = 31698 # total combined unique Gene IDs that were measured for infection-relevant differential expression from the three studies mentioned above

K = 4287 # total unique Gene IDs that were DEGs in at least one of the studies (could be up or down)

T = 350 # number of predicted target Gene IDs in this study

O = 120 # number of predicted target Gene IDs in this study that match to a DEG found in at least one of the studies

p_val <- phyper(O-1, K, N-K, T, lower.tail=FALSE) # using one tailed as we assume there would be enrichment based on the hypothesis
p_val #2.73399e-23


## next need to test if same enrichment from the non-svRNA TFO set - can randomly select from the predicted Gene ID targets for twelve remaining TFOs

set.seed(123) 

# Target gene list from the nonsvRNA_TFOs - there are 1,543 of these
nonsvRNA_TFO_geneTargets = read.table(here("data/nonsvRNA_TFO_targets_geneID.txt"), header =F)
target_genes = nonsvRNA_TFO_geneTargets$V1

# Parameters
num_samples <- 1000   # number of times to sample
sample_size <- 350    # number of genes IDs per sample (same as the # samples in the svRNA-TFO set)

# Background DE gene list - 4395 combined GeneIDs = # need to match against this list
background_DEgenes = read.table(here("data/All_DE_4287_bckgrnd_geneID.txt"), header =F)
background_genes = background_DEgenes$V1


# Initialize vector to store results
overlap_counts <- numeric(num_samples)

# Loop
for (i in 1:num_samples) {
  sampled_genes <- sample(target_genes, sample_size, replace = FALSE)
  overlap_counts[i] <- sum(sampled_genes %in% background_genes)
}

# Convert results to data frame
results_df <- data.frame(
  iteration = 1:num_samples,
  overlap_count = overlap_counts
)

# View first few rows
head(results_df)
tail(results_df)

# Optional: summary statistics
summary(results_df$overlap_count)

# -----------------------------
# Compute statistics for the simulated distribution
# -----------------------------
mean_overlap <- mean(results_df$overlap_count) #54.25
var_overlap <- var(results_df$overlap_count)
sd_overlap <- sd(results_df$overlap_count)

se_overlap = sd_overlap/sqrt(1000)

#cat("Mean overlap:", mean_overlap, "\n")
#cat("Variance:", var_overlap, "\n")
#cat("SD:", sd_overlap, "\n")

# -----------------------------
# Step 4: Plot histogram with observed overlap
# -----------------------------

# Observed overlap
observed_overlap <- 120

# Expected overlap
expected_overlap <- round (K/N * T,1) #47.3

# Fold enrichment
fold_enrich <- round(O/expected_overlap, 1) # ~2.5

data <- data.frame(
     Category = c("Expected","Permutat_ctrl","_Observed"),
     DE_genes = c(47.3,54.25, 120),  # Updated expected overlap for T = 350 #updated for mean_overlap
# Add variance or SD only for Paramute_ctrl
SD = c(NA, sd_overlap, NA)  # replace 5.2 with your actual SD for paramute_ctrl
)

# Plot with error bar representing the Standard Error of the Mean for Paramute_ctrl (measures the uncertainty in the estimate of the mean)

ggplot(data, aes(x=Category, y=DE_genes, fill=Category)) +
  geom_bar(stat="identity", width=0.5) +
  geom_errorbar(aes(ymin = DE_genes - SD, ymax = DE_genes + SD),
                width = 0.2, color="black", na.rm = TRUE) +
        geom_text(aes(label=round(DE_genes,1)), vjust=-0.5, size=5) +
        ylim(0, max(data$DE_genes)*1.2) +
        labs(title="***",
             y="Number of DE Genes",
             x="") +
        annotate("text", x=1, y=max(data$DE_genes)*1.15, 
               label=paste0("Fold enrichment ~ ", fold_enrich), size=4) +
        theme_minimal() +
        theme(legend.position="none", text = element_text(size=12))


# -----------------------------
# Compare observed to simulated distribution
# -----------------------------

# Monte Carlo Test
# overlap_counts = vector of 1000 random overlaps
p_val_empirical <- mean(overlap_counts >= observed_overlap)
p_val_empirical
# 0

# Parametric Z-score (assumes normal)
z_score <- (observed_overlap - mean_overlap) / sd_overlap

# same test but compare to the expected
z_score_ <- (expected_overlap - mean_overlap)/ sd_overlap

# One-tailed p-value (enrichment: observed > mean?)
p_val2 <- 1 - pnorm(z_score) #0 p-value significant there is enrichment

# One-tailed p-value (enrichment: expected = mean?)
p_val3 <- 1 - pnorm(z_score_) #0.8799807 no significant difference


cat("Observed overlap:", observed_overlap, "\n") # 120
cat("Mean of random overlaps:", round(mean_overlap,2), "\n") # 54.25
cat("SD of random overlaps:", round(sd_overlap,2), "\n") # 5.92
cat("Z-score:", round(z_score,2), "\n") # 11.12
cat("One-tailed p-value:", signif(p_val2, 4), "\n") # p-value effectively = 0

# graphic not used in paper - shows the observed overlap far from the background that was modelled using parametric test 1000 (i.e., Z-score of 11.11 is very far from mean)
ggplot(results_df, aes(x=overlap_count)) +
  geom_histogram(binwidth=1, fill="skyblue", color="black", alpha=0.7) +
  geom_vline(xintercept = observed_overlap, color="red", linetype="dashed", size=1.2) +
  labs(title="Null distribution of overlaps (random sampling)",
       x="Number of overlaps",
       y="Frequency") +
  annotate("text", x = observed_overlap + 5, y = max(table(overlap_counts))*0.9,
           label = paste0("Observed overlap = ", observed_overlap),
           color = "red", size=5) +
  theme_minimal()


###### end of script

sessionInfo()
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 26100)

# Matrix products: default

#locale:
#  [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8   
# [3] LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=English_Canada.utf8    

# attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] VennDiagram_1.7.3   futile.logger_1.4.3 RColorBrewer_1.1-3  here_1.0.1         
# [5] ggplot2_3.4.4      

