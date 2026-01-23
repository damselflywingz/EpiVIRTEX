## ---------------------------
##
## Script name: 20240308_random_breakpoint_script.R
##
## Purpose: Permutation test for proximity of triplex-forming oligonucleotide in SARS-CoV-2 Wuhan-Hu-1 and recombination breakpoints and hotspots identified by Lytras et al. 2022, https://doi.org/10.1093/gbe/evac018
##
## Author: Amber R. Paulson, Ph.D.
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
## ---------------------------
## Notes:
##   
## Nov 24, 2025 - updates made for the deeper scan results 
##
##
##
## ---------------------------


##-------------------------------------------------------------------------
##                    Part I - ORF-constrained analysis                   -
##-------------------------------------------------------------------------

# !Important Note on implementation: 
# The R Project should be opened with RStudio first, and then next proceed to open and run the script.

############################################ Function to calculate genomic distances
calculate_genomic_distances <- function(annotation_file1, annotation_file2) {
  # Read annotation files
  annotations1 <- file1
  annotations2 <- file2
  
  # Extract chromosome, start, and end positions
  chr1 <- annotations1$chromosome
  start1 <- annotations1$start
  end1 <- annotations1$end
  
  chr2 <- annotations2$chromosome
  start2 <- annotations2$start
  end2 <- annotations2$end
  
  # Calculate distances
  distances <- matrix(NA, nrow = nrow(annotations1), ncol = nrow(annotations2))
  for (i in 1:nrow(annotations1)) {
    for (j in 1:nrow(annotations2)) {
      if (chr1[i] == chr2[j]) {
        distances[i, j] <- min(abs(start1[i] - start2[j]), abs(end1[i] - end2[j]))
      } else {
        distances[i, j] <- NA
      }
    }
  }
  
  return(distances)
} # end of function


# in addition to a genome-wide analysis, each ORF with a svRNA-TFO of concern was assessed for potentail linkage to recombination breakpoint or hotspot

## nsp2

library(here)

nsp2_TFO = read.table(here("data/nsp2.txt"), header=T) #for nsp2 ToC is located from 824 - 835

nsp2.x = 1:1914
nsp2_sample_count = 191 # use 10 % or 1000 max
set.seed(20251124)
nsp2_random = sample(nsp2.x,nsp2_sample_count,replace=FALSE)

nsp2.chrom = c(rep("nsp2",nsp2_sample_count))
nsp2.start = nsp2_random
nsp2.end = nsp2_random
nsp2_ran_annot = data.frame(nsp2.chrom,nsp2.start,nsp2.end)
names(nsp2_ran_annot) = c("chromosome","start","end")
head(nsp2_ran_annot)

file1 = nsp2_TFO
file2 = nsp2_ran_annot

nsp2_191_distances = calculate_genomic_distances(file1, file2)
#print(nsp2_1000_distances)
mean(nsp2_191_distances) #476.0314

head(nsp2_ran_annot)
nsp2_df = data.frame(nsp2_191_distances)
nsp2_df2 =data.frame(Distances = c(t(nsp2_df)), stringsAsFactors=FALSE)

nsp2_df3 = data.frame(nsp2.chrom, nsp2_df2)
names(nsp2_df3) = c("TFO","distance")
head(nsp2_df3)

nsp2_df3_norm =nsp2_df3$distance/1.914

nsp2_df4 = nsp2_df3
nsp2_df4$distance = nsp2_df3_norm
mean(nsp2_df4$distance) #248.7102

## nsp3

# nsp3 has four breakpoints so need to modify the function to identify the shortest distance against four simulated breakpoints

nsp3.x = 1:5835
nsp3_sample_count = 584

nsp3_ran_annot = data.frame()

for (i in 1:nsp3_sample_count){
  output = sample(nsp3.x,4,replace=FALSE)
  nsp3_ran_annot = rbind(nsp3_ran_annot, output)
}

head(nsp3_ran_annot) # when sampling is completed, the col.names are randomly generated

nsp3_ran_annot_combine = data.frame(start = c(nsp3_ran_annot$X582L,nsp3_ran_annot$X4150L, nsp3_ran_annot$X2759L, nsp3_ran_annot$X1181L))
head(nsp3_ran_annot_combine)
dim(nsp3_ran_annot_combine)


nsp3.chrom = c(rep("nsp3",2336))
nsp3.start = nsp3_ran_annot_combine
nsp3.end = nsp3_ran_annot_combine
nsp3_ran_annot_combine2 = data.frame(nsp3.chrom,nsp3.start,nsp3.end)
names(nsp3_ran_annot_combine2) = c("chromosome","start","end")
head(nsp3_ran_annot_combine2) # 2336 locations 

# the nsp TFO region
nsp3_TFO = read.table(here("data/nsp3_TFO.txt"),header=T)

file1 = nsp3_TFO
file2 = nsp3_ran_annot_combine2

nsp3_2336_distances = calculate_genomic_distances(file1, file2) # use the above function to calculate the distances for every breakpoint
#print(nsp3_4000_distances)
#mean(nsp3_4000_distances) 

nsp3_df_ = data.frame(nsp3_2336_distances)
nsp3_df2_ =data.frame(Distances = c(t(nsp3_df_)), stringsAsFactors=FALSE)

nsp3_df3 = data.frame(nsp3.chrom, nsp3_df2_)
names(nsp3_df3) = c("TFO","distance")
head(nsp3_df3)

##################### Function to calculate minimum every four rows
calculate_min_every_four_rows <- function(df, column_name) {
  # Initialize new dataframe
  new_df <- data.frame(minimum = numeric())
  
  # Calculate minimum every four rows
  for (i in seq(1, nrow(df), by = 4)) {
    min_value <- min(df[[column_name]][i:(i + 3)])
    new_df <- rbind(new_df, data.frame(minimum = min_value))
  }
  
  row.names(new_df) <- NULL  # Reset row names
  
  return(new_df)
} #end function

# Usage 
column_name <- "distance"
nsp3_min_df <- calculate_min_every_four_rows(nsp3_df3, column_name)
head(nsp3_df3)
head(nsp3_min_df) #works

nsp3.chrom2 = c(rep("nsp3",584))
nsp3_df4 = data.frame(nsp3.chrom2, nsp3_min_df)
names(nsp3_df4) = c("TFO","distance")
head(nsp3_df4) # works
mean(nsp3_df4$distance) # 904.1421

nsp3_df4_norm =nsp3_df4$distance/5.835
nsp3_df5 = nsp3_df4
nsp3_df5$distance = nsp3_df4_norm
mean(nsp3_df5$distance) # 154.9515

## spike (S-ORF) has 4 breakpoints and 3822 bp + 1 hotspot region of 475 bp - need to combine the codes from N_3UTR and nsp3 above

# Nov 24, 2025 - spike - now with deeper scan has two ToCs in a region with 4 breakpoints, and 2822 bp + 1 hotspot region of 475 bp
# approach will be to half the number or random and combine for s.1 and s.2

s.x = 1:3822
s_sample_count = 382 # use 10 % or 1000 max
s_sample_count_half = 382/2

# need to repeat the below twice
set.seed(20251124)
s_ran_annot.1 = data.frame()

for (i in 1:s_sample_count_half){
  output = sample(s.x,4,replace=FALSE)
  s_ran_annot.1 = rbind(s_ran_annot.1, output)
}

head(s_ran_annot.1) # replace random col.names if required

s_ran_annot.1_combine = data.frame(location = c(s_ran_annot.1$X2295L,s_ran_annot.1$X557L,s_ran_annot.1$X3590L,s_ran_annot.1$X146L))
head(s_ran_annot.1_combine)
dim(s_ran_annot.1_combine)

s.chrom = c(rep("S",764))
s.start = s_ran_annot.1_combine
s.end = s_ran_annot.1_combine
s_ran_annot.1_combine2 = data.frame(s.chrom,s.start,s.end)
names(s_ran_annot.1_combine2) = c("chromosome","start","end")
head(s_ran_annot.1_combine2)

# add a recombination hot spot region after every fourth row

# Define the size of the larger region
larger_region_size <- 3822

# Define the size of the region to select
selected_region_size <- 475

# Generate a random start position within the larger region
set.seed(20251124)
random_start <- sample(1:(larger_region_size - selected_region_size + 1), 1)

# Calculate the end position of the selected region
random_end <- random_start + selected_region_size - 1

# Print the randomly selected region
#cat("Randomly selected region: ", random_start, "-", random_end, "\n")

# set a dataframe
s_ran_annot.1_hs = data.frame()


# run 382 loops, sample (191 for the 2 spike analysis Nov 24 2025) 

for (i in 1:s_sample_count_half){
  random_start = sample(1:(larger_region_size - selected_region_size + 1), 1, replace = FALSE)
  s_ran_annot.1_hs = rbind(s_ran_annot.1_hs, random_start)
}

head(s_ran_annot.1_hs) #change col.name
srandom.1_start2 = s_ran_annot.1_hs$X557L
srandom.1_end = s_ran_annot.1_hs$X557L + selected_region_size - 1

s_ran_annot.1_hs$end = srandom.1_end

s.chrom2.1 = c(rep("S",382))

s_ran_annot.1_hs_combine2 = data.frame(s.chrom2.1,srandom.1_start2,srandom.1_end)
names(s_ran_annot.1_hs_combine2) = c("chromosome","start","end")

# now need to the combine the 382 hot spots with the 1528 breakpoints
head(s_ran_annot.1_combine2)
head(s_ran_annot.1_hs_combine2)

############################ Function to insert rows from df2 into df1 after every fourth row
insert_rows_every_fourth <- function(df1, df2) {
  # Calculate number of rows to insert
  num_insert_rows <- nrow(df1) %/% 4
  
  # Initialize index for insertion
  insert_index <- seq(4, by = 5, length.out = num_insert_rows)
  
  # Insert rows from df2 into df1
  for (i in seq_along(insert_index)) {
    df1 <- rbind(df1[1:insert_index[i], ], df2[i, ], df1[(insert_index[i] + 1):nrow(df1), ])
  }
  
  row.names(df1) <- NULL  # Reset row names
  
  return(df1)
} #end function

# Usage example
s_ran_annot.1_combine3 <- insert_rows_every_fourth(s_ran_annot.1_combine2, s_ran_annot.1_hs_combine2)

tail(s_ran_annot.1_combine3) # notice the function adds two extra rows, but the rest of the dataframe is correct
s_ran_annot.1_combine4 = s_ran_annot.1_combine3[1:955,]

# first svRNA-TFO_S.1 same as before
spike_TFO = read.table(here("data/spike_TFO.txt"),header=T)

file1 = spike_TFO
file2 = s_ran_annot.1_combine4

s.1_955_distances = calculate_genomic_distances(file1, file2)
#print(s_5000_distances)
#mean(s_5000_distances) 

s.1_df = data.frame(s.1_955_distances)
s.1_df2 =data.frame(distances = c(t(s.1_df)), stringsAsFactors=FALSE)

s.1.chrom2 = c(rep("S",955))
s.1_df3 = data.frame(s.1.chrom2, s.1_df2)
names(s.1_df3) = c("TFO","distance")
head(s.1_df3)


######################## Function to calculate minimum every fifth rows
calculate_min_every_five_rows <- function(df, column_name) {
  # Initialize new dataframe
  new_df <- data.frame(minimum = numeric())
  
  # Calculate minimum every five rows
  for (i in seq(1, nrow(df), by = 5)) {
    min_value <- min(df[[column_name]][i:(i + 4)])
    new_df <- rbind(new_df, data.frame(minimum = min_value))
  }
  
  row.names(new_df) <- NULL  # Reset row names
  
  return(new_df)
} #end function

# Usage 
column_name <- "distance"
spike.1_min_df <- calculate_min_every_five_rows(s.1_df3, column_name)
head(s.1_df3)
head(spike.1_min_df) #works

s.1.chrom2 = c(rep("S",191))
s.1_df4 = data.frame(s.1.chrom2, spike.1_min_df)
names(s.1_df4) = c("TFO","distance")
head(s.1_df4) # works
mean(s.1_df4$distance) # 306.555

s.1_df4_norm =s.1_df4$distance/3.822

s.1_df5 = s.1_df4
s.1_df5$distance = s.1_df4_norm
mean(s.1_df5$distance) # 80.208

## second spike TFO (svRNA-TFO-S.2)

# need to repeat the below twice
set.seed(202511242)
s_ran_annot.2 = data.frame()

for (i in 1:s_sample_count_half){
  output = sample(s.x,4,replace=FALSE)
  s_ran_annot.2 = rbind(s_ran_annot.2, output)
}

head(s_ran_annot.2) #replace random col.names if required

s_ran_annot.2_combine = data.frame(location = c(s_ran_annot.2$X1334L,s_ran_annot.2$X1959L,s_ran_annot.2$X3705L,s_ran_annot.2$X2858L))
head(s_ran_annot.2_combine)
dim(s_ran_annot.2_combine)

s.chrom = c(rep("S",764))
s.start = s_ran_annot.2_combine
s.end = s_ran_annot.2_combine
s_ran_annot.2_combine2 = data.frame(s.chrom,s.start,s.end)
names(s_ran_annot.2_combine2) = c("chromosome","start","end")
head(s_ran_annot.2_combine2)

# add a recombination hot spot region after every fourth row

# Define the size of the larger region
larger_region_size <- 3822

# Define the size of the region to select
selected_region_size <- 475

# Generate a random start position within the larger region
set.seed(202511242)
random_start <- sample(1:(larger_region_size - selected_region_size + 1), 1)

# Calculate the end position of the selected region
random_end <- random_start + selected_region_size - 1

# Print the randomly selected region
#cat("Randomly selected region: ", random_start, "-", random_end, "\n")

# set a dataframe
s_ran_annot.2_hs = data.frame()

# run 382 loops, sample (191 for the 2 spike analysis Nov 24 2025) 

for (i in 1:s_sample_count_half){
  random_start = sample(1:(larger_region_size - selected_region_size + 1), 1, replace = FALSE)
  s_ran_annot.2_hs = rbind(s_ran_annot.2_hs, random_start)
}

head(s_ran_annot.2_hs) #change col.name
srandom.2_start2 = s_ran_annot.2_hs$X1959L
srandom.2_end = s_ran_annot.2_hs$X1959L + selected_region_size - 1

s_ran_annot.2_hs$end = srandom.2_end

s.chrom2.2 = c(rep("S",382))

s_ran_annot.2_hs_combine2 = data.frame(s.chrom2.2,srandom.2_start2,srandom.2_end)
names(s_ran_annot.2_hs_combine2) = c("chromosome","start","end")

# now need to the combine the 382 hot spots with the 1528 breakpoints
head(s_ran_annot.2_combine2)
head(s_ran_annot.2_hs_combine2)

# Usage
s_ran_annot.2_combine3 <- insert_rows_every_fourth(s_ran_annot.2_combine2, s_ran_annot.2_hs_combine2)

tail(s_ran_annot.2_combine3) # notice the function adds two extra rows, but the rest of the dataframe is correct
s_ran_annot.2_combine4 = s_ran_annot.2_combine3[1:955,]

#second svRNA-TFO-S.2
spike2_TFO = read.table(here("data/spike2_TFO.txt"),header=T)

file1 = spike2_TFO
file2 = s_ran_annot.2_combine4

s.2_955_distances = calculate_genomic_distances(file1, file2)
#print(s_5000_distances)
#mean(s_5000_distances) 

s.2_df = data.frame(s.2_955_distances)
s.2_df2 =data.frame(distances = c(t(s.2_df)), stringsAsFactors=FALSE)

s.2.chrom2 = c(rep("S",955))
s.2_df3 = data.frame(s.2.chrom2, s.2_df2)
names(s.2_df3) = c("TFO","distance")
head(s.2_df3)

# Calculate minimum distance every five rows
column_name <- "distance"
spike.2_min_df <- calculate_min_every_five_rows(s.2_df3, column_name)
head(s.2_df3)
head(spike.2_min_df) #works

s.2.chrom2 = c(rep("S",191))
s.2_df4 = data.frame(s.2.chrom2, spike.2_min_df)
names(s.2_df4) = c("TFO","distance")
head(s.2_df4)  works
mean(s.2_df4$distance) # 271.712

s.2_df4_norm =s.2_df4$distance/3.822

s.2_df5 = s.2_df4
s.2_df5$distance = s.2_df4_norm
mean(s.2_df5$distance) # 71.09159

# need to merge the two normalized datasets then calculate the average etc.

spike_TFO_combined_df5 = rbind(s.1_df5,s.2_df5)
# 382
mean(spike_TFO_combined_df5$distance) # 75.64979

# E

E_TFO = read.table(here("data/E_TFO.txt"), header=T) 

set.seed(20251124)
E.x = 1:228
E_random = sample(E.x,23, replace=FALSE) 

E.chrom = c(rep("E",23))
E.start = E_random
E.end = E_random
E_ran_annot = data.frame(E.chrom,E.start,E.end)
names(E_ran_annot) = c("chromosome","start","end")
head(E_ran_annot)

file1 = E_TFO
file2 = E_ran_annot

E_23_distances = calculate_genomic_distances(file1, file2)

mean(E_23_distances) # 71.56522

head(E_ran_annot)

E_df = data.frame(E_23_distances)
E_df2 =data.frame(Distances = c(t(E_df)), stringsAsFactors=FALSE)

E_df3 = data.frame(E.chrom, E_df2)
names(E_df3) = c("TFO","distance")
head(E_df3)

E_df3_norm =E_df3$distance/.228
E_df4 = E_df3
E_df4$distance = E_df3_norm

mean(E_df4$distance) # 313.8825

# the nucleocapside (N-ORF) region is a bit different in there is a 475 bp recombination hot spot within a 1630 bp region (N + 3'UTR)
# TFO of concern located at 1114-1126 in this region

# Define the size of the larger region
larger_region_size <- 1630

# Define the size of the region to select
selected_region_size <- 475

# Generate a random start position within the larger region
set.seed(20251124)
random_start <- sample(1:(larger_region_size - selected_region_size + 1), 1)

# Calculate the end position of the selected region
random_end <- random_start + selected_region_size - 1

# Print the randomly selected region
#cat("Randomly selected region: ", random_start, "-", random_end, "\n")

# set a dataframe

n_ran_annot = data.frame()

# run 163 loops, sample 
set.seed(20251124)

for (i in 1:163){
  random_start = sample(1:(larger_region_size - selected_region_size + 1), 1, replace = FALSE)
  n_ran_annot = rbind(n_ran_annot, random_start)
}

head(n_ran_annot) #change col.name as required

random_start2=n_ran_annot$X247L
random_end <- n_ran_annot$X247L + selected_region_size - 1

n_ran_annot$end = random_end

n.chrom = c(rep("N_3UTR",163))

n_ran_annot_combine2 = data.frame(n.chrom,random_start2,random_end)
names(n_ran_annot_combine2) = c("chromosome","start","end")

N3UTR_TFO = read.table(here("data/N3UTR_TFO.txt"),header=T)

file1 = N3UTR_TFO
file2 = n_ran_annot_combine2

N3UTR_163_distances = calculate_genomic_distances(file1, file2)
#print(N3UTR_1000_distances)
mean(N3UTR_163_distances) #235.6503

N3UTR_df = data.frame(N3UTR_163_distances)
N3UTR_df2 =data.frame(Distances = c(t(N3UTR_df)), stringsAsFactors=FALSE)

N_df3 = data.frame(n.chrom, N3UTR_df2)
names(N_df3) = c("TFO","distance")
head(N_df3)

N_df3_norm =N_df3$distance/1.630
N_df4 = N_df3
N_df4$distance = N_df3_norm
mean(N_df4$distance) # 144.5707

# Actual TFO distances to breakpoints

actual = read.table(here("data/actual_TFO.txt"), header=1)
actual_norm = read.table(here("data/actual_TFO_norm.txt"), header=1)


##-------------------------------------------------------------------------
##              Part I - ORF-constrained statistical analysis             -  
##-------------------------------------------------------------------------

all_rand_dist_norm_2S = rbind(nsp2_df4, nsp3_df5,spike_TFO_combined_df5, E_df4, N_df4, actual_norm)
head(all_rand_dist_norm_2S)
levels(all_rand_dist_norm_2S$TFO)

# need the levels for the statistical analysis below

all_rand_dist_norm_2S$TFO = as.factor(all_rand_dist_norm_2S$TFO)
levels(all_rand_dist_norm_2S$TFO)

#re-order

all_rand_dist_norm_2S$TFO = factor(all_rand_dist_norm_2S$TFO, levels = c("nsp2","nsp3","S","E","N_3UTR","Other_TFO","ToC"))
levels(all_rand_dist_norm_2S$TFO)

write.table(all_rand_dist_norm, here("data/20251124_distances_SCV2_norm.csv"), row.names=FALSE, col.names =TRUE, sep = ",")

P = ggplot(all_rand_dist_norm_2S, aes(x=TFO, y=distance)) +
  geom_violin(bw = 30, scale = "width", fill = "lightblue") + theme_bw()

PP = P + stat_summary(fun.data="mean_cl_normal", geom="point", color="black")+
  stat_summary(fun.data="mean_se", geom="errorbar", color="black", width = 0.2)

# print the numeric values from the stat_summary
ggplot_build(PP)$`data`[[3]]


#library("lattice")

oneway.test(distance ~ TFO, data = all_rand_dist_norm_2S, var.equal=TRUE)

res_aov_norm = aov(distance ~ TFO, data=all_rand_dist_norm_2S)
summary(res_aov_norm)

library(multcomp)

post_test_norm = glht(res_aov_norm,
                      linfct = mcp(TFO = "Tukey"))
summary(post_test_norm)

#kruskal wallace test

kruskal.test(distance~TFO, data=all_rand_dist_norm_2S)
pairwise.wilcox.test(all_rand_dist_norm_2S$distance, all_rand_dist_norm_2S$TFO,
                     p.adjust.method = "BH")

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  all_rand_dist_norm$distance and all_rand_dist_norm$TFO 

#nsp2    nsp3    S       E       N_3UTR  Other_TFO
#nsp3      2.8e-14 -       -       -       -       -                
#S         < 2e-16 6.0e-12 -       -       -       -                
#E         0.42445 0.00154 1.7e-07 -       -       -        
#N_3UTR    7.7e-10 0.42445 6.5e-11 0.00128 -       -        
#Other_TFO 0.45764 0.23369 0.00387 0.42445 0.12928 -        
#ToC       0.00055 0.00937 0.04379 0.00285 0.00333 0.00699  

#P value adjustment method: BH 


##-------------------------------------------------------------------------
##                   Part II - Whole genome analysis                      -
##-------------------------------------------------------------------------

# Based on Lytras et al. 2022 need to model 21 x breakpoints, and 2 x 475 bp hotpsots, calculate nearest for 6 x TFOs, then average across all distances

g.x = 1:29903
g_sample_count = 200 # use 10 % or 1200 max
sed.seed(20251124)
g_ran_annot = data.frame()

for (i in 1:g_sample_count){
  output = sample(g.x,21,replace=FALSE)
  g_ran_annot = rbind(g_ran_annot, output)
}

g_ran_annot
colnames(g_ran_annot) = c("ran1","ran2","ran3","ran4","ran5","ran6","ran7","ran8","ran9","ran10","ran11","ran12","ran13","ran14","ran15","ran16","ran17","ran18","ran19","ran20","ran21")

g_ran_bp_combine = data.frame(location = c(g_ran_annot$ran1,g_ran_annot$ran2,g_ran_annot$ran3,g_ran_annot$ran4,g_ran_annot$ran5,g_ran_annot$ran6,g_ran_annot$ran7,g_ran_annot$ran8,g_ran_annot$ran9,g_ran_annot$ran10,g_ran_annot$ran11,g_ran_annot$ran12,g_ran_annot$ran13,g_ran_annot$ran14,g_ran_annot$ran15,g_ran_annot$ran16,g_ran_annot$ran17,g_ran_annot$ran18,g_ran_annot$ran19,g_ran_annot$ran20,g_ran_annot$ran21))
head(g_ran_bp_combine)
dim(g_ran_bp_combine)

g.bp.chrom = c(rep("g",4200))
g.bp.start = g_ran_bp_combine
g.bp.end = g_ran_bp_combine
g_ran_bp_combine2 = data.frame(g.bp.chrom,g.bp.start,g.bp.end)
names(g_ran_bp_combine2) = c("chromosome","start","end")
head(g_ran_bp_combine2)


### add a recombination hot spot region after every 22nd and 23rd row

# Define the size of the larger region
larger_region_size <- 29903

# Define the size of the region to select
selected_region_size <- 475

# Generate a random start position within the larger region
set.seed(20251124)
random_start <- sample(1:(larger_region_size - selected_region_size + 1), 2, replace=FALSE)

# Calculate the end position of the selected region
random_end <- random_start + selected_region_size - 1

# Print the randomly selected region
#cat("Randomly selected region: ", random_start, "-", random_end, "\n")

# set a dataframe

g_ran_annot_hs = data.frame()

# run 200 loops, sample 
set.seed(20251124)
for (i in 1:g_sample_count){
  random_start = sample(1:(larger_region_size - selected_region_size + 1), 2, replace = FALSE)
  g_ran_annot_hs = rbind(g_ran_annot_hs, random_start)
}

# g_ran_annot_hs
grandom_start2_1 = g_ran_annot_hs$X22775L
grandom_start2_2 = g_ran_annot_hs$X16941L

grandom_end_1 = g_ran_annot_hs$X22775L + selected_region_size - 1
grandom_end_2 = g_ran_annot_hs$X16941L + selected_region_size - 1

g_ran_annot_hs$end_1 = grandom_end_1
g_ran_annot_hs$end_2 = grandom_end_2

colnames(g_ran_annot_hs) = c("start_1","start_2","end_1","end_2")

# the hotspot tables must be separated and then inter-collated every 22nd and 23rd row into the breakpoint bed file

# first hotspots
g.hs.chrom2 = c(rep("g",200))
g.hs.start_1 = g_ran_annot_hs$start_1
g.hs.end_1 = g_ran_annot_hs$end_1
g_hs_22_row = data.frame(g.hs.chrom2,g.hs.start_1,g.hs.end_1)
names(g_hs_22_row) = c("chromosome","start","end")
head(g_hs_22_row)


# second hotspots
g.hs.start_2 = g_ran_annot_hs$start_2
g.hs.end_2 = g_ran_annot_hs$end_2
g_hs_23_row = data.frame(g.hs.chrom2,g.hs.start_2,g.hs.end_2)
names(g_hs_23_row) = c("chromosome","start","end")
head(g_hs_23_row)

# now need to the combine the 1000 x 2 hot spots with the 21,000 breakpoints
head(g_ran_bp_combine2)
head(g_hs_22_row)
head(g_hs_23_row)

######################## Function to insert rows from df2 into df1 after every 21st row
insert_rows_every_21st <- function(df1, df2) {
  # Calculate number of rows to insert
  num_insert_rows <- nrow(df1) %/% 21
  
  # Initialize index for insertion
  insert_index <- seq(21, by = 22, length.out = num_insert_rows)
  
  # Insert rows from df2 into df1
  for (i in seq_along(insert_index)) {
    df1 <- rbind(df1[1:insert_index[i], ], df2[i, ], df1[(insert_index[i] + 1):nrow(df1), ])
  }
  
  row.names(df1) <- NULL  # Reset row names
  
  return(df1)
} #end function

# Usage
g_ran_bp_hs_combine1 <- insert_rows_every_21st(g_ran_bp_combine2, g_hs_22_row)
tail(g_ran_bp_hs_combine1) #notice the function adds two extra rows, but the rest of the dataframe is correct
g_ran_bp_hs_combine2 = g_ran_bp_hs_combine1[1:4400,]

########################## Function to insert rows from df2 into df1 after every 22nd row
insert_rows_every_22nd <- function(df1, df2) {
  # Calculate number of rows to insert
  num_insert_rows <- nrow(df1) %/% 22
  
  # Initialize index for insertion
  insert_index <- seq(22, by = 23, length.out = num_insert_rows)
  
  # Insert rows from df2 into df1
  for (i in seq_along(insert_index)) {
    df1 <- rbind(df1[1:insert_index[i], ], df2[i, ], df1[(insert_index[i] + 1):nrow(df1), ])
  }
  
  row.names(df1) <- NULL  # Reset row names
  
  return(df1)
} #end function

# g_ran_bp_hs_combine3 = data.frame()
g_ran_bp_hs_combine3 <- insert_rows_every_22nd(g_ran_bp_hs_combine2, g_hs_23_row)
tail(g_ran_bp_hs_combine3) # notice the function adds two extra rows, but the rest of the dataframe is correct
g_ran_bp_hs_combine4 = g_ran_bp_hs_combine3[1:4600,]

# read in the genomic locations of the TFOs of interest in the Wuhan-Hu-1 genome.

#library(here)
g_TFO = read.table(here("data/g_TFO.txt"),header=T)

file1 = g_TFO
file2 = g_ran_bp_hs_combine4

g_4600_distances = calculate_genomic_distances(file1, file2)

g_df = data.frame(g_4600_distances)
g_df2 =data.frame(distances = c(t(g_df)), stringsAsFactors=FALSE)

g.chrom3 = c(rep("g",4600)) 

# then can re-label
g_df3 = data.frame(g.chrom3, g_df2)
names(g_df3) = c("TFO","distance")
head(g_df3)


########################## Function to calculate minimum every 23rd rows
calculate_min_every_23_rows <- function(df, column_name) {
  # Initialize new dataframe
  new_df <- data.frame(minimum = numeric())
  
  # Calculate minimum every 23 rows
  for (i in seq(1, nrow(df), by = 23)) {
    min_value <- min(df[[column_name]][i:(i + 22)])
    new_df <- rbind(new_df, data.frame(minimum = min_value))
  }
  
  row.names(new_df) <- NULL  # Reset row names
  
  return(new_df)
} #end function

# Usage 
column_name <- "distance"
g_min_df <- calculate_min_every_23_rows(g_df3, column_name)
head(g_df3)
head(g_min_df) # works

g.chrom4 = c(rep("nsp2",200),rep("nsp3",200),rep("S",400),rep("E",200),rep("N_3UTR",200))
g_df4 = data.frame(g.chrom4, g_min_df)
names(g_df4) = c("TFO","distance")
head(g_df4) # works
mean(g_df4$distance) # 709.925

actual = read.table(here("data/actual_TFO_fullgenome.txt"), header=1)


##-------------------------------------------------------------------------
##            Part II - Genome-wide statistical analysis                  -
##-------------------------------------------------------------------------

# combine all together

g.chrom5 = c(rep("random",1200))
g_df5 = data.frame(g.chrom5, g_min_df)
names(g_df5) = c("TFO","distance")


#statistical analysis of distances generated

g_rand_all2 = rbind(g_df5, actual)
head(g_rand_all2)
levels(g_rand_all2$TFO)


ggplot(g_rand_all2, aes(x=TFO, y=distance)) + 
  geom_boxplot(width=0.1) + theme_bw()

P =ggplot(g_rand_all2, aes(x=TFO, y=distance)) +
  geom_violin() + theme_bw()

P + stat_summary(fun.data="mean_sdl",geom="pointrange", color="red")

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

P + stat_summary(fun.data=data_summary)

# need the levels for the statistical analysis below

g_rand_all2$TFO = as.character(g_rand_all2$TFO)
head(g_rand_all2)

g_rand_all2$TFO = as.factor(g_rand_all2$TFO)
levels(g_rand_all2$TFO)

#re-order
g_rand_all2$TFO = factor(g_rand_all2$TFO, levels = c("random","Other_TFO","ToC"))
levels(g_rand_all2$TFO)


# https://statsandr.com/blog/anova-in-r/

summary(g_rand_all2)
boxplot(distance ~ TFO, data= g_rand_all2)


bp = ggplot(g_rand_all2, aes(x =TFO, y = distance,))+
  geom_boxplot(outlier.size = .01, lwd =0.01, width =.5, fill=c("gray","orange","orange")) +
  xlab("randomized and actual breakpoints/hotspots") + ylab("distance per kbp") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(size=9))


pdf("20240330_boxplot_pdf_genome2.pdf", width=3.5, height=3)
bp
dev.off()

# july 27, 2024 printing out th results and making violin plots

write.table(g_rand_all2, here("data/20240727_breakpoint_results_random.txt"), sep = ",", row.names = FALSE, col.names = TRUE)

#write.table(g_rand_all, "20240322_distances_SCV2_genome.csv", row.names=FALSE, col.names =TRUE, sep = ",")

Q = ggplot(g_rand_all2, aes(x=TFO, y=distance)) +
  geom_violin() + theme_bw()

# from Hmisc, wrapper for ggplot2 to work with summary_stat:  smean.sdl computes the mean plus or minus a constant times the standard deviation. 

Q  + stat_summary(fun.data="mean_sdl",geom="pointrange", color="black")

R = ggplot(g_rand_all2, aes(x=TFO, y=distance))
R + geom_jitter()

P_ = ggplot(g_rand_all2, aes(x=TFO, y=distance)) +
  geom_violin(bw = 100, scale = "width", fill = "lightblue") + theme_bw()

# add the error bar for SEM and report the values
P__2 = P_ + stat_summary(fun.data="mean_cl_normal", geom="point", color="black")+
  stat_summary(fun.data="mean_se", geom="errorbar", color="black", width = 0.2)

# print the numeric values from the stat_summary
ggplot_build(P__2)$`data`[[3]]


oneway.test(distance ~ TFO, data = g_rand_all2, var.equal=TRUE)

res_aov = aov(distance ~ TFO, data=g_rand_all2)
summary(res_aov)


#library(multcomp)

linfct = mcp(TFO = "Tukey")

post_test = glht(res_aov,
                 linfct = mcp(TFO = "Tukey"))
summary(post_test)


#kruskal wallace test

kruskal.test(distance~TFO, data=g_rand_all2)
pairwise.wilcox.test(g_rand_all2$distance, g_rand_all2$TFO,
                     p.adjust.method = "BH")

# end script

#sessionInfo()
#R version 4.2.1 (2022-06-23 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 26100)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C     
#[5] LC_TIME=English_Canada.utf8    

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] lattice_0.20-45 ggplot2_3.4.4   here_1.0.1      multcomp_1.4-25 TH.data_1.1-2   MASS_7.3-57     survival_3.3-1  mvtnorm_1.2-3  