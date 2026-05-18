## ---------------------------
##
## Script name: 20240905_SynPlot2_Align_tree.R
##
## Purpose: This is an R script for translated nucleotide alignment, IQTREE-2 tree-building, tip extraction to generate the "pairs.txt" file for Synplot2, and visualization of the findings combined with small RNA-seq coverage in RPM 
##
## Author: Amber R. Paulson, PhD
##
## Date Created: 2024-12-15
## Github: https://github.com/damselflywingz/EpiVIRTEX/
##
##
## Copyright (c) Amber Paulson, 2025 - governed by CC-BY-NC 4.0 - Attribution-Non-Commercial 4.0 International 
## see https://creativecommons.org/licenses/by-nc/4.0/
##
##
## ---------------------------
##
## !Important Note: 
## 
## The R Project should be opened with RStudio first, and then next proceed to open and run the script.  
## IQTree and Synplot2 tools are implemented on the linux server environment by command line. 
##
##
## ---------------------------


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings) 
library(here) # this works with Rstudio to set path to the location to the RProject file
library(plyr)
detach("package:plyr", unload=TRUE)


###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 1 Alignment                             ###
###                                                                     ###
###########################################################################
###########################################################################

##-------------------------------------------------------------------------
##                    Alignment from fasta                                -
##-------------------------------------------------------------------------

# load DNAStringSet, then convert to RNAStringset and align
# See methods in manuscript for accessions list
# February 19, 2024 - added more sequences (100 total, but there is 1 duplicate - so 99) buidling on Lytras et al. 2021, including early lineage A and B, PCoVs, BtCoVs from Temmam et al. Also included the RecMRCA predicted SARS-CoV-2 ancestor

# N ORF
#Full_NORF_seq3 = readDNAStringSet(here("data/20240218_alignment_frBLAST_full_N_coronaviruses_9_addedTemmam_linAB_RecCA_.txt"))
Full_NORF_seq3_ = readDNAStringSet(here("data/20250105_alignment_frBLAST_full_N_coronaviruses_9_addedTemmam_linAB_RecCA_.txt")) #99

#library(ade4)
u_dna4 <- unique(Full_NORF_seq3_) # 76
u_rna4 <- RNAStringSet(u_dna4) # convert to RNA

# April 19, 2024 - nsp2 ORF - identified that one of the accessions (MT072865.1_nsp2) has ambiguous bases so removed
#Full_nsp2_seq = readDNAStringSet(here("data/20240419_nsp2_full_98seq.fasta")) #noticed that NC_014470.1 was missing compared to other alignments (BM48-31)
Full_nsp2_seq = readDNAStringSet(here("data/20241102_nsp2_full_98seq_CHECKED.fasta")) #98
u_dna4_nsp2 <- unique(Full_nsp2_seq) # 76
u_rna4_nsp2 <- RNAStringSet(u_dna4_nsp2) # convert to RNA


# April 21, 2024 - nsp3 ORF - identified that one of the accessions (MT072865.1) has ambiguous bases so removed
# April 21, 2024 - nsp3 ORF - identified poor alignment with the recMRCA - most recent recreated ancestor - removed from the nsp3 analysis
Full_nsp3_seq = readDNAStringSet(here("data/20240421_nsp3_full_97seq.fasta"))
u_dna4_nsp3 <- unique(Full_nsp3_seq) # 85
u_rna4_nsp3 <- RNAStringSet(u_dna4_nsp3) # convert to RNA

## May 10, 2024 - S - identified that one of the accessions (MT072865.1) has ambiguous bases so removed
## May 10, 2024 - S - when the S-oRF is aligned - KF294457.1, DQ412043.1, DQ648857.1, KJ473815.1, and NC_014470 are quite different then the remaining 93 (creating big gaps and issues, opted to remove and work from most conserved sequences)
Full_S_seq = readDNAStringSet(here("data/20240510_S_full_93seq.fasta"))
u_dna4_S <- unique(Full_S_seq) # 80
u_rna4_S <- RNAStringSet(u_dna4_S) # convert to RNA

## July 27, 2024 - E - 
Full_E_seq = readDNAStringSet(here("data/20240727_E_full_100seq_.fasta"))
u_dna4_E <- unique(Full_E_seq) # 54
u_rna4_E <- RNAStringSet(u_dna4_E) # convert to RNA


#BiocManager::install("DECIPHER")
library(DECIPHER)

# start with unaligned RNAstringset

# Nucleocapsid
AA_u_rna4_alignment = AlignTranslation(u_rna4,
                                       sense="+",
                                       direction ="5' to 3'",
                                       readingFrame = 1,
                                       type = "RNAStringSet",
                                       geneticCode=GENETIC_CODE)

BrowseSeqs(AA_u_rna4_alignment, colWidth = 150, highlight=1) # to view
BrowseSeqs(AA_u_rna4_alignment, colWidth = 150, htmlFile=here("20250125_align_N_Sarbeco_RecCA_RNA.html"), openURL = FALSE, highlight=1) # to save

# nsp2

AA_u_rna4_alignment_nsp2 = AlignTranslation(u_rna4_nsp2,
                                       sense="+",
                                       direction ="5' to 3'",
                                       readingFrame = 1,
                                       type = "RNAStringSet",
                                       geneticCode=GENETIC_CODE)

BrowseSeqs(AA_u_rna4_alignment_nsp2, colWidth = 150, highlight=1) # to view
BrowseSeqs(AA_u_rna4_alignment_nsp2, colWidth = 150, htmlFile=here("20241102_align_nsp2_Sarbeco_RecCA_RNA.html"), openURL = FALSE, highlight=1) # to save


# nsp3

AA_u_rna4_alignment_nsp3 = AlignTranslation(u_rna4_nsp3,
                                            sense="+",
                                            direction ="5' to 3'",
                                            readingFrame = 1,
                                            type = "RNAStringSet",
                                            geneticCode=GENETIC_CODE)

BrowseSeqs(AA_u_rna4_alignment_nsp3, colWidth = 150, highlight=1) # to view
BrowseSeqs(AA_u_rna4_alignment_nsp3, colWidth = 150, htmlFile=here("20240421_align_nsp3_Sarbeco_RecCA_RNA.html"), openURL = FALSE, highlight=1) # to save


# S-ORF

AA_u_rna4_alignment_S = AlignTranslation(u_rna4_S,
                                            sense="+",
                                            direction ="5' to 3'",
                                            readingFrame = 1,
                                            type = "RNAStringSet",
                                            geneticCode=GENETIC_CODE)

BrowseSeqs(AA_u_rna4_alignment_S, colWidth = 150, highlight=1) # to view
BrowseSeqs(AA_u_rna4_alignment_S, colWidth = 150, htmlFile=here("20240510_align_S_Sarbeco_RecCA_RNA.html"), openURL = FALSE, highlight=1) # to save

# E-ORF

AA_u_rna4_alignment_E = AlignTranslation(u_rna4_E,
                                         sense="+",
                                         direction ="5' to 3'",
                                         readingFrame = 1,
                                         type = "RNAStringSet",
                                         geneticCode=GENETIC_CODE)

BrowseSeqs(AA_u_rna4_alignment_E, colWidth = 150, highlight=1) # to view
BrowseSeqs(AA_u_rna4_alignment_E, colWidth = 150, htmlFile=here("20240727_align_E_Sarbeco_RecCA_RNA.html"), openURL = FALSE, highlight=1) # to save



##-------------------------------------------------------------------------
##                        Export alignment                                -
##-------------------------------------------------------------------------


writeXStringSet(AA_u_rna4_alignment_nsp2, file=here("data/aln_nsp2_checked_.fasta"))
writeXStringSet(AA_u_rna4_alignment_nsp3, file=here("data/aln_nsp3.fasta"))
writeXStringSet(AA_u_rna4_alignment_S, file=here("data/aln_S.fasta"))
writeXStringSet(AA_u_rna4_alignment_E, file=here("data/aln_E.fasta"))
writeXStringSet(AA_u_rna4_alignment, file = here("data/20250105_N_RecCA_codon_align_transl_.fasta")) # removed MT072864.1 which has an ambiguous base


## alternatively used ape in the past as well (see below)

#library(ape)

#outfile_alignment_transl = file.path(here(), "20231215_N_codon_align_transl.fasta")
#write.dna(AA_u_rna2_alignment, outfile_alignment_transl, format = "fasta", nbcol=-1, colsep="")

#outfile_alignment_transl = file.path(here(), "20240114_N_codon_align_transl.fasta")
#write.dna(AA_u_rna3_alignment, outfile_alignment_transl, format = "fasta", nbcol=-1, colsep="")

#outfile_alignment_transl_SCV2_rel = file.path(here(), "20231222_N_codon_align_transl_SCV2_rel.fasta")
#write.dna(AA_u_rna2_alignment_SCV2_rel, outfile_alignment_transl_SCV2_rel, format = "fasta", nbcol=-1, colsep="")

#outfile_alignment_transl = file.path(here(), "20240408_N_RecCA_codon_align_transl.fasta")
#write.dna(AA_u_rna4_alignment, outfile_alignment_transl, format = "fasta", nbcol=-1, colsep="")

AA_u_rna2_alignment_tranl_masked = readRNAStringSet(here("20231222_N_codon_align_transl_masked_CORRECT.fasta"))
BrowseSeqs(AA_u_rna2_alignment_tranl_masked, colWidth = 150, highlight=1)

AA_u_rna3_alignment_tranl = readRNAStringSet(here("20240114_N_codon_align_transl.fasta"))
BrowseSeqs(AA_u_rna3_alignment_tranl_masked, colWidth = 150, highlight=1)

AA_u_rna2_alignment_tranl_SCV2_rel = readRNAStringSet(here("20231222_N_codon_align_transl_SCV2_rel.fasta"))
BrowseSeqs(AA_u_rna2_alignment_tranl_SCV2_rel, colWidth = 150, highlight=1)
# for AY274119.3, replaced the final "A" with "-", and for DQ648857, replaced the final "CA" with "--"

# AY264119 and DQ648857 are short on the end need to clip this alignment at 1218 so we can keep all the taxa
# for AY274119.3, replaced the final "A" with "-"; and for DQ648857, replaced the final "CA" with "--" and this was the version imported to the server for Synplot2
# this is a known issue, and is reported by Frith 2014, maker of the package

AA_u_rna4_alignment_tranl = readRNAStringSet(here("20240408_N_RecCA_codon_align_transl_.fasta"))
BrowseSeqs(AA_u_rna4_alignment_tranl, colWidth = 150, highlight=1)

# did the same editing to AY264119 and DQ648857 above


###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 2 Phylogeny                             ###
###                                                                     ###
###########################################################################
###########################################################################

##-------------------------------------------------------------------------
##              Recombination-aware fragmented alignments                -
##-------------------------------------------------------------------------

# nsp2 April 19, 2024
# from here, split the alignment into Region 1 (1-879 of the alignment (878 Lytras breakpoint + 1 for codon); 880 - 1917)

# Nov 2, 2024 - updates to the pipeline - IQTREE2 needs to be read in unique sequences in the alignment to avoid the duplicates
# correct alignment "checked"

Region1_nsp2_seq = readRNAStringSet(here("data/20241108_aln_nsp2_checked2_region1_.fasta")) #76
Region1_nsp2_seq_ = unique(Region1_nsp2_seq) #69

Region2_nsp2_seq = readRNAStringSet(here("data/20241108_aln_nsp2_checked2_region2_.fasta")) #76
Region2_nsp2_seq_ = unique(Region2_nsp2_seq) #71
# need an unrooted tree (see the example tree.pdf), https://www.biostars.org/p/62668/

writeXStringSet(Region1_nsp2_seq_, file=here("data/20241108_aln_nsp2_unique_region1_2IQ.fasta")) #replaced "-" with "_" in header Jan 19, 2025 to visualize with ggmsa
writeXStringSet(Region2_nsp2_seq_, file=here("data/20241108_aln_nsp2_unique_region2_2IQ.fasta"))

# nsp3 April 25, 2024
# from here, split alignment into five regions based on the Lytras breakpoints, but sticking to the Synplot requirement of the alignments being a length that is divisible by 3
# recMRCA nsp3 could not be aligned, so it was aligned secondarily using Clustal Omega in Geneious

#Region1_nsp3_seq = readRNAStringSet(here("data/aln_nsp3_region1_recMRCA.fasta")) #pulled error because of the gaps
Region1_nsp3_seq = readBStringSet(here("data/20241110_aln_nsp3_region1_1-384recMRCA_use__.fasta")) #86
Region1_nsp3_seq_ = unique(Region1_nsp3_seq) #66

Region2_nsp3_seq = readBStringSet(here("data/20241111_aln_nsp3_region2_385-984recMRCA_use.fasta")) #86
Region2_nsp3_seq_ = unique(Region2_nsp3_seq) #64

Region3_nsp3_seq = readBStringSet(here("data/aln_nsp3_region3_985-2356_recMRCA_use.fasta")) #86
Region4_nsp3_seq = readBStringSet(here("data/aln_nsp3_region4_2357-5592_recMRCA_use.fasta")) #86
Region5_nsp3_seq = readBStringSet(here("data/aln_nsp3_region5_5593-5937_recMRCA_use.fasta")) #86


Region3_nsp3_seq_ = unique(Region3_nsp3_seq) #72
Region4_nsp3_seq_ = unique(Region4_nsp3_seq) #80
Region5_nsp3_seq_ = unique(Region5_nsp3_seq) #57

writeXStringSet(Region1_nsp3_seq_, file=here("data/20241108_aln_nsp3_unique_region1_use.fasta"))
writeXStringSet(Region2_nsp3_seq_, file=here("data/20241111_aln_nsp3_unique_region2.fasta"))
writeXStringSet(Region3_nsp3_seq_, file=here("data/20241102_aln_nsp3_unique_region3.fasta"))
writeXStringSet(Region4_nsp3_seq_, file=here("data/20241102_aln_nsp3_unique_region4.fasta"))
writeXStringSet(Region5_nsp3_seq_, file=here("data/20241102_aln_nsp3_unique_region5.fasta"))


# S-ORF May 10, 2024
# from here, split alignment into five regions based on Lytras breakpoints, but sticking to the Synplot requirement of the alignments being a length that is divisible by 3 by adjusting forward where required

Region1_S_seq = readBStringSet(here("data/aln_S_region1_1-945.fasta")) #80
Region2_S_seq = readBStringSet(here("data/aln_S_region2_946-1884.fasta")) #80
Region3_S_seq = readBStringSet(here("data/aln_S_region3_1885-2637.fasta")) #80
Region4_S_seq = readBStringSet(here("data/aln_S_region4_2638-3336.fasta")) #80
Region5_S_seq = readBStringSet(here("data/aln_S_region5_3337-3876.fasta")) #80

Region1_S_seq_ = unique(Region1_S_seq) #68
Region2_S_seq_ = unique(Region2_S_seq) #64
Region3_S_seq_ = unique(Region3_S_seq) #61
Region4_S_seq_ = unique(Region4_S_seq) #58
Region5_S_seq_ = unique(Region5_S_seq) #63

writeXStringSet(Region1_S_seq_, file=here("data/20241102_aln_S_unique_region1.fasta"))
writeXStringSet(Region2_S_seq_, file=here("data/20241102_aln_S_unique_region2.fasta"))
writeXStringSet(Region3_S_seq_, file=here("data/20241102_aln_S_unique_region3.fasta"))
writeXStringSet(Region4_S_seq_, file=here("data/20241102_aln_S_unique_region4.fasta"))
writeXStringSet(Region5_S_seq_, file=here("data/20241102_aln_S_unique_region5.fasta"))


##-------------------------------------------------------------------------
##                       Import phylogenies and visualize                 -
##-------------------------------------------------------------------------

#Notes:

#Jan 5, 2025 - picked up ran all IQTREE2 across regions to be used for Synplot2 as well as for E, N, and a 146 nt truncated alignment for N-ORF with the hotspot region removed on the 3' end (see Lytras eta l. 2022)

#E July 27, 2024, #Feb 8, 2025 had to update to remove "-" in the headers
Region_E_seq = readBStringSet(here("data/20250208_aln_E_shortheads_.fasta")) # these are already unique

#Region_N_full = readBStringSet(here("data/20240408_N_RecCA_codon_align_transl_.fasta")) # need to replace all "-" with "_" in the fasta headers to work with ggmsa below
Region_N_full = readBStringSet(here("data/20250105_N_RecCA_codon_align_transl_.fasta")) # need to replace all "-" with "_" in the fasta headers to work with ggmsa below

Region_N_trunc = readBStringSet(here("data/20250105_TRUNCATED_N_RecCA_codon_align_transl_1120.fasta"))

#Nov 2, 2024 - Import from IQTREE 2 version 1.61 

#library(ape)

#### N ORF - Newick tree with support exported from Geneious
# Nov 24, 2024 - for visualization with ggtree to include alignments want to also import the bootstrap values but ape does not support Node labels
# http://blog.phytools.org/2021/08/old-new-function-to-read-nexus-tree.html
# if the contree is imported to Geneious, the can be exported as a newick with support values

library(phytools)

#N
iqtr_N_ORF_full = read.newick(here("data/iqtree2_out/20250105_N_RecCA_codon_align_transl_2R_.fasta.contree"))
iqtr_N_ORF_trunc = read.newick(here("data/iqtree2_out/20250105_N_TRUNCATED_RecCA_codon_align_transl_2R.fasta.contree"))

#nsp2
iqtr_nsp2_ORF_reg1 = read.newick(here("data/iqtree2_out/20250131_aln_nsp2_unique_region1_2R.fasta.contree"))
iqtr_nsp2_ORF_reg2 = read.newick(here("data/iqtree2_out/20241108_aln_nsp2_unique_region2_2IQT_2R.fasta.contree"))

#S
iqtr_s_ORF_reg2 = read.newick(here("data/iqtree2_out/20250202_aln_S_unique_region2_2IQ_FIX.fasta_2R.contree"))
iqtr_s_ORF_reg1 = read.newick(here("data/iqtree2_out/20241102_aln_S_unique_region1_2R.fasta.contree"))
iqtr_s_ORF_reg3 = read.newick(here("data/iqtree2_out/20241102_aln_S_unique_region3_2R.fasta.contree"))
iqtr_s_ORF_reg4 = read.newick(here("data/iqtree2_out/20241102_aln_S_unique_region4_2R.fasta.contree"))
iqtr_s_ORF_reg5 = read.newick(here("data/iqtree2_out/20241102_aln_S_unique_region5_2R.fasta.contree"))

#E_ORF - using the Feb 8, 2025 version, with short headers, exported from Geneious 
iqtr_e_ORF = read.newick(here("data/iqtree2_out/20250208_aln_E_shortheads_.fasta_2R.contree"))

#nsp3
iqtr_nsp3_reg1 = read.newick(here("data/iqtree2_out/20241108_aln_nsp3_unique_region1_use_2R.fasta.contree"))
iqtr_nsp3_reg2 = read.newick(here("data/iqtree2_out/20241111_aln_nsp3_unique_region2_2R.fasta.contree"))
iqtr_nsp3_reg3 = read.newick(here("data/iqtree2_out/20250705_aln_nsp3_unique_region3.fasta_2R.contree"))
iqtr_nsp3_reg4 = read.newick(here("data/iqtree2_out/20250705_aln_nsp3_unique_region4_2R.fasta.contree"))
iqtr_nsp3_reg5 = read.newick(here("data/iqtree2_out/20250705_aln_nsp3_unique_region5_2R.fasta.contree"))


# visualize the unrooted trees for the regions containing TFOs in nsp2, nsp3, S, E, and N, plus the nsp3 region2 which is encodes svRNA_3.2 in SCV1

library(ggtree)
library(ggplot2)
#detach("package:phyloseq",unload=TRUE)

#N
q_iqtr_N_ORF = ggtree(iqtr_N_ORF_full,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_N_ORF + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))

#this works pretty well can get the order, also see:  https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function

q_iqtr_N_ORF_trunc = ggtree(iqtr_N_ORF_trunc,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_N_ORF_trunc + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))

# The N-ORF truncated was generated to see if major changes with the 475 bp hotspot region removed, but not included in the manuscript

#nsp2 region1
q_iqtr_nsp2_ORF_reg1 = ggtree(iqtr_nsp2_ORF_reg1,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_nsp2_ORF_reg1 + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))

#S region2
q_iqtr_s_ORF_reg2 = ggtree(iqtr_s_ORF_reg2,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_s_ORF_reg2 + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))

#E ORF - this was only phylogeny not split according to Lytras et al. 2022 because of short length (228 nt)
q_iqtr_e_ORF = ggtree(iqtr_e_ORF,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_e_ORF + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))

#nsp3 region1
q_iqtr_nsp3_reg1 = ggtree(iqtr_nsp3_reg1,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_nsp3_reg1 + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))

#nsp3 region2
q_iqtr_nsp3_reg2 = ggtree(iqtr_nsp3_reg2,layout="equal_angle") +geom_tiplab(size=.5, aes(angle=angle), offset = -.005) 
q_iqtr_nsp3_reg2 + coord_cartesian(clip = 'off') +
  theme_tree(plot.margin=margin(80,120,50,50))


### visualize the tree for supplemental figure to include short alignments

#ggtree(iqtr_N_ORF_phytool) + geom_tiplab(size=1) + geom_nodelab(size=2, hjust = 2) # works but need cut off for bootstraps
# https://yulab-smu.top/treedata-book/chapter13.html

ggtree(iqtr_N_ORF_full) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p4 = ggtree(iqtr_N_ORF_full) + geom_tiplab(size=1)
p5 <- p4 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

# not included in manuscript - to determine if major changes with the hotspot region removed
ggtree(iqtr_N_ORF_trunc) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p4_t = ggtree(iqtr_N_ORF_trunc) + geom_tiplab(size=1)
p5_t <- p4_t + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)


# something is weird with the tree tip labels as some have single quote and others don't
N_tip_labels = iqtr_N_ORF_full$tip.label
N_tip_labels_ <- gsub("'", '', N_tip_labels)
iqtr_N_ORF_full$tip.label = N_tip_labels_ 

N_tip_labels_t = iqtr_N_ORF_trunc$tip.label
N_tip_labels_t_ <- gsub("'", '', N_tip_labels_t)
iqtr_N_ORF_trunc$tip.label = N_tip_labels_t_ 


# need to root this tree for visualization only

#library(ape)

iqtr_N_ORF_root = root(iqtr_N_ORF_full, outgroup = "KY352407.1_27665_28924_BtCoV_BtKY72",resolve.root = TRUE)
p6 = ggtree(iqtr_N_ORF_root,linewidth = 0.1) + geom_tiplab(size=1)
p7 <- p6 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p8 = p7 + geom_treescale(x=0, fontsize =1)


iqtr_N_ORF_root_t = root(iqtr_N_ORF_trunc, outgroup = "KY352407.1_27665_28924_BtCoV_BtKY72",resolve.root = TRUE)
p6_t = ggtree(iqtr_N_ORF_root_t,linewidth = 0.1) + geom_tiplab(size=1)
p7_t <- p6_t + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p8_t = p7_t + geom_treescale(x=0, fontsize =1)


#nsp2 reg1 for supplement

ggtree(iqtr_nsp2_ORF_reg1) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p9 = ggtree(iqtr_nsp2_ORF_reg1) + geom_tiplab(size=1)
p10 <- p9 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

#fix tip labels
nsp2_reg1_tip_labels = iqtr_nsp2_ORF_reg1$tip.label
nsp2_reg1_tip_labels_ <- gsub("'", '', nsp2_reg1_tip_labels)
iqtr_nsp2_ORF_reg1$tip.label = nsp2_reg1_tip_labels_ 


# root
iqtr_nsp2_ORF_reg1_root = root(iqtr_nsp2_ORF_reg1, outgroup = "KY352407.1_BtKY72",resolve.root = TRUE)
p11 = ggtree(iqtr_nsp2_ORF_reg1_root,linewidth = 0.1) + geom_tiplab(size=1)
p12 <- p11 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p13 = p12 + geom_treescale(x=0, fontsize =1)

#s reg2 for supplement
ggtree(iqtr_s_ORF_reg2) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p9_s = ggtree(iqtr_s_ORF_reg2) + geom_tiplab(size=1)
p10_s <- p9_s + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

#fix tip labels
s_reg2_tip_labels = iqtr_s_ORF_reg2$tip.label
s_reg2_tip_labels_ <- gsub("'", '', s_reg2_tip_labels)
iqtr_s_ORF_reg2$tip.label = s_reg2_tip_labels_ 


# root
iqtr_s_ORF_reg2_root = root(iqtr_s_ORF_reg2, outgroup = "KY352407.1_BtKY72",resolve.root = TRUE)
p11_s = ggtree(iqtr_s_ORF_reg2_root,linewidth = 0.1) + geom_tiplab(size=1)
p12_s <- p11_s + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p13_s = p12_s + geom_treescale(x=0, fontsize =1)


#Dec 27, 2025 update adding S region 3 for supplement
ggtree(iqtr_s_ORF_reg3) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p9a_s = ggtree(iqtr_s_ORF_reg3) + geom_tiplab(size=1)
p10a_s <- p9a_s + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

#fix tip labels
s_reg3_tip_labels = iqtr_s_ORF_reg3$tip.label
s_reg3_tip_labels_ <- gsub("'", '', s_reg3_tip_labels)
iqtr_s_ORF_reg3$tip.label = s_reg3_tip_labels_ 

#root
iqtr_s_ORF_reg3_root = root(iqtr_s_ORF_reg3, outgroup = "KY352407.1_coronavirus_strain_BtKY72",resolve.root = TRUE)
p11a_s = ggtree(iqtr_s_ORF_reg3_root,linewidth = 0.1) + geom_tiplab(size=1)
p12a_s <- p11a_s + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p13a_s = p12a_s + geom_treescale(x=0, fontsize =1)


#E orf for supplement
ggtree(iqtr_e_ORF) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p9_e = ggtree(iqtr_e_ORF) + geom_tiplab(size=1)
p10_e <- p9_e + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

#fix tip labels
e_ORF_tip_labels = iqtr_e_ORF$tip.label
e_ORF_tip_labels_ <- gsub("'", '', e_ORF_tip_labels)
iqtr_e_ORF$tip.label = e_ORF_tip_labels_ 


# root
#iqtr_e_ORF_root = root(iqtr_e_ORF, outgroup = "KY352407.1_BtKY72",resolve.root = TRUE)
#p11_e = ggtree(iqtr_e_ORF_root,linewidth = 0.1) + geom_tiplab(size=1)
#p12_e <- p11_e + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
#p13_e = p12_e + geom_treescale(x=0, fontsize =1)

# trying an alternate roote for E

iqtr_e_ORF_root2 = root(iqtr_e_ORF, outgroup = "MZ081376.1_RmYN05",resolve.root = TRUE)
p11_er2 = ggtree(iqtr_e_ORF_root2,linewidth = 0.1) + geom_tiplab(size=1)
p12_er2 <- p11_er2 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p13_er2 = p12_er2 + geom_treescale(x=0, fontsize =1)


#Feb 9, 2025
#nsp3 reg1 for supplement

ggtree(iqtr_nsp3_reg1) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p9_nsp3_r1 = ggtree(iqtr_nsp3_reg1) + geom_tiplab(size=1)
p10_nsp3_r1 <- p9_nsp3_r1 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

#fix tip labels
nsp3_reg1_tip_labels = iqtr_nsp3_reg1$tip.label
nsp3_reg1_tip_labels_ <- gsub("'", '', nsp3_reg1_tip_labels)
iqtr_nsp3_reg1$tip.label = nsp3_reg1_tip_labels_ 


# root
iqtr_nsp3_reg1_root = root(iqtr_nsp3_reg1, outgroup = "KY352407.1_BtKY72",resolve.root = TRUE)
p11_nsp3_r1 = ggtree(iqtr_nsp3_reg1_root,linewidth = 0.1) + geom_tiplab(size=1)
p12_nsp3_r1 <- p11_nsp3_r1 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p13_nsp3_r1 = p12_nsp3_r1 + geom_treescale(x=0, fontsize =1)


#nsp3 reg2 for supplement

ggtree(iqtr_nsp3_reg2) + geom_tiplab(size=1) + geom_nodelab(size=2,hjust=2)
p9_nsp3_r2 = ggtree(iqtr_nsp3_reg2) + geom_tiplab(size=1)
p10_nsp3_r2 <- p9_nsp3_r2 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 2,color = "black",  nudge_y = 0.6, hjust = 1)

#fix tip labels
nsp3_reg2_tip_labels = iqtr_nsp3_reg2$tip.label
nsp3_reg2_tip_labels_ <- gsub("'", '', nsp3_reg2_tip_labels)
iqtr_nsp3_reg2$tip.label = nsp3_reg2_tip_labels_ 


# root
iqtr_nsp3_reg2_root = root(iqtr_nsp3_reg2, outgroup = "KY352407.1_BtKY72",resolve.root = TRUE)
p11_nsp3_r2 = ggtree(iqtr_nsp3_reg2_root,linewidth = 0.1) + geom_tiplab(size=1)
p12_nsp3_r2 <- p11_nsp3_r2 + geom_nodelab(mapping = aes(x = branch, label = label,subset = !is.na(as.numeric(label)) & as.numeric(label) > 64),size = 1,color = "black",  nudge_y = 0.2, hjust = 1)
p13_nsp3_r2 = p12_nsp3_r2 + geom_treescale(x=0, fontsize =1)


### July 4, 2025

#N region

N_full_align_RNA_ = Region_N_full # see above all "-" in the headers need to be replaced with "_"

#library(ggmsa)
#https://rdrr.io/cran/ggmsa/f/vignettes/ggmsa.Rmd <----

# The ggmsa() output is compatible with geom_facet() and ggtreeExtra::geom_fruit() and can be used to visualize a tree, multiple sequence alignments, and different types of associated data to explore their underlying linkages/associations.
#https://github.com/YuLab-SMU/ggtree/issues/308 could not find
#dd = tidy_msa(N_TFO_align)
dd_long = tidy_msa(N_full_align_RNA_, start = 1054, end = 1080) #from an RNAStringSet
dd_long2 = tidy_msa(N_full_align_RNA_, 1120, 1132) 

#p8 + geom_facet(geom=geom_msa, data=dd, panel = 'TFO', color="Chemistry_AA")
p8 + geom_facet(geom=geom_msa, data=dd_long2, panel = 'Cao et al. hairpin', color="Chemistry_NT")
p9 = p8 + geom_facet(geom=geom_msa, data=dd_long, panel = 'hairpin', color="Chemistry_NT", char_width=0.5) + geom_facet(geom=geom_msa, data=dd_long2, panel = 'TFO', color="Chemistry_NT", char_width=0.5)

p9_ = facet_widths(p9, widths = c(1.5, .34, .16))
# note the Cao et al. hairpin was not included in the submitted version of the manuscript.

#nsp2 Region1 alignment for supplement:

Region_nsp2_reg1_ = readBStringSet(here("data/iqtree2_out/20241108_aln_nsp2_unique_region1_2IQ_.fasta")) # see above all "-" in the headers need to be replaced with "_"
nsp2_TFO_region = tidy_msa(Region_nsp2_reg1_, 821, 832) 
p14= p13 + geom_facet(geom=geom_msa, data=nsp2_TFO_region, panel = 'TFO', color="Chemistry_NT")
#p9 = p8 + geom_facet(geom=geom_msa, data=dd_long, panel = 'hairpin', color="Chemistry_NT", char_width=0.5) + geom_facet(geom=geom_msa, data=dd_long2, panel = 'TFO', color="Chemistry_NT", char_width=0.5)
p14_ = facet_widths(p14, widths = c(1.5, .5))

#S region2 alignment for supplement:

Region_s_reg2_ = readBStringSet(here("data/iqtree2_out/20250202_aln_S_unique_region2_2IQ_FIX.fasta")) # see above all "-" in the headers need to be replaced with "_"
S_TFO_region = tidy_msa(Region_s_reg2_, 876, 887) 
p14_s= p13_s + geom_facet(geom=geom_msa, data=S_TFO_region, panel = 'TFO', color="Chemistry_NT")
#p9 = p8 + geom_facet(geom=geom_msa, data=dd_long, panel = 'hairpin', color="Chemistry_NT", char_width=0.5) + geom_facet(geom=geom_msa, data=dd_long2, panel = 'TFO', color="Chemistry_NT", char_width=0.5)
p14_s_ = facet_widths(p14_s, widths = c(1.5, .5))

#Dec 27, 2025 update adding S region3 alignment for supplement:

Region_s_reg3 = readBStringSet(here("data/iqtree2_out/20241102_aln_S_unique_region3.fasta")) # see above all "-" in the headers need to be replaced with "_"
S_TFO_region3 = tidy_msa(Region_s_reg3, 558, 570) 
p14a_s= p13a_s + geom_facet(geom=geom_msa, data=S_TFO_region3, panel = 'TFO', color="Chemistry_NT")
#p9 = p8 + geom_facet(geom=geom_msa, data=dd_long, panel = 'hairpin', color="Chemistry_NT", char_width=0.5) + geom_facet(geom=geom_msa, data=dd_long2, panel = 'TFO', color="Chemistry_NT", char_width=0.5)
p14a_s_ = facet_widths(p14a_s, widths = c(1.5, .5))

#E ORF

Region_E_seq #54 sequences
E_TFO_region = tidy_msa(Region_E_seq, 52, 63) 
p14_e= p13_e + geom_facet(geom=geom_msa, data=E_TFO_region, panel = 'TFO', color="Chemistry_NT")
p14_e_ = facet_widths(p14_e, widths = c(1.8, .2))


#E ORF alternative root

p14_er2= p13_er2 + geom_facet(geom=geom_msa, data=E_TFO_region, panel = 'TFO', color="Chemistry_NT")
p14_er2_ = facet_widths(p14_er2, widths = c(1.8, .2))

#Feb 8, 2025, April 12, 2025 continued

#nsp3 regions for supplement
Region_nsp3_reg1_ = readBStringSet(here("data/iqtree2_out/20241108_aln_nsp3_unique_region1_use_fix.fasta")) #66 converted "T" to "U" in recMRCA
Region_nsp3_reg1_RNAtovis = readRNAStringSet(here("data/iqtree2_out/20241108_aln_nsp3_unique_region1_use_fix.fasta")) #66 converted "T" to "U" in recMRCA

BrowseSeqs(Region_nsp3_reg1_, openURL = TRUE)

#Region1_nsp3_seq_ #66

nsp3_svRNA3.1_region = tidy_msa(Region_nsp3_reg1_, 334, 351) 
nsp3_TFO_region = tidy_msa(Region_nsp3_reg1_, 358, 369) 
#p14_nsp3_r1 = p13_nsp3_r1 + geom_facet(geom=geom_msa, data=nsp3_svRNA3.1_region, panel = 'svRNA_3.1', color="Chemistry_NT")
p15_nsp3_r1 = p13_nsp3_r1 + geom_facet(geom=geom_msa, data=nsp3_svRNA3.1_region, panel = 'svRNA_3.1', color="Chemistry_NT") + geom_facet(geom=geom_msa, data=nsp3_TFO_region, panel = 'TFO', color="Chemistry_NT")
#p14_nsp3_r1_ = facet_widths(p14_nsp3_r1, widths = c(1.5, .5))
p15_nsp3_r1_ = facet_widths(p15_nsp3_r1, widths = c(1.6, .25, .15))


#Region2_nsp3_seq_ #64 
nsp3_svRNA3.2_region = tidy_msa(Region2_nsp3_seq_, 85, 105) 
#BrowseSeqs(Region2_nsp3_seq_, openURL = TRUE) #SCV2 group has insertion, COVZ45 and COVXC21 (batCoVs) expanded further"GAAGAA"
#p13_nsp3_r2
#p14_nsp3_r2 = p13_nsp3_r2 + geom_facet(geom=geom_msa, data=nsp3_svRNA3.2_region, panel = 'svRNA_3.2', color="Chemistry_NT")
p14_nsp3_r2_ = facet_widths(p14_nsp3_r2, widths = c(1.6, .4))



##-------------------------------------------------------------------------
##              Extract the tip labels to make the pairs.txt              -
##-------------------------------------------------------------------------


####July 5, 2025 - re-extracing tip labels from the IQTREE 2 phylogenies to use for Synplot2 analysis
###### Get tip labels from the unrooted trees, rooted trees are not compatible with Synplot2 and for visualization purposes only

library(dplyr)
library(tidyr) 

# N ORF

iqtr_N_ORF_full

tiplables_iqtr_N_ORF_full = iqtr_N_ORF_full$tip.label

is_tip_iqtr_N_ORF_full = iqtr_N_ORF_full$edge[,2] <= length(iqtr_N_ORF_full$tip.label)
order_tips_iqtr_N_ORF_full <- iqtr_N_ORF_full$edge[is_tip_iqtr_N_ORF_full, 2]

tiblables_correct_iqtr_N_ORF_full = iqtr_N_ORF_full$tip.label[order_tips_iqtr_N_ORF_full]

df_iqtr_N_ORF_full = as.data.frame(tiblables_correct_iqtr_N_ORF_full)

# need to shift column 2 down one row and bring the last row to row #1 so it will go around the phylogeny (see Synplot2 tree.pdf)

# beware rename can sometimes give errors because dplyr is not compatible with plyr see: https://stackoverflow.com/questions/26371279/dplyr-0-3-0-2-rename-idiom-unstable-when-reshape-package-is-loaded
# detach("package:dplyr", character.only = TRUE)

df_iqtr_N_ORF_full = df_iqtr_N_ORF_full %>%
  rename(col1 = tiblables_correct_iqtr_N_ORF_full) #need to rename column to avoid having same object names

df1_iqtr_N_ORF_full = df_iqtr_N_ORF_full %>%
  mutate(col2 = col1)

text_to_add_iqtr_N_ORF_full = tiblables_correct_iqtr_N_ORF_full[1]

df1_iqtr_N_ORF_full = df1_iqtr_N_ORF_full %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df1_iqtr_N_ORF_full = df1_iqtr_N_ORF_full %>% 
  replace_na(list(col2 =text_to_add_iqtr_N_ORF_full)) ## this is correct but need to add ".fasta" and separate by space

df1_iqtr_N_ORF_full = df1_iqtr_N_ORF_full %>%
  mutate(col1 = paste(df1_iqtr_N_ORF_full$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_iqtr_N_ORF_full$col2,".fasta", sep=""))

names(df1_iqtr_N_ORF_full) = NULL

#need to remove the "'" symbol which was added to the tip labels in the nexus format
#https://stackoverflow.com/questions/32857873/remove-a-character-from-the-entire-data-frame - modified

df1_iqtr_N_ORF_full[] <- lapply(df1_iqtr_N_ORF_full, gsub, pattern="'", replacement="")

write.table(df1_iqtr_N_ORF_full, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_iqtr_N_ORF_full.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2

# Note: this needs to be modified to remove "\r" and replace with "" in notepad++ (see the Synplot2 README)
# once the "\r" is removed, export this to the server to test with Synplot2

# April 19, 2024 
# Updated July 5, 2025

# this is for the nsp2

tiplables_nsp2_reg1= iqtr_nsp2_ORF_reg1$tip.label


#fix tip labels nsp2 reg2, reg1 is fixed for visualization
nsp2_reg2_tip_labels = iqtr_nsp2_ORF_reg2$tip.label
nsp2_reg2_tip_labels_ <- gsub("'", '', nsp2_reg2_tip_labels)
iqtr_nsp2_ORF_reg2$tip.label = nsp2_reg2_tip_labels_ 

tiplables_nsp2_reg2= iqtr_nsp2_ORF_reg2$tip.label


is_tip_nsp2_reg1 = iqtr_nsp2_ORF_reg1$edge[,2] <= length(iqtr_nsp2_ORF_reg1$tip.label)
order_tips_nsp2_reg1 <- iqtr_nsp2_ORF_reg1$edge[is_tip_nsp2_reg1, 2]

is_tip_nsp2_reg2 = iqtr_nsp2_ORF_reg2$edge[,2] <= length(iqtr_nsp2_ORF_reg2$tip.label)
order_tips_nsp2_reg2 <- iqtr_nsp2_ORF_reg2$edge[is_tip_nsp2_reg2, 2]

tiblables_correct_nsp2_reg1 = iqtr_nsp2_ORF_reg1$tip.label[order_tips_nsp2_reg1]
tiblables_correct_nsp2_reg2 = iqtr_nsp2_ORF_reg2$tip.label[order_tips_nsp2_reg2]

df_nsp2_reg1 = as.data.frame(tiblables_correct_nsp2_reg1)
df_nsp2_reg2 = as.data.frame(tiblables_correct_nsp2_reg2)

# need to shift column 2 down one row and bring the last row to row #1 so it will go around the phylogeny (see Synplot2 tree.pdf)
# rename can sometimes give errors because dplyr is not compatible with plyr see: https://stackoverflow.com/questions/26371279/dplyr-0-3-0-2-rename-idiom-unstable-when-reshape-package-is-loaded

#detach("package:dplyr", character.only = TRUE)
#library(dplyr)

df_nsp2_reg1 = df_nsp2_reg1 %>%
  rename(col1 = tiblables_correct_nsp2_reg1) #need to rename column to avoid having same object names
df1_nsp2_reg1 = df_nsp2_reg1 %>%
  mutate(col2 = col1)
text_to_add_nsp2_reg1 = tiblables_correct_nsp2_reg1[1]
df1_nsp2_reg1 = df1_nsp2_reg1 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 


df_nsp2_reg2 = df_nsp2_reg2 %>%
  rename(col1 = tiblables_correct_nsp2_reg2) #need to rename column to avoid having same object names

df1_nsp2_reg2 = df_nsp2_reg2 %>%
  mutate(col2 = col1)

text_to_add_nsp2_reg2 = tiblables_correct_nsp2_reg2[1]

df1_nsp2_reg2 = df1_nsp2_reg2 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

#library(tidyr) 
df1_nsp2_reg1 = df1_nsp2_reg1 %>% 
  replace_na(list(col2 =text_to_add_nsp2_reg1)) ## this is correct but need to add ".fasta" and separate by space

df1_nsp2_reg2 = df1_nsp2_reg2 %>% 
  replace_na(list(col2 =text_to_add_nsp2_reg2))


#tiblelables_correct2 = df1$tiblelables_correct2
#df2 = as.data.frame(tiblelables_correct2)

df1_nsp2_reg1 = df1_nsp2_reg1 %>%
  mutate(col1 = paste(df1_nsp2_reg1$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp2_reg1$col2,".fasta", sep=""))
df1_nsp2_reg2 = df1_nsp2_reg2 %>%
  mutate(col1 = paste(df1_nsp2_reg2$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp2_reg2$col2,".fasta", sep=""))

names(df1_nsp2_reg1) = NULL
names(df1_nsp2_reg2) = NULL

head(df1_nsp2_reg1)
tail(df1_nsp2_reg1)
head(df1_nsp2_reg2)
tail(df1_nsp2_reg2)

write.table(df1_nsp2_reg1, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp2_region1.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_nsp2_reg2, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp2_region2.txt"), quote=FALSE, row.names=FALSE,sep=" ")

# Note: this needs to be modified to remove "\r" and replace with "" in notepad++ (see the Synplot2 README)
# once the "\r" is removed, export this to the server to test with Synplot2

# April 25, 2024 - nsp3 five regions
# July 5, 2025 - updated the script

tiplables_nsp3_reg1= iqtr_nsp3_reg1$tip.label
tiplables_nsp3_reg2= iqtr_nsp3_reg2$tip.label

#fix tip labels nsp3 regions 3-5, regions 1 and 2 are already fixed for visualization above
nsp3_reg3_tip_labels = iqtr_nsp3_reg3$tip.label
nsp3_reg3_tip_labels_ <- gsub("'", '', nsp3_reg3_tip_labels)
iqtr_nsp3_reg3$tip.label = nsp3_reg3_tip_labels_ 

nsp3_reg4_tip_labels = iqtr_nsp3_reg4$tip.label
nsp3_reg4_tip_labels_ <- gsub("'", '', nsp3_reg4_tip_labels)
iqtr_nsp3_reg4$tip.label = nsp3_reg4_tip_labels_ 

nsp3_reg5_tip_labels = iqtr_nsp3_reg5$tip.label
nsp3_reg5_tip_labels_ <- gsub("'", '', nsp3_reg5_tip_labels)
iqtr_nsp3_reg5$tip.label = nsp3_reg5_tip_labels_ 


tiplables_nsp3_reg3= iqtr_nsp3_reg3$tip.label
tiplables_nsp3_reg4= iqtr_nsp3_reg4$tip.label
tiplables_nsp3_reg5= iqtr_nsp3_reg5$tip.label


is_tip_nsp3_reg1 = iqtr_nsp3_reg1$edge[,2] <= length(iqtr_nsp3_reg1$tip.label)
order_tips_nsp3_reg1 <- iqtr_nsp3_reg1$edge[is_tip_nsp3_reg1, 2]

is_tip_nsp3_reg2 = iqtr_nsp3_reg2$edge[,2] <= length(iqtr_nsp3_reg2$tip.label)
order_tips_nsp3_reg2 <- iqtr_nsp3_reg2$edge[is_tip_nsp3_reg2, 2]

is_tip_nsp3_reg3 = iqtr_nsp3_reg3$edge[,2] <= length(iqtr_nsp3_reg3$tip.label)
order_tips_nsp3_reg3 <- iqtr_nsp3_reg3$edge[is_tip_nsp3_reg3, 2]

is_tip_nsp3_reg4 = iqtr_nsp3_reg4$edge[,2] <= length(iqtr_nsp3_reg4$tip.label)
order_tips_nsp3_reg4 <- iqtr_nsp3_reg4$edge[is_tip_nsp3_reg4, 2]

is_tip_nsp3_reg5 = iqtr_nsp3_reg5$edge[,2] <= length(iqtr_nsp3_reg5$tip.label)
order_tips_nsp3_reg5 <- iqtr_nsp3_reg5$edge[is_tip_nsp3_reg5, 2]


tiblables_correct_nsp3_reg1 = iqtr_nsp3_reg1$tip.label[order_tips_nsp3_reg1]
tiblables_correct_nsp3_reg2 = iqtr_nsp3_reg2$tip.label[order_tips_nsp3_reg2]
tiblables_correct_nsp3_reg3 = iqtr_nsp3_reg3$tip.label[order_tips_nsp3_reg3]
tiblables_correct_nsp3_reg4 = iqtr_nsp3_reg4$tip.label[order_tips_nsp3_reg4]
tiblables_correct_nsp3_reg5 = iqtr_nsp3_reg5$tip.label[order_tips_nsp3_reg5]


df_nsp3_reg1 = as.data.frame(tiblables_correct_nsp3_reg1)
df_nsp3_reg2 = as.data.frame(tiblables_correct_nsp3_reg2)
df_nsp3_reg3 = as.data.frame(tiblables_correct_nsp3_reg3)
df_nsp3_reg4 = as.data.frame(tiblables_correct_nsp3_reg4)
df_nsp3_reg5 = as.data.frame(tiblables_correct_nsp3_reg5)

# need to shift column 2 down one row and bring the last row to row #1 so it will go around the phylogeny (see Synplot2 tree.pdf)
# rename can sometimes give errors because dplyr is not compatible with plyr see: https://stackoverflow.com/questions/26371279/dplyr-0-3-0-2-rename-idiom-unstable-when-reshape-package-is-loaded

df_nsp3_reg1 = df_nsp3_reg1 %>%
  rename(col1 = tiblables_correct_nsp3_reg1) #need to rename column to avoid having same object names
df1_nsp3_reg1 = df_nsp3_reg1 %>%
  mutate(col2 = col1)
text_to_add_nsp3_reg1 = tiblables_correct_nsp3_reg1[1]
df1_nsp3_reg1 = df1_nsp3_reg1 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_nsp3_reg2 = df_nsp3_reg2 %>%
  rename(col1 = tiblables_correct_nsp3_reg2) #need to rename column to avoid having same object names
df1_nsp3_reg2 = df_nsp3_reg2 %>%
  mutate(col2 = col1)
text_to_add_nsp3_reg2 = tiblables_correct_nsp3_reg2[1]
df1_nsp3_reg2 = df1_nsp3_reg2 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_nsp3_reg3 = df_nsp3_reg3 %>%
  rename(col1 = tiblables_correct_nsp3_reg3) #need to rename column to avoid having same object names
df1_nsp3_reg3 = df_nsp3_reg3 %>%
  mutate(col2 = col1)
text_to_add_nsp3_reg3 = tiblables_correct_nsp3_reg3[1]
df1_nsp3_reg3 = df1_nsp3_reg3 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_nsp3_reg4 = df_nsp3_reg4 %>%
  rename(col1 = tiblables_correct_nsp3_reg4) #need to rename column to avoid having same object names
df1_nsp3_reg4 = df_nsp3_reg4 %>%
  mutate(col2 = col1)
text_to_add_nsp3_reg4 = tiblables_correct_nsp3_reg4[1]
df1_nsp3_reg4 = df1_nsp3_reg4 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_nsp3_reg5 = df_nsp3_reg5 %>%
  rename(col1 = tiblables_correct_nsp3_reg5) #need to rename column to avoid having same object names
df1_nsp3_reg5 = df_nsp3_reg5 %>%
  mutate(col2 = col1)
text_to_add_nsp3_reg5 = tiblables_correct_nsp3_reg5[1]
df1_nsp3_reg5 = df1_nsp3_reg5 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 


#library(tidyr) 
df1_nsp3_reg1 = df1_nsp3_reg1 %>% 
  replace_na(list(col2 =text_to_add_nsp3_reg1)) ## this is correct but need to add ".fasta" and separate by space

df1_nsp3_reg2 = df1_nsp3_reg2 %>% 
  replace_na(list(col2 =text_to_add_nsp3_reg2)) ## this is correct but need to add ".fasta" and separate by space

df1_nsp3_reg3 = df1_nsp3_reg3 %>% 
  replace_na(list(col2 =text_to_add_nsp3_reg3)) ## this is correct but need to add ".fasta" and separate by space

df1_nsp3_reg4 = df1_nsp3_reg4 %>% 
  replace_na(list(col2 =text_to_add_nsp3_reg4)) ## this is correct but need to add ".fasta" and separate by space

df1_nsp3_reg5 = df1_nsp3_reg5 %>% 
  replace_na(list(col2 =text_to_add_nsp3_reg5)) ## this is correct but need to add ".fasta" and separate by space

#tiblelables_correct2 = df1$tiblelables_correct2
#df2 = as.data.frame(tiblelables_correct2)

df1_nsp3_reg1 = df1_nsp3_reg1 %>%
  mutate(col1 = paste(df1_nsp3_reg1$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp3_reg1$col2,".fasta", sep=""))

df1_nsp3_reg2 = df1_nsp3_reg2 %>%
  mutate(col1 = paste(df1_nsp3_reg2$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp3_reg2$col2,".fasta", sep=""))

df1_nsp3_reg3 = df1_nsp3_reg3 %>%
  mutate(col1 = paste(df1_nsp3_reg3$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp3_reg3$col2,".fasta", sep=""))

df1_nsp3_reg4 = df1_nsp3_reg4 %>%
  mutate(col1 = paste(df1_nsp3_reg4$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp3_reg4$col2,".fasta", sep=""))

df1_nsp3_reg5 = df1_nsp3_reg5 %>%
  mutate(col1 = paste(df1_nsp3_reg5$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_nsp3_reg5$col2,".fasta", sep=""))

names(df1_nsp3_reg1) = NULL
names(df1_nsp3_reg2) = NULL
names(df1_nsp3_reg3) = NULL
names(df1_nsp3_reg4) = NULL
names(df1_nsp3_reg5) = NULL

#head(df1_nsp3_reg1)
#tail(df1_nsp3_reg1)

write.table(df1_nsp3_reg1, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp3_region1.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_nsp3_reg2, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp3_region2.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_nsp3_reg3, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp3_region3.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_nsp3_reg4, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp3_region4.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_nsp3_reg5, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_nsp3_region5.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2

# Note: this needs to be modified to remove "\r" and replace with "" in notepad++ (see the Synplot2 README)
# once the "\r" is removed, export this to the server to test with Synplot2


# May 10, 2024 - S ORF five regions
# July 6, 2025 - update script below of IQTREE 2 phylogenies

tiplables_S_reg2= iqtr_s_ORF_reg2$tip.label


#fix tip labels S_reg2 is fixed already above
s_reg1_tip_labels = iqtr_s_ORF_reg1$tip.label
s_reg1_tip_labels_ <- gsub("'", '', s_reg1_tip_labels)
iqtr_s_ORF_reg1$tip.label = s_reg1_tip_labels_ 

s_reg3_tip_labels = iqtr_s_ORF_reg3$tip.label
s_reg3_tip_labels_ <- gsub("'", '', s_reg3_tip_labels)
iqtr_s_ORF_reg3$tip.label = s_reg3_tip_labels_ 

s_reg4_tip_labels = iqtr_s_ORF_reg4$tip.label
s_reg4_tip_labels_ <- gsub("'", '', s_reg4_tip_labels)
iqtr_s_ORF_reg4$tip.label = s_reg4_tip_labels_ 

s_reg5_tip_labels = iqtr_s_ORF_reg5$tip.label
s_reg5_tip_labels_ <- gsub("'", '', s_reg5_tip_labels)
iqtr_s_ORF_reg5$tip.label = s_reg5_tip_labels_ 

tiplables_S_reg1= iqtr_s_ORF_reg1$tip.label
tiplables_S_reg3= iqtr_s_ORF_reg3$tip.label
tiplables_S_reg4= iqtr_s_ORF_reg4$tip.label
tiplables_S_reg5= iqtr_s_ORF_reg5$tip.label


is_tip_S_reg1 = iqtr_s_ORF_reg1$edge[,2] <= length(iqtr_s_ORF_reg1$tip.label)
order_tips_S_reg1 <- iqtr_s_ORF_reg1$edge[is_tip_S_reg1, 2]

is_tip_S_reg2 = iqtr_s_ORF_reg2$edge[,2] <= length(iqtr_s_ORF_reg2$tip.label)
order_tips_S_reg2 <- iqtr_s_ORF_reg2$edge[is_tip_S_reg2, 2]

is_tip_S_reg3 = iqtr_s_ORF_reg3$edge[,2] <= length(iqtr_s_ORF_reg3$tip.label)
order_tips_S_reg3 <- iqtr_s_ORF_reg3$edge[is_tip_S_reg3, 2]

is_tip_S_reg4 = iqtr_s_ORF_reg4$edge[,2] <= length(iqtr_s_ORF_reg4$tip.label)
order_tips_S_reg4 <- iqtr_s_ORF_reg4$edge[is_tip_S_reg4, 2]

is_tip_S_reg5 = iqtr_s_ORF_reg5$edge[,2] <= length(iqtr_s_ORF_reg5$tip.label)
order_tips_S_reg5 <- iqtr_s_ORF_reg5$edge[is_tip_S_reg5, 2]


tiblables_correct_S_reg1 = iqtr_s_ORF_reg1$tip.label[order_tips_S_reg1]
tiblables_correct_S_reg2 = iqtr_s_ORF_reg2$tip.label[order_tips_S_reg2]
tiblables_correct_S_reg3 = iqtr_s_ORF_reg3$tip.label[order_tips_S_reg3]
tiblables_correct_S_reg4 = iqtr_s_ORF_reg4$tip.label[order_tips_S_reg4]
tiblables_correct_S_reg5 = iqtr_s_ORF_reg5$tip.label[order_tips_S_reg5]


df_S_reg1 = as.data.frame(tiblables_correct_S_reg1)
df_S_reg2 = as.data.frame(tiblables_correct_S_reg2)
df_S_reg3 = as.data.frame(tiblables_correct_S_reg3)
df_S_reg4 = as.data.frame(tiblables_correct_S_reg4)
df_S_reg5 = as.data.frame(tiblables_correct_S_reg5)

# need to shift column 2 down one row and bring the last row to row #1 so it will go around the phylogeny (see Synplot2 tree.pdf)
# rename can sometimes give errors because dplyr is not compatible with plyr see: https://stackoverflow.com/questions/26371279/dplyr-0-3-0-2-rename-idiom-unstable-when-reshape-package-is-loaded

df_S_reg1 = df_S_reg1 %>%
  rename(col1 = tiblables_correct_S_reg1) #need to rename column to avoid having same object names
df1_S_reg1 = df_S_reg1 %>%
  mutate(col2 = col1)
text_to_add_S_reg1 = tiblables_correct_S_reg1[1]
df1_S_reg1 = df1_S_reg1 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_S_reg2 = df_S_reg2 %>%
  rename(col1 = tiblables_correct_S_reg2) #need to rename column to avoid having same object names
df1_S_reg2 = df_S_reg2 %>%
  mutate(col2 = col1)
text_to_add_S_reg2 = tiblables_correct_S_reg2[1]
df1_S_reg2 = df1_S_reg2 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_S_reg3 = df_S_reg3 %>%
  rename(col1 = tiblables_correct_S_reg3) #need to rename column to avoid having same object names
df1_S_reg3 = df_S_reg3 %>%
  mutate(col2 = col1)
text_to_add_S_reg3 = tiblables_correct_S_reg3[1]
df1_S_reg3 = df1_S_reg3 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_S_reg4 = df_S_reg4 %>%
  rename(col1 = tiblables_correct_S_reg4) #need to rename column to avoid having same object names
df1_S_reg4 = df_S_reg4 %>%
  mutate(col2 = col1)
text_to_add_S_reg4 = tiblables_correct_S_reg4[1]
df1_S_reg4 = df1_S_reg4 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

df_S_reg5 = df_S_reg5 %>%
  rename(col1 = tiblables_correct_S_reg5) #need to rename column to avoid having same object names
df1_S_reg5 = df_S_reg5 %>%
  mutate(col2 = col1)
text_to_add_S_reg5 = tiblables_correct_S_reg5[1]
df1_S_reg5 = df1_S_reg5 %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

#library(tidyr) 
df1_S_reg1 = df1_S_reg1 %>% 
  replace_na(list(col2 =text_to_add_S_reg1)) ## this is correct but need to add ".fasta" and separate by space

df1_S_reg2 = df1_S_reg2 %>% 
  replace_na(list(col2 =text_to_add_S_reg2)) ## this is correct but need to add ".fasta" and separate by space

df1_S_reg3 = df1_S_reg3 %>% 
  replace_na(list(col2 =text_to_add_S_reg3)) ## this is correct but need to add ".fasta" and separate by space

df1_S_reg4 = df1_S_reg4 %>% 
  replace_na(list(col2 =text_to_add_S_reg4)) ## this is correct but need to add ".fasta" and separate by space

df1_S_reg5 = df1_S_reg5 %>% 
  replace_na(list(col2 =text_to_add_S_reg5)) ## this is correct but need to add ".fasta" and separate by space

#tiblelables_correct2 = df1$tiblelables_correct2
#df2 = as.data.frame(tiblelables_correct2)

df1_S_reg1 = df1_S_reg1 %>%
  mutate(col1 = paste(df1_S_reg1$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_S_reg1$col2,".fasta", sep=""))

df1_S_reg2 = df1_S_reg2 %>%
  mutate(col1 = paste(df1_S_reg2$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_S_reg2$col2,".fasta", sep=""))

df1_S_reg3 = df1_S_reg3 %>%
  mutate(col1 = paste(df1_S_reg3$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_S_reg3$col2,".fasta", sep=""))

df1_S_reg4 = df1_S_reg4 %>%
  mutate(col1 = paste(df1_S_reg4$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_S_reg4$col2,".fasta", sep=""))

df1_S_reg5 = df1_S_reg5 %>%
  mutate(col1 = paste(df1_S_reg5$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_S_reg5$col2,".fasta", sep=""))

names(df1_S_reg1) = NULL
names(df1_S_reg2) = NULL
names(df1_S_reg3) = NULL
names(df1_S_reg4) = NULL
names(df1_S_reg5) = NULL

#head(df1_nsp3_reg1)
#tail(df1_nsp3_reg1)

write.table(df1_S_reg1, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_S_region1.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_S_reg2, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_S_region2.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_S_reg3, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_S_region3.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_S_reg4, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_S_region4.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2
write.table(df1_S_reg5, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_S_region5.txt"), quote=FALSE, row.names=FALSE,sep=" ") # this works but the "\r" must be removed in notepad++ before import to Synplot2

# Note: this needs to be modified to remove "\r" and replace with "" in notepad++ (see the Synplot2 README)
# once the "\r" is removed, export this to the server to test with Synplot2


# E ORF
# July 5, 2025 update for IQTREE 2 

tiplables_E= iqtr_e_ORF$tip.label
#tiplabels have been fixed above for the visualization

is_tip_E = iqtr_e_ORF$edge[,2] <= length(iqtr_e_ORF$tip.label)
order_tips_E <- iqtr_e_ORF$edge[is_tip_E, 2]

tiblables_correct_E = iqtr_e_ORF$tip.label[order_tips_E]

df_E = as.data.frame(tiblables_correct_E)

# need to shift column 2 down one row and bring the last row to row #1 so it will go around the phylogeny (see Synplot2 tree.pdf)
# rename can sometimes give errors because dplyr is not compatible with plyr see: https://stackoverflow.com/questions/26371279/dplyr-0-3-0-2-rename-idiom-unstable-when-reshape-package-is-loaded

#library(dplyr)


df_E = df_E %>%
  rename(col1 = tiblables_correct_E) #need to rename column to avoid having same object names

df1_E = df_E %>%
  mutate(col2 = col1)

text_to_add_E = tiblables_correct_E[1]

df1_E = df1_E %>% 
  mutate(col2 = lead(col2, n = 1, default = NA)) 

#library(tidyr) 
df1_E = df1_E %>% 
  replace_na(list(col2 =text_to_add_E)) ## this is correct but need to add ".fasta" and separate by space

#tiblelables_correct2 = df1$tiblelables_correct2
#df2 = as.data.frame(tiblelables_correct2)

df1_E = df1_E %>%
  mutate(col1 = paste(df1_E$col1,".fasta",sep="")) %>%
  mutate(col2 = paste(df1_E$col2,".fasta", sep=""))

names(df1_E) = NULL

head(df1_E)
tail(df1_E)

write.table(df1_E, here("data/2_Synplot/20250705_Synplot2_input_codon_align_transl_E.txt"), quote=FALSE, row.names=FALSE,sep=" ") 
# this works but the "\r" must be removed in notepad++ before import to Synplot2


###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 3 Analysis by Synplot2                  ###
###                                                                     ###
###########################################################################
###########################################################################

##-------------------------------------------------------------------------
##                    Parse the fasta file                               -
##-------------------------------------------------------------------------

#July 5, 2025 - updated Synplot2 run with N and S for updated for IQTREE-2

# this step is done in command line
# use seqkit command below to parse out the fast file by header (id)
# seqkit split --by-id multi_fasta_file.fasta

# then go into the file where the split fastas are, and run the following to return the file names to match the pairs.txt file
# note the value "50" should be adjusted appropriately to the number of characters in the file name that need to be retained to match the pairs.txt file

# for f in *; do mv "$f" "${f:50}"; done 

# then run Synplot2 (see the README for installation instructions)
# move the pairs file into the "*fasta.split/" directory and run

#synplot2 pairs.txt

# once the program runs, we want the "synplot_plot.txt", which is 1) codon #, 2) observed # substitutions for all pairs, 3) expected # substitutions for all pairs, and 4) standard deviation (estimated from the null)

# July 5, 2025 N-ORF Codon pruning AY274119.3_28120_29335_SCV1_Tor2, final A pruned at 1,219, and DQ_648857.1 final CA at 1228, 1229


##-------------------------------------------------------------------------
##                      Visualize the data                                -
##-------------------------------------------------------------------------


# the Synplot2 README has instructions for making a plot in R - this script includes ggplot2 options and combining with the coverage from small RNAseq dataset
# the caTools package can calculate moving averages

#install.packages("caTools")
library(caTools)

#July 2025 updated for IQTREE-2 and re-do fo the Synplot2

f_new <- read.table(here("data/20250706_synplot2_plot_N_ORF.txt"))
win = 20 

# April 19, 2024 - nsp2 - region 1 and region 2 results need to be joined and positions numbered accordingly
#g1 = read.table(here("data/nsp2_region1_synplot_plot.txt"))
#g2 = read.table(here("data/nsp2_region1_synplot_plot.txt"))
#g_all = rbind(g1,g2)

# August 11, 2025 - update for nsp2 - for IQTREE-2 analysis
g1_new = read.table(here("data/synplot_plot_nsp2_reg1.txt"))
g2_new = read.table(here("data/synplot_plot_nsp2_reg2.txt"))
g_all_new = rbind(g1_new, g2_new)


# April 25, 2024 - nsp3 - regions 1-5 need to be joined and positions numbered accordingly
# Aug 11, 2024 - nsp3 - regions 1-5 need to be joined and positions numbered accordingly, updated for IQTREE-2/Synplot2 results

h1_new = read.table(here("data/synplot_plot_nsp3_reg1.txt"))
h2_new = read.table(here("data/synplot_plot_nsp3_reg2.txt"))
h3_new = read.table(here("data/synplot_plot_nsp3_reg3.txt"))
h4_new = read.table(here("data/synplot_plot_nsp3_reg4.txt"))
h5_new = read.table(here("data/synplot_plot_nsp3_reg5.txt"))
h_all_new = rbind(h1_new, h2_new, h3_new, h4_new, h5_new)

## the dimensions add up - when combining recMRCA it takes away 39 gaps in region2 and 3 gaps in region3. In region 4 we added 3 gaps during hand-pruning, but recMRCA took away 3.
## full nsp3 alignment is 5937 (5835 nucleotides and 102 gaps) - after separating regions and aligning in recMRCA total is 5895 nucleotides (5835 + 60 gaps, after aligning recMRCA and hand-pruning 42 gaps less)


# May 10, 2024 - S - regions 1-5 need to be joined and positions numbered accordingly

# July 5, 2025 - S- regions updated

i1_new = read.table(here("data/synplot_plot_S_reg1.txt"))
i2_new = read.table(here("data/synplot_plot_S_reg2.txt"))
i3_new = read.table(here("data/synplot_plot_S_reg3.txt"))
i4_new = read.table(here("data/synplot_plot_S_reg4.txt"))
i5_new = read.table(here("data/synplot_plot_S_reg5.txt"))
i_all_new = rbind(i1_new, i2_new, i3_new, i4_new, i5_new)

# July 28, 2024 - E

# August 11, 2025 updated for IQTREE-2/Synplot2

j_new = read.table(here("data/synplot_plot_E_reg1.txt"))


# the Cooper et al. paper is per nucleotide and this is per codon, so need to make triplicate each row but numbered sequentially
#July 6, 2025 updated

fff_new = f_new %>% 
  slice(rep(1:n(), each = 3)) #works, need to replace V1 with correct sequential order (3 per codon for nucleotides)

nucl_ = dim(fff_new)[1] # retrieve the number of nucleotides from the dimension of the df
V1_replace_ = 1:nucl_

fff_new['V1'] = V1_replace_

# below is for nsp2

ggg = g_all %>% 
  slice(rep(1:n(), each = 3)) #works, need to replace V1 with correct sequential order (3 per codon for nucleotides)

nucl_ = dim(ggg)[1] # retrieve the number of nucleotides from the dimension of the df
V1_replace_ = 1:nucl_

ggg['V1'] = V1_replace_

## compared to KF294457.1, SCV2 and SCV1 are missing 3 bases at 1771, 1772, 1773 that can be removed to give 1914 nucleotides in the data for plotting
# SCV2
rows_to_remove = c(1771,1772,1773)
ggg_filtered = ggg[-which(ggg$V1 %in% rows_to_remove),] #dim (ggg_filtered) #1914


# below is for nsp2 updated on Aug 11, 2025
#library(tidyr)
#library(dplyr)

#for this alignment the Wuhan and Tor2 first three bases "GCA" are cut off from the alignment so need to add three blank spots, and then remove the additional ones from above
new_row = data.frame(V1 = 1, V2=0, V3=0, V4=0, stringsAsFactors=FALSE)
gg_all_new = rbind(new_row, g_all_new)
head(gg_all_new)

ggg_new_ = gg_all_new %>%
  slice(rep(1:n(), each = 3)) #works, need to replace V1 with correct sequential order (3 per codon for nucleotides)

nucl_ = dim(ggg_new_)[1] # retrieve the number of nucleotides from the dimension of the df
V1_replace_ = 1:nucl_

ggg_new_['V1'] = V1_replace_

## compared to KF294457.1, SCV2 and SCV1 are missing 3 bases at 1771, 1772, 1773 that can be removed to give 1914 nucleotides in the data for plotting
# SCV2
rows_to_remove = c(1771,1772,1773)
ggg_filtered_new = ggg_new_[-which(ggg_new_$V1 %in% rows_to_remove),] #dim (ggg_filtered_new) #1914
ggg_filtered_S1_new = ggg_filtered_new


#below is for nsp3
#updates for Dec 26, 2025 - labeling corrections

#below is for nsp3 updated August 11, 2025 for IQTREE-2/Synplot2

hhh_new = h_all_new %>%
  slice(rep(1:n(), each = 3))

nucl_ = dim(hhh_new)[1] # retrieve the number of nucleotides from the dimension of the df
V1_replace_ = 1:nucl_

hhh_new['V1'] = V1_replace_


## need to adjust the data to align with SCV2, an SCV1 (adjusting for gaps)

# SCV2 nsp3 #the nine bases that were hand-pruned from wuhan hu-1 to get Synplot2 to work do not have to be removed from the genome plot as the bases still exist there
numbers <- c(seq(476, 478, by = 1),seq(1199,1210, by=1),seq(1213,1237,by=1))
rows_to_remove_S2_nsp3 = c(1,2,3,379,380,381,382,383,384,607,608,609,1279,1280,1281,2215,2216,2217,3403,3404, numbers)
hhh_filtered_new = hhh_new[-which(hhh_new$V1 %in% rows_to_remove_S2_nsp3),] #dim (hhh_filtered_new) #5835 - this is the correct length for SCV2 with no gaps

## SCV1 nsp3 #the six bases that were hand-pruned from Tor2 to get Synplot2 to work do not have to be removed from the genome plot

numbers <- c(seq(476, 553, by = 1),seq(1198,1218, by=1),seq(1225,1248,by=1))
rows_to_remove_S1_nsp3 = c(382,383,384,3402,3403,3404,numbers)
hhh_filtered_S1_new = hhh_new[-which(hhh_new$V1 %in% rows_to_remove_S1_nsp3),] #dim (hhh_filtered_S1_new) #5766 - this is the correct length for SCV1 with no gaps

# when analyzing for nucleotides need to multiply by 3
win_nucl=win*3

# May 10 2024 - below is for S ORF
# July 6, 2025 - updates below for S-ORF

iii_new = i_all_new %>%
  slice(rep(1:n(), each = 3))

nucl_ = dim(iii_new)[1] # retrieve the number of nucleotides from the dimension of the df
V1_replace_ = 1:nucl_

iii_new['V1'] = V1_replace_

## need to adjust the data to align with SCV2, an SCV1 (adjusting for gaps)

# SCV2 S-ORF
numbers <- c(seq(1, 24, by = 1),seq(577, 585, by=1),seq(673, 684,by=1),seq(1156, 1158,by=1),seq(1969, 1974,by=1))
rows_to_remove_S2_S = numbers
iii_filtered_new = iii_new[-which(iii_new$V1 %in% rows_to_remove_S2_S),] #dim (iii_filtered) #3822 - this is the correct length for SCV2 with no gaps

## SCV1 S-ORF
numbers <- c(seq(1, 12, by = 1),seq(232, 252, by = 1), seq(463, 474, by = 1), seq(577, 585, by=1), seq(673, 684,by=1), seq(784, 801, by = 1), seq(1156, 1158,by=1), seq(1492, 1494, by = 1), seq(1969, 1974,by=1), seq(2095, 2106, by = 1))
rows_to_remove_S1_S = numbers
iii_filtered_S1_new = iii_new[-which(iii_new$V1 %in% rows_to_remove_S1_S),] #dim (iii_filtered_S1) #3768 - this is the correct length for SCV1 with no gaps

# below is for E ORF
# below is for E for Aug 11, 2025
jjj_new = j_new %>% 
  slice(rep(1:n(), each = 3)) #works, need to replace V1 with correct sequential order (3 per codon for nucleotides)

nucl_ = dim(jjj_new)[1] # retrieve the number of nucleotides from the dimension of the df
V1_replace_ = 1:nucl_

jjj_new['V1'] = V1_replace_

# SCV2 E-ORF
numbers <- c(208, 209, 210)
rows_to_remove_S2_E = numbers
jjj_filtered_new = jjj_new[-which(jjj_new$V1 %in% rows_to_remove_S2_E),] #dim (jjj_filtered_new) #228 - this is the correct length for SCV2 with no gaps

##-------------------------------------------------------------------------
##                  Integrate the small RNA-seq coverage                  -
##-------------------------------------------------------------------------

# Below the small RNA-seq coverage data is brought in to be combined with the Synplot2 Obs/Exp results
# See more detail in "Updated20240121_2022_11_15_script_smallRNA_calu_3_proc_align_quant.R", plotting the Synplot2 results with the 12 h A B replicates from the Wyler et al. 2021 dataset of small RNA-seq

#SCV2 - Wuhan Hu-1 reference genome
A_12h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_12h_A_hist_rpm.txt"),header=F)
B_12h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_12h_B_hist_rpm.txt"),header=F)
colnames(A_12h) = c("chr","position","count")
colnames(B_12h) = c("chr","position","count")

#SCV1 - Tor2 reference genome
A_12h_S1 = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_12h_A_hist_rpm.txt"),header=F)
B_12h_S1 = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_12h_B_hist_rpm.txt"),header=F)
colnames(A_12h_S1) = c("chr","position","count")
colnames(B_12h_S1) = c("chr","position","count")


#library(ggplot2)
#library(scales)
library(dplyr)
#
A_12ha =A_12h %>%
  mutate(replicate = "A")

B_12ha =B_12h %>%
  mutate(replicate = "B")

AB_12h = rbind(A_12ha, B_12ha) #combine rep A and B

A_12ha_S1 =A_12h_S1 %>%
  mutate(replicate = "A")

B_12ha_S1 =B_12h_S1 %>%
  mutate(replicate = "B")

AB_12h_S1 = rbind(A_12ha_S1, B_12ha_S1) #combine rep A and B

# if want to present mean of A and B - noting this was not used in the paper but there is a genome-wide plot of coverage below

#SCV1:

merge_12h_S1 = merge(A_12h_S1, B_12h_S1, by = "position")
col_order = c("chr.x", "position", "count.x","count.y")
merge_12h_S1_ = merge_12h_S1[, col_order]
colnames(merge_12h_S1_)[1] = "chr" 
merge_12h_S1_$mean.count = apply(merge_12h_S1_[,3:4],1,mean)

#SCV2:

merge_12h_S2 = merge(A_12h, B_12h, by = "position")
col_order = c("chr.x", "position", "count.x","count.y")
merge_12h_S2_ = merge_12h_S2[, col_order]
colnames(merge_12h_S2_)[1] = "chr" 
merge_12h_S2_$mean.count = apply(merge_12h_S2_[,3:4],1,mean)

# create the basic genome-wide plots not used in the paper

#S2_12h_mean_plot=ggplot(merge_12h_S2_,aes(x=position,y=mean.count))+
#geom_area(position = 'identity', alpha = 0.7) +
#  theme_bw()

#S1_12h_mean_plot =ggplot(merge_12h_S1_,aes(x=position,y=mean.count))+
#geom_area(position = 'identity', alpha = 0.7) +
#  theme_bw()

# add highlights - # https://stackoverflow.com/questions/32543176/highlight-areas-within-certain-x-range-in-ggplot2

# December 15, 2025 updates to the S1_rects and S2_rects
S1_rects = read.delim(here("data/20251215_SCV1_highlights.txt")) #read in the highlight annotations - newer version see here: G:\My Drive\epigenetic_interference_NOT_SHARED\Draft Manuscript\Paulson_AR_SARS_ncRNA_analysis\data 
S2_rects = read.delim(here("data/20251215_SCV2_highlights.txt")) #read in the highlight annotations - newer version see here: G:\My Drive\epigenetic_interference_NOT_SHARED\Draft Manuscript\Paulson_AR_SARS_ncRNA_analysis\data


library(scales)
colors_type = c("miRNA" = "slateblue1","svRNA_TFO" = "violetred1","HIS_svRNA" = "orange", "TFO" = "grey")

# the below plots are based on mean of the 2 small-RNA seq replicates genome-wide, but not included in the manuscript
S1_avg_plot = ggplot() +
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  geom_area(data = merge_12h_S1_, aes(x=position,y=mean.count), position = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c(colors_type)) +
  labs(title="SARS-CoV-1", y = "Mean RPM", size = 10) +
  ylim(0,5000) +
  theme_bw() +
  theme(text = element_text(size=10)) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = .1))


S2_avg_plot = ggplot() +
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  geom_area(data = merge_12h_S2_, aes(x=position,y=mean.count), position = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c(colors_type)) +
  labs(title="SARS-CoV-2", y = "Mean RPM", size = 10) +
  ylim(0,5000) +
  theme_bw() +
  theme(text = element_text(size=10)) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = .1))

# see closer to the bottom for the multiplot at L1909 below

##-------------------------------------------------------------------------
##        Sub-set for visualizing specific genomic regions                -
##-------------------------------------------------------------------------


# the coverage data can be subset for specific positions, here the N-ORF + 3 UTR is subset

N_3UTR_ = subset(AB_12h, position>28273 & position<29904) #only goes up to 29903
N_3UTR_S1 = subset(AB_12h_S1, position>28119 & position<29741) #only goes up to 29740

N_3UTR_avg = subset(merge_12h_S2_, position>28273 & position<29904) #only goes up to 29903)
N_3UTR_avg_S1 = subset(merge_12h_S1_, position>28119 & position<29741) #only goes up to 29740


# the coverage data can be subset for specific positions, here the nsp2 is subset

nsp2_ = subset(AB_12h, position>805 & position<2720) 
nsp2_S1_ = subset(AB_12h_S1, position>804 & position<2719) 


# the coverage data can be subset for specific positions, here the nsp3 is subset

nsp3_ = subset(AB_12h, position>2719 & position<8555)
nsp3_S1_ = subset(AB_12h_S1, position>2718 & position<8485) 

# the coverage data can be subset for specific positions, here the S-ORF is subset

S_ = subset(AB_12h, position>21562 & position<25385)
S_S1_ = subset(AB_12h_S1, position>21491 & position<25260) 

# the coverage data can be subset for specific positions, here the E-ORF is subset

E_ = subset(AB_12h, position>26244 & position<26473)
E_S1_ = subset(AB_12h_S1, position>26116 & position<26348) 

# N ORF

head(fff_new)
f2_new = fff_new

# nsp2

g2_new = ggg_filtered_new
g2_S1_new = ggg_filtered_S1_new

# nsp3

h2_new = hhh_filtered_new
h2_S1_new = hhh_filtered_S1_new

# S ORF

i2_new = iii_filtered_new
i2_S1_new = iii_filtered_S1_new

# E ORF

j2_new = jjj_filtered_new
j2_S1_new = jjj_new

# converted the V1 to start at the correct genomic locations for the N-ORF
new_V1 = 28274:29491
new_V2 = 28274:29551
new_V3 = 28274:29548

new_V1_nsp2 = 806:2719
new_V1_nsp2_S1 = 805:2718

new_V1_nsp3 = 2720:8554
new_V1_nsp3_S1 = 2719:8484

new_V1_S = 21563:25384
new_V1_S_S1 = 21492:25259

new_V1_E = 26245:26472
new_V1_E_S1 = 26117:26347

d2['V1'] = new_V1
d_m['V1'] = new_V2
e2['V1'] = new_V3

g2_S1_new['V1'] = new_V1_nsp2_S1
g2_new['V1'] = new_V1_nsp2


h2_new['V1'] = new_V1_nsp3
h2_S1_new['V1'] = new_V1_nsp3_S1

i2_new['V1'] = new_V1_S
i2_S1_new['V1'] = new_V1_S_S1

j2_new['V1'] = new_V1_E
j2_S1_new['V1'] = new_V1_E_S1

# create a column in the dataframe that includes the #Exp/obs using the running mean
# see using the caTools package to calculate running averages per window size
# library(caTools)

SSC_NC= runmean(d2$V2,win_nucl)/runmean(d2$V3,win_nucl) #Exp/obs
SSC_NC_= runmean(d_m$V2,win_nucl)/runmean(d_m$V3,win_nucl) #Exp/obs
SSC_NC_sr = runmean(e2$V2,win_nucl)/runmean(e2$V3,win_nucl) #Exp/obs

# SCV2 and SCV1 can use the same Exp/obs data, adjusting for the start and end positions, as well as removing gapped positions from the alignment
SSC_NC_N2 = runmean(f2$V2,win_nucl)/runmean(f2$V3,win_nucl) #Exp/obs
SSC_NC_N2_new = runmean(f2_new$V2,win_nucl)/runmean(f2_new$V3,win_nucl) #Exp/obs

head(SSC_NC)
head(SSC_NC_)
head(SSC_NC_sr)
head(SSC_NC_N2)
head(SSC_NC_N2_new)

d2['V5'] = SSC_NC
tail(d2)

d_m['V5'] = SSC_NC_
tail(d_m)

e2['V5'] = SSC_NC_sr
tail(e2)

f2['V5'] = SSC_NC_N2
tail(f2)

f2_new['V5'] = SSC_NC_N2_new
tail(f2_new)

# nsp2 data

SSC_nsp2_new = runmean(g2_new$V2,win_nucl)/runmean(g2_new$V3,win_nucl)
SSC_nsp2_S1_new = runmean(g2_S1_new$V2,win_nucl)/runmean(g2_S1_new$V3,win_nucl)

g2['V5'] = SSC_nsp2
tail(g2)

g2_S1['V5'] = SSC_nsp2_S1
tail(g2_S1)

g2_new['V5'] = SSC_nsp2_new
tail(g2_new)

g2_S1_new['V5'] = SSC_nsp2_S1_new
tail(g2_S1_new)

# nsp3 data

SSC_nsp3_new = runmean(h2_new$V2,win_nucl)/runmean(h2_new$V3,win_nucl)
SSC_nsp3_S1_new = runmean(h2_S1_new$V2,win_nucl)/runmean(h2_S1_new$V3,win_nucl)

h2_new['V5'] = SSC_nsp3_new
tail(h2_new)

h2_S1_new['V5'] = SSC_nsp3_S1_new
tail(h2_S1_new)

# S ORF data

SSC_S_new = runmean(i2_new$V2,win_nucl)/runmean(i2_new$V3,win_nucl)
SSC_S_S1_new = runmean(i2_S1_new$V2,win_nucl)/runmean(i2_S1_new$V3,win_nucl)

i2_new['V5'] = SSC_S_new
i2_S1_new['V5'] = SSC_S_S1_new

#E ORF data

SSC_E_new = runmean(j2_new$V2,win_nucl)/runmean(j2_new$V3,win_nucl)
SSC_E_S1_new = runmean(j2_S1_new$V2,win_nucl)/runmean(j2_S1_new$V3,win_nucl)

j2_new['V5'] = SSC_E_new
tail(j2_new)

j2_S1_new['V5'] = SSC_E_S1_new
tail(j2_S1_new)

# Next step in visualization is for N ORF

# since we want to align the small RNA-seq count data against the SCV2 and SCV1 reference genomes, need to selectively omit positions 22, 23, 24, 640, 641, 642, 1243, 1244, 1245, 1246, 1247, 1248, 1273, 1274, 1275, 1276, 1277, 1278 not included in the Synplot2

# SCV2
rows_to_remove = c(22, 23, 24, 640, 641, 642, 1243, 1244, 1245, 1246, 1247, 1248,1273,1274,1275,1276,1277,1278)
f2_filtered = f2[-which(f2$V1 %in% rows_to_remove),] #dim (f2_filtered) #1260
f2_new_filtered = f2_new[-which(f2_new$V1 %in% rows_to_remove),] #dim (f2_new_filtered) #1260

# SCV1
numbers <- seq(1278, 1220, by = -1)
rows_to_remove_S1 = c(640, 641, 642, numbers)
f2_filtered_S1 = f2[-which(f2$V1 %in% rows_to_remove_S1),] #dim (f2_filtered_S1) #1216 this is correct # 
f2_new_filtered_S1 = f2_new[-which(f2_new$V1 %in% rows_to_remove_S1),] #dim (f2_filtered_S1) #1216 this is correct #  

# renumber genomic positions relative to SCV2
new_V4 = 28274:29533 
f2_filtered['V1'] = new_V4
f2_new_filtered['V1'] = new_V4

# renumber genomic positions relative to SCV1
new_V4_S1 = 28120:29335
f2_filtered_S1['V1'] = new_V4_S1
f2_new_filtered_S1['V1'] = new_V4_S1


# SCV2:
#in order to bind this into the data set there will need to be 1630 rows, but the 370 last rows have no Synplot2 data

num_rows = 370

data_sheet <- data.frame(
  col1 = rep(NA, num_rows),
  col2 = rep(0, num_rows),
  col3 = rep(0, num_rows),
  col4 = rep(0, num_rows),
  col5 = rep(0, num_rows)
)

data_sheet$col1 <- seq(29534, 29903)
colnames(data_sheet) = c("V1","V2","V3","V4","V5")

f2_mod = rbind(f2_filtered, data_sheet)
f2_new_mod = rbind(f2_new_filtered, data_sheet)
dim(f2_mod) #1630 rows

# SCV1:
# in order to bind this into the small RNA-seq count data there will need to be 1621 rows, but the last 405 rows have no Synplot2 data

num_rows_S1 = 405

data_sheet_S1 <- data.frame(
  col1 = rep(NA, num_rows_S1),
  col2 = rep(0, num_rows_S1),
  col3 = rep(0, num_rows_S1),
  col4 = rep(0, num_rows_S1),
  col5 = rep(0, num_rows_S1)
)

data_sheet_S1$col1 <- seq(29336, 29740)
colnames(data_sheet_S1) = c("V1","V2","V3","V4","V5")

f2_mod_S1 = rbind(f2_filtered_S1, data_sheet_S1)
f2_new_mod_S1 = rbind(f2_new_filtered_S1, data_sheet_S1)
dim(f2_mod_S1) #1621 rows


# below makes a vector with duplicated values because of the A and B replicate
d2_V5_2 = c(d2$V5, d2$V5)
d_m_V5_2 = c(d_m$V5, d_m$V5)
e2_V5_2 = c(e2$V5, e2$V5)
f2_V5_2 = c(f2$V5, f2$V5)
f2_V5_2_mod = c(f2_mod$V5, f2_mod$V5)
f2_V5_2_mod_S1 = c(f2_mod_S1$V5, f2_mod_S1$V5)

# July 6, 2025 update
f2_new_V5_2 = c(f2_new$V5, f2_new$V5)
f2_new_V5_2_mod = c(f2_new_mod$V5, f2_new_mod$V5)
f2_new_V5_2_mod_S1 = c(f2_new_mod_S1$V5, f2_new_mod_S1$V5)

# nsp2
g2_V5_2_new = c(g2_new$V5, g2_new$V5)
g2_V5_2_S1_new = c(g2_S1_new$V5, g2_S1_new$V5)

# nsp3
h2_V5_2_new = c(h2_new$V5, h2_new$V5)
h2_V5_2_S1_new = c(h2_S1_new$V5, h2_S1_new$V5) 

# S ORF
i2_V5_2_new = c(i2_new$V5, i2_new$V5)
i2_V5_2_S1_new = c(i2_S1_new$V5, i2_S1_new$V5) 

# E ORF
j2_V5_2_new = c(j2_new$V5, j2_new$V5)
j2_V5_2_S1_new = c(j2_S1_new$V5, j2_S1_new$V5) 

Nucleo_plus = cbind(Nucleo,d2_V5_2) # this one binds for truncated nucleo region
Nucleo_plus_ = cbind(Nucleo_,d_m_V5_2) # this on binds for full N to be used with masked
Nucleo_plus_sr = cbind(Nucleo_plus_,e2_V5_2) # this one binds the same full N, but includes e2 (SC2-related)

Nucleo_N2_mod = cbind(N_3UTR_,f2_V5_2_mod) #this binds with the small RNA-seq count data, including 3UTR region
Nucleo_N2_mod_S1 = cbind(N_3UTR_S1,f2_V5_2_mod_S1) # 3,242 long N-ORF and 3'UTR region

# July 6, 2025 update
Nucleo_N2_mod_new = cbind(N_3UTR_,f2_new_V5_2_mod) #this binds with the small RNA-seq count data, including 3UTR region
Nucleo_N2_mod_S1_new = cbind(N_3UTR_S1,f2_new_V5_2_mod_S1) # 3,242 long N-ORF and 3'UTR region

# Aug 11, 2025 update
Nsp2_N2_new = cbind(nsp2_,g2_V5_2_new) # 3,828 long nsp2
Nsp2_N2_S1_new = cbind(nsp2_S1_,g2_V5_2_S1_new)

Nsp3_N2_new = cbind(nsp3_,h2_V5_2_new) # 11670 long nsp3
Nsp3_N2_S1_new = cbind(nsp3_S1_,h2_V5_2_S1_new)

S_N2_new = cbind(S_,i2_V5_2_new) # 7644 long S - is 2x 3822
S_N2_S1_new = cbind(S_S1_,i2_V5_2_S1_new)

E_N2_new = cbind(E_,j2_V5_2_new) 
E_N2_S1_new = cbind(E_S1_,j2_V5_2_S1_new)

# scale the Exp/Obs data to include on the second y-axis
#head(Nucleo_plus)
#range(Nucleo_plus$d2_V5_2)

# July 6, 2025 update
Nucleo_N2_mod_new['f2_new_V5_2_mod'] = 1000*(Nucleo_N2_mod_new$f2_new_V5_2_mod)
Nucleo_N2_mod_S1_new['f2_new_V5_2_mod_S1'] = 1000*(Nucleo_N2_mod_S1_new$f2_new_V5_2_mod_S1)

Nsp2_N2_new['g2_V5_2'] = 1000*(Nsp2_N2_new$g2_V5_2)
Nsp2_N2_S1_new['g2_V5_2_S1'] = 1000*(Nsp2_N2_S1_new$g2_V5_2_S1)

Nsp3_N2_new['h2_V5_2'] = 1000*(Nsp3_N2_new$h2_V5_2)
Nsp3_N2_S1_new['h2_V5_2_S1'] = 1000*(Nsp3_N2_S1_new$h2_V5_2_S1) 

S_N2_new['i2_V5_2_new'] = 1000*(S_N2_new$i2_V5_2_new)
S_N2_S1_new['i2_V5_2_S1_new'] = 1000*(S_N2_S1_new$i2_V5_2_S1_new) 

E_N2_new['j2_V5_2_new'] = 1000*(E_N2_new$j2_V5_2_new)
E_N2_S1_new['j2_V5_2_S1_new'] = 1000*(E_N2_S1_new$j2_V5_2_S1_new) 

# initially needed to replace 0s om f2_V5_2 to "" to avoid plotting, but used the alpha option to hide ggplot
# Replace 0 values in the "Value" column with ""
# Nucleo_N2_mod$f2_V5_2_mod <- ifelse(Nucleo_N2_mod$f2_V5_2_mod == 0, "NA", Nucleo_N2_mod$f2_V5_2_mod)

# add comma to the labels in the x-axis, https://stackoverflow.com/questions/74583618/how-to-show-p-value-and-axis-value-correctly-in-r

mycomma <- function(z) {
  out <- scales::label_comma()(z)
  sub("^([0-9]),([0-9]{3})$", "\\1\\2", out)
}

# using the plots below, which include modified legends and the highlighting

Nucleo_12h_plot_N2 = ggplot(Nucleo_N2_mod, aes(x=position))+
  geom_line(aes(x=position,y = f2_V5_2_mod, alpha = I(ifelse(position > 29532, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28250,30000))+
  theme_bw()

PlotNC_bottom_20_N2 = Nucleo_12h_plot_N2 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="bottom")


#library(scales)
#colors_type = c("miRNA" = "slateblue1","svRNA_TFO" = "violetred1","HIS_svRNA" = "orange")

# Nucleocapsid
# SCV2 Rep A and B:
Nucleo_12h_plot_N2 = ggplot(Nucleo_N2_mod, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = f2_V5_2_mod, alpha = I(ifelse(position > 29532, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28274,30000))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

PlotNC_bottom_20_N2 = Nucleo_12h_plot_N2 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

PlotNC_bottom_20_N2_ = PlotNC_bottom_20_N2 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV2 Avg (not included in manuscript went with individual replicates showing)
Nucleo_12h_plot_N2_avg = ggplot(Nucleo_N2_mod_avg, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = f2_V5_2_mod, alpha = I(ifelse(position > 29532, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=mean.count), position = 'identity', alpha = 0.8) +
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28270,29910))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

PlotNC_bottom_20_N2_avg = Nucleo_12h_plot_N2_avg + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

PlotNC_bottom_20_N2_avg_ = PlotNC_bottom_20_N2_avg +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV1 - A and B replicate:
Nucleo_12h_plot_N2_S1 = ggplot(Nucleo_N2_mod_S1, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = f2_V5_2_mod_S1, alpha = I(ifelse(position > 29334, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28115,29745))+
  ggtitle("SARS-CoV") +
  theme_bw()

PlotNC_bottom_20_N2_S1 = Nucleo_12h_plot_N2_S1 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

PlotNC_bottom_20_N2_S1_ = PlotNC_bottom_20_N2_S1 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV1 Avg: (not included in manuscript, went with individual replicates showing)
Nucleo_12h_plot_N2_avg_S1 = ggplot(Nucleo_N2_mod_avg_S1, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = f2_V5_2_mod_S1, alpha = I(ifelse(position > 29334, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=mean.count), position = 'identity', alpha = 0.8) +
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28115,29745))+
  ggtitle("SARS-CoV") +
  theme_bw()

PlotNC_bottom_20_N2_avg_S1 = Nucleo_12h_plot_N2_avg_S1 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

PlotNC_bottom_20_N2_avg_S1_ = PlotNC_bottom_20_N2_avg_S1 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))


##-------------------------------------------------------------------------
##               Combining with the multiple plot function                -
##-------------------------------------------------------------------------


# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# to save the pdf to file - this is used in the supplement
#pdf(here("20250706_multi_S2_S1_mean_genome.pdf"),height=3.5, width=7.5)
multiplot(S2_avg_plot, S1_avg_plot, cols=1) # here is the genome-wide multiplot with average small-RNA seq coverage not included in manuscript
#dev.off()

#library(scales)
#colors_type = c("miRNA" = "slateblue1","svRNA_TFO" = "violetred1","HIS_svRNA" = "orange")


# Nucleocapsid
# SCV2 Rep A and B:
Nucleo_12h_plot_N2_new = ggplot(Nucleo_N2_mod_new, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = f2_new_V5_2_mod, alpha = I(ifelse(position > 29532, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28274,30000))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

PlotNC_bottom_20_N2_new = Nucleo_12h_plot_N2_new + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

PlotNC_bottom_20_N2_new_ = PlotNC_bottom_20_N2_new +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))


# SCV1 - A and B replicate:
Nucleo_12h_plot_N2_S1_new = ggplot(Nucleo_N2_mod_S1_new, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = f2_new_V5_2_mod_S1, alpha = I(ifelse(position > 29334, 0, 1))), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(28115,29745))+
  ggtitle("SARS-CoV-1") +
  theme_bw()

PlotNC_bottom_20_N2_S1_new = Nucleo_12h_plot_N2_S1_new + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

PlotNC_bottom_20_N2_S1_new_ = PlotNC_bottom_20_N2_S1_new +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# to save the pdf to file - ended up presenting the replicates for the zoom plots
#pdf(here("20250706_multi_win20_N_3UTR_S1_S2_new.pdf"),height=4, width=7.5)
#multiplot(PlotNC_bottom_20_N2_new_, PlotNC_bottom_20_N2_S1_new_, cols=1) # this is the N + 3 UTR zoom figure, 76 unique sequences, 
#dev.off()

# Nsp2
# SCV2 Rep A and B:
Nsp2_12h_plot_N2 = ggplot(Nsp2_N2_new, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = g2_V5_2), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(806,2719))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

Plot_nsp2_bottom_20_N2 = Nsp2_12h_plot_N2 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_nsp2_bottom_20_N2 = Plot_nsp2_bottom_20_N2 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV1 - A and B replicate:
Nsp2_12h_plot_N2_S1 = ggplot(Nsp2_N2_S1_new, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = g2_V5_2_S1), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(805,2718))+
  ggtitle("SARS-CoV-1") +
  theme_bw()

Plot_nsp2_bottom_20_N2_S1 = Nsp2_12h_plot_N2_S1 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_nsp2_bottom_20_N2_S1 = Plot_nsp2_bottom_20_N2_S1 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

multiplot(Plot_nsp2_bottom_20_N2, Plot_nsp2_bottom_20_N2_S1, cols=1) 


# Nsp3 - Dec 26, 2025 updated
# SCV2 Rep A and B:
Nsp3_12h_plot_N2 = ggplot(Nsp3_N2_new, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = h2_V5_2), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(2720,8554))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

Plot_nsp3_bottom_20_N2 = Nsp3_12h_plot_N2 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_nsp3_bottom_20_N2 = Plot_nsp3_bottom_20_N2 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV1 - A and B replicate:
Nsp3_12h_plot_N2_S1 = ggplot(Nsp3_N2_S1_new, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = h2_V5_2_S1), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(2719,8484))+
  ggtitle("SARS-CoV-1") +
  theme_bw()

Plot_nsp3_bottom_20_N2_S1 = Nsp3_12h_plot_N2_S1 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_nsp3_bottom_20_N2_S1 = Plot_nsp3_bottom_20_N2_S1 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

multiplot(Plot_nsp3_bottom_20_N2, Plot_nsp3_bottom_20_N2_S1, cols=1) 


# S-ORF # Dec 26, 2025 updated
# SCV2 Rep A and B:
S_12h_plot_N2 = ggplot(S_N2_new, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = i2_V5_2_new), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(21563,25384))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

Plot_S_bottom_20_N2 = S_12h_plot_N2 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_S_bottom_20_N2 = Plot_S_bottom_20_N2 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV1 - A and B replicate:
S_12h_plot_N2_S1 = ggplot(S_N2_S1_new, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = i2_V5_2_S1_new), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(21492,25259))+
  ggtitle("SARS-CoV-1") +
  theme_bw()

Plot_S_bottom_20_N2_S1 = S_12h_plot_N2_S1 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_S_bottom_20_N2_S1 = Plot_S_bottom_20_N2_S1 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# Updated Dec 26, 2025

# E ORF
# SCV2 Rep A and B:
E_12h_plot_N2 = ggplot(E_N2_new, aes(x=position))+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = j2_V5_2_new), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(26245,26472))+
  ggtitle("SARS-CoV-2") +
  theme_bw()

Plot_E_bottom_20_N2 = E_12h_plot_N2 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_E_bottom_20_N2 = Plot_E_bottom_20_N2 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

# SCV1 - A and B replicate:
E_12h_plot_N2_S1 = ggplot(E_N2_S1_new, aes(x=position))+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 2000, fill = Type), alpha=0.5) +
  scale_fill_manual(values = c(colors_type)) +
  geom_line(aes(x=position,y = j2_V5_2_S1_new), linewidth = .2) +
  geom_area(aes(x=position,y=count,colour = replicate, fill= replicate), position = 'identity', alpha=0.1, linewidth=0.2)+
  scale_y_continuous(name = "RPM", labels = scales::comma,sec.axis = sec_axis(~ ./1000, name = "Obs/Exp", labels = scales::comma))+
  scale_x_continuous(labels=mycomma, limits = c(26117,26347))+
  ggtitle("SARS-CoV-1") +
  theme_bw()

Plot_E_bottom_20_N2_S1 = E_12h_plot_N2_S1 + geom_hline(yintercept=1000, linetype="dashed", color="red",linewidth=0.2) + theme(legend.position="top")

Plot_E_bottom_20_N2_S1 = Plot_E_bottom_20_N2_S1 +
  theme(legend.key.size = unit(0.2, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.margin=margin(c(0,0,0,0)))

multiplot(Plot_E_bottom_20_N2, Plot_E_bottom_20_N2_S1, cols=1) 


# End of script

#> sessionInfo() #20231230
#R version 4.2.1 (2022-06-23 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19045)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
#[5] LC_TIME=English_Canada.utf8    

#attached base packages:
#  [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scales_1.2.1        caTools_1.18.2      tidyr_1.3.0         dplyr_1.1.2         ggplot2_3.4.4       ggtree_3.6.2       
#[7] phyloseq_1.42.0     phangorn_2.11.1     ape_5.7-1           DECIPHER_2.26.0     RSQLite_2.3.3       here_1.0.1         
#[13] Biostrings_2.66.0   GenomeInfoDb_1.34.9 XVector_0.38.0      IRanges_2.32.0      S4Vectors_0.36.2    BiocGenerics_0.44.0        
