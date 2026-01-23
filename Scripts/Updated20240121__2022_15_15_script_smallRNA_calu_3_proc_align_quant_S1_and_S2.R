## ------------------------------------------------------------------------
##
## Script name: Updated20240121_2022_11_15_script_smallRNA_calu_3_proc_align_quant.R
##
## Purpose of script: processing, alignment and counting of small RNA sequencing data from SARS-CoV-2, and SARS-CoV-1 infected Calu-3 cells from Wyler et al. 2021 (and also reported in Merino et al. 2020), with some inspiration from Kim et al. 2021.
##
##
## Author: Amber R. Paulson, PhD
##
## Date Created: 2022-11-15
##
## Copyright (c) Amber Paulson, 2025 - governed by CC-BY-NC 4.0 - Attribution-Non-Commercial 4.0 International 
## see https://creativecommons.org/licenses/by-nc/4.0/
##
##
## Email: Amber[dot]Rose[dot]Paulson[at]gmail[dot]com
## Github: https://github.com/damselflywingz/Mod-RNA/
##
## internal note - original script and output is located: D:\Covid_hackthon\2021_mar13_smallRNA_transcriptome
##
## ------------------------------------------------------------------------
##
## Notes: Sections 1 - 4 are in Bash and implemented in a Linux environment
##
##        Sections 5 and above are in R
##
## ------------------------------------------------------------------------



##-------------------------------------------------------------------------
##                1. SRA tool and retrieve fastq files                    -
##-------------------------------------------------------------------------

# SARS-CoV-2-related small RNA-seq libraries from Wyler et al. 2021 were downloaded on January 15, 2022
# SARS-CoV-related small RNA-seq libraries were downloaded on January 21, 2024

# November 15, 2022 must install the SRA tool to get the fastq files

# https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit 

wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -vxzf sratoolkit.tar.gz

export PATH=$HOME/sratoolkit.3.0.1-ubuntu64/bin:$PATH

# may need to add to path, check version if it is not executable
# export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin

# opted to use the most current module when prompted (sra-toolkit/3.0.0 with StdEnv/2020 gcc/9.3.0)

prefetch SRR11550032
fasterq-dump SRR11550032

# downloaded SRR11550025 - SRR11550032 and SRR11550015-SRR11550018 (all but the SARS-CoV-1 libraries) January 15, 2022
# downloaded SRR11550019 - SRR11550024 on January 21, 2024
# consolidated all *.fastq files into /Fastq_files/


##-------------------------------------------------------------------------
##                      2. Flexbar trimming                               -
##-------------------------------------------------------------------------

module load flexbar #using modeules pre-loaded on Digital Research Alliance Canada's Centre for Advanced Computing (Queen's University)
# v3.0.3

# two-pass method provided in Wyler et al. 2021

flexbar -r SRR11550025.fastq -m 15 -q TAIL -qf i1.8 -qt 30 -x 3 -as TGGAATTCTCGGGTGCCAAGG -t Flexbar_out2/SRR11550025_flexbar

# this retained 89 % of input bases and 99 % of reads, for example

# to remove the poly A tail set the minimum overlap to 10, and the error rate to zero to have strict removal only, this makes a lot of the reads too short to keep, there may be reads with 9 polyA or less, but they could be from the 3' UTR region

flexbar -r Flexbar_out2/SRR11550025_flexbar.fastq -m 15 -as AAAAAAAAAA -ao 10 -at 0 -t Flexbar_out3/SRR11550025_flexbar3

#flexbar -r Flexbar_out2/SRR11550024_flexbar.fastq -m 15 -as AAAAAAAAAA -ao 10 -at RIGHT -t Flexbar_out3/SRR11550024_flexbar3 # later version because "-at 0 needs to set TRIM-END default is RIGHT did not end up using this it results in slightly different results from v 3.0.3

# this retains 54 % of the sequences, for example


##-------------------------------------------------------------------------
##                  3. Alignment to reference                             -
##-------------------------------------------------------------------------

module load bowtie2
# v2.4.1 (SCV2), #v2.4.1 (SCV)

# make the bowtie2 index
bowtie2-build -f GCF_009858895.2_ASM985889v3_genomic.fna GCF_009858895.2_ASM985889v3_genomic #SCV2
bowtie2-build -f AY310120.1_SARS_coronavirus_FRA.fa AY310120.1_SARS_coronavirus_FRA

# next step is to align with bowtie2
bowtie2 -t -x GCF_009858895.2_ASM985889v3_genomic -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550025_flexbar3.fastq 20221115_bowtie/calu3_smallRNA_S2_12h_A.bam
# this had 1.78 % alignment

# working inside of the bowtie directory
bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550026_flexbar3.fastq calu3_smallRNA_S2_12h_B.bam
# 2.09 % alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550027_flexbar3.fastq calu3_smallRNA_S2_24h_A.bam
# 2.44% alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550028_flexbar3.fastq calu3_smallRNA_S2_24h_B.bam
# 1.88 % alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550029_flexbar3.fastq calu3_smallRNA_S2_4h_A.bam
# 0.05 % overall alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550030_flexbar3.fastq calu3_smallRNA_S2_4h_B.bam
# 0.05 % overall alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550031_flexbar3.fastq calu3_smallRNA_untr_4h_A.bam
# 0.01 % overall alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550032_flexbar3.fastq calu3_smallRNA_untr_4h_B.bam
# 0.01 % overall alignment

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550015_flexbar3.fastq calu3_smallRNA_mock_24h_A.bam
# 0.01%

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550016_flexbar3.fastq calu3_smallRNA_mock_24h_B.bam
# 0.02%

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550017_flexbar3.fastq calu3_smallRNA_mock_4h_A.bam
# 0.01%

bowtie2 -t -x ../GCF_009858895.2_ASM985889v3_genomic -q ../2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550018_flexbar3.fastq calu3_smallRNA_mock_4h_B.bam
# 0.01%

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 20241121_dowload_GSM4477919-24/fastq_files/Flexbar_out3/SRR11550024_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_S1_4h_B.bam
# 0.05 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 20241121_dowload_GSM4477919-24/fastq_files/Flexbar_out3/SRR11550023_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_S1_4h_A.bam
# 0.04 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 20241121_dowload_GSM4477919-24/fastq_files/Flexbar_out3/SRR11550022_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_S1_24h_B.bam
# 3.54 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 20241121_dowload_GSM4477919-24/fastq_files/Flexbar_out3/SRR11550021_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_S1_24h_A.bam
# 3.75 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 20241121_dowload_GSM4477919-24/fastq_files/Flexbar_out3/SRR11550020_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_S1_12h_B.bam
# 1.40 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 20241121_dowload_GSM4477919-24/fastq_files/Flexbar_out3/SRR11550019_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_S1_12h_A.bam
# 0.86 overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550031_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_untr_S1_4h_A.bam
# 0.02 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550032_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_untr_S1_4h_B.bam
# 0.03 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550015_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_mock_S1_24h_A.bam
# 0.02 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550016_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_mock_S1_24h_B.bam
# 0.04 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550017_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_mock_S1_4h_A.bam
# 0.02 % overall alignment rate

bowtie2 -t -x AY310120.1_SARS_coronavirus_FRA -q 2022115_downloadGSM4477932-49/Fastq_files/Flexbar_out3/SRR11550018_flexbar3.fastq 20240121_bowtie/calu3_smallRNA_mock_S1_4h_B.bam
# 0.01 % overall alignment rate


# sort the bam alignment file as suggested https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
# v1.9 with nixpkgs/1609

# conda install bioconda::bioconductor-rsamtools # use this if installing on local server

samtools sort calu3_smallRNA_S2_12h_A.bam -o calu3_smallRNA_S2_12h_A_sort.bam
samtools sort calu3_smallRNA_S2_12h_B.bam -o calu3_smallRNA_S2_12h_B_sort.bam

samtools sort calu3_smallRNA_S2_24h_A.bam -o calu3_smallRNA_S2_24h_A_sort.bam
samtools sort calu3_smallRNA_S2_24h_B.bam -o calu3_smallRNA_S2_24h_B_sort.bam

samtools sort calu3_smallRNA_S2_4h_A.bam -o calu3_smallRNA_S2_4h_A_sort.bam
samtools sort calu3_smallRNA_S2_4h_B.bam -o calu3_smallRNA_S2_4h_B_sort.bam

samtools sort calu3_smallRNA_untr_4h_A.bam -o calu3_smallRNA_untr_4h_A_sort.bam
samtools sort calu3_smallRNA_untr_4h_B.bam -o calu3_smallRNA_untr_4h_B_sort.bam

samtools sort calu3_smallRNA_mock_4h_A.bam -o calu3_smallRNA_mock_4h_A_sort.bam
samtools sort calu3_smallRNA_mock_4h_B.bam -o calu3_smallRNA_mock_4h_B_sort.bam

samtools sort calu3_smallRNA_mock_24h_A.bam -o calu3_smallRNA_mock_24h_A_sort.bam
samtools sort calu3_smallRNA_mock_24h_B.bam -o calu3_smallRNA_mock_24h_B_sort.bam

samtools sort calu3_smallRNA_S1_4h_A.bam -o calu3_smallRNA_S1_4h_A_sort.bam 
samtools sort calu3_smallRNA_S1_4h_B.bam -o calu3_smallRNA_S1_4h_B_sort.bam

samtools sort calu3_smallRNA_S1_12h_A.bam -o calu3_smallRNA_S1_12h_A_sort.bam 
samtools sort calu3_smallRNA_S1_12h_B.bam -o calu3_smallRNA_S1_12h_B_sort.bam

samtools sort calu3_smallRNA_S1_24h_A.bam -o calu3_smallRNA_S1_24h_A_sort.bam 
samtools sort calu3_smallRNA_S1_24h_B.bam -o calu3_smallRNA_S1_24h_B_sort.bam

samtools sort calu3_smallRNA_untr_S1_4h_A.bam -o calu3_smallRNA_untr_S1_4h_A_sort.bam
samtools sort calu3_smallRNA_untr_S1_4h_B.bam -o calu3_smallRNA_untr_S1_4h_B_sort.bam

samtools sort calu3_smallRNA_mock_S1_4h_A.bam -o calu3_smallRNA_mock_S1_4h_A_sort.bam
samtools sort calu3_smallRNA_mock_S1_4h_B.bam -o calu3_smallRNA_mock_S1_4h_B_sort.bam

samtools sort calu3_smallRNA_mock_S1_24h_A.bam -o calu3_smallRNA_mock_S1_24h_A_sort.bam
samtools sort calu3_smallRNA_mock_S1_24h_B.bam -o calu3_smallRNA_mock_S1_24h_B_sort.bam

##-------------------------------------------------------------------------
##                4.Scaling to RPM and counts                             -
##-------------------------------------------------------------------------


# scaling factor - https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/#:~:text=Count%20up%20the%20total%20reads,you%20reads%20per%20million%20(RPM)
# https://bioinformatics.ccr.cancer.gov/btep/questions/what-is-the-difference-between-rpkm-fpkm-and-tpm#:~:text=When%20you%20use%20TPM%2C%20the,harder%20to%20compare%20samples%20directly.  
# Kim et al 2021 report RPM (reads per million) - to calculate in bedtools -scale, divide total reads/1,000,000 and then use that to scale
# bedtools v2.27.1
# conda install bioconda::bedtools #this is for local servers 

bedtools genomecov -ibam calu3_smallRNA_S2_12h_A_sort.bam -d -scale 0.15961824 > calu3_smallRNA_S2_12h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_S2_12h_A_sort.bam -d > calu3_smallRNA_S2_12h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S2_12h_B_sort.bam -d -scale 0.18182003 > calu3_smallRNA_S2_12h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_S2_12h_B_sort.bam -d > calu3_smallRNA_S2_12h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S2_24h_B_sort.bam -d -scale 0.09343582 > calu3_smallRNA_S2_24h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_S2_24h_B_sort.bam -d > calu3_smallRNA_S2_24h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S2_24h_A_sort.bam -d -scale 0.07963927 > calu3_smallRNA_S2_24h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_S2_24h_A_sort.bam -d > calu3_smallRNA_S2_24h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S2_4h_B_sort.bam -d -scale 0.15348709 > calu3_smallRNA_S2_4h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_S2_4h_B_sort.bam -d > calu3_smallRNA_S2_4h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S2_4h_A_sort.bam -d -scale 0.15377542 > calu3_smallRNA_S2_4h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_S2_4h_A_sort.bam -d > calu3_smallRNA_S2_4h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_24h_A_sort.bam -d -scale 0.23123414 > calu3_smallRNA_mock_24h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_24h_A_sort.bam -d > calu3_smallRNA_mock_24h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_24h_B_sort.bam -d -scale 0.18886369 > calu3_smallRNA_mock_24h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_24h_B_sort.bam -d > calu3_smallRNA_mock_24h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_4h_A_sort.bam -d -scale 0.09805383 > calu3_smallRNA_mock_4h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_4h_A_sort.bam -d > calu3_smallRNA_mock_4h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_4h_B_sort.bam -d -scale 0.18886369 > calu3_smallRNA_mock_4h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_4h_B_sort.bam -d > calu3_smallRNA_mock_4h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_untr_4h_A_sort.bam -d -scale 0.08735529 > calu3_smallRNA_untr_4h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_untr_4h_A_sort.bam -d > calu3_smallRNA_untr_4h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_untr_4h_B_sort.bam -d -scale 0.15196309 > calu3_smallRNA_untr_4h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_untr_4h_B_sort.bam -d > calu3_smallRNA_untr_4h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S1_4h_B_sort.bam -d -scale 0.10490424 > calu3_smallRNA_S1_4h_B_hist_rpm.txt 
bedtools genomecov -ibam calu3_smallRNA_S1_4h_B_sort.bam -d > calu3_smallRNA_S1_4h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S1_4h_A_sort.bam -d -scale 0.12257915 > calu3_smallRNA_S1_4h_A_hist_rpm.txt 
bedtools genomecov -ibam calu3_smallRNA_S1_4h_A_sort.bam -d > calu3_smallRNA_S1_4h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S1_12h_B_sort.bam -d -scale 0.06748350 > calu3_smallRNA_S1_12h_B_hist_rpm.txt 
bedtools genomecov -ibam calu3_smallRNA_S1_12h_B_sort.bam -d > calu3_smallRNA_S1_12h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S1_12h_A_sort.bam -d -scale 0.18528002 > calu3_smallRNA_S1_12h_A_hist_rpm.txt 
bedtools genomecov -ibam calu3_smallRNA_S1_12h_A_sort.bam -d > calu3_smallRNA_S1_12h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S1_24h_B_sort.bam -d -scale 0.17347059 > calu3_smallRNA_S1_24h_B_hist_rpm.txt 
bedtools genomecov -ibam calu3_smallRNA_S1_24h_B_sort.bam -d > calu3_smallRNA_S1_24h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_S1_24h_A_sort.bam -d -scale 0.14319561 > calu3_smallRNA_S1_24h_A_hist_rpm.txt 
bedtools genomecov -ibam calu3_smallRNA_S1_24h_A_sort.bam -d > calu3_smallRNA_S1_24h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_S1_24h_A_sort.bam -d -scale 0.23123414 > calu3_smallRNA_mock_S1_24h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_S1_24h_A_sort.bam -d > calu3_smallRNA_mock_S1_24h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_S1_24h_B_sort.bam -d -scale 0.18886369 > calu3_smallRNA_mock_S1_24h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_S1_24h_B_sort.bam -d > calu3_smallRNA_mock_S1_24h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_S1_4h_A_sort.bam -d -scale 0.09805383 > calu3_smallRNA_mock_S1_4h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_S1_4h_A_sort.bam -d > calu3_smallRNA_mock_S1_4h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_mock_S1_4h_B_sort.bam -d -scale 0.18886369 > calu3_smallRNA_mock_S1_4h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_mock_S1_4h_B_sort.bam -d > calu3_smallRNA_mock_S1_4h_B_hist.txt

bedtools genomecov -ibam calu3_smallRNA_untr_S1_4h_A_sort.bam -d -scale 0.08735529 > calu3_smallRNA_untr_S1_4h_A_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_untr_S1_4h_A_sort.bam -d > calu3_smallRNA_untr_S1_4h_A_hist.txt

bedtools genomecov -ibam calu3_smallRNA_untr_S1_4h_B_sort.bam -d -scale 0.15196309 > calu3_smallRNA_untr_S1_4h_B_hist_rpm.txt
bedtools genomecov -ibam calu3_smallRNA_untr_S1_4h_B_sort.bam -d > calu3_smallRNA_untr_S1_4h_B_hist.txt

##-------------------------------------------------------------------------
##                          5. visualization                              -
##-------------------------------------------------------------------------

# according to Wyler et al. 2021. notes that Calu-3 cells infected with SARS-CoV-2 show clear signs of cell death at 24 hpi (Figure S1I)
# which likely influences the viral yield at the latest time point, and that data from that time point, including e.g. the level of viral RNA (Figure 1A) 
# the 24 hpi timepoint should be therefore interpreted with caution.
# export the coverage files from the cluster the next sections are written for R

# internal note  for location of the data downloaded from servr ("D:/Covid_hackthon/2021_mar13_smallRNA_transcriptome/2022_nov16_genomecov_RPM_scale/all_DL/")

#library(here)

# SCV2 data

A_12h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_12h_A_hist_rpm.txt"),header=F)
B_12h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_12h_B_hist_rpm.txt"),header=F)

A_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_24h_A_hist_rpm.txt"),header=F)
B_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_24h_B_hist_rpm.txt"),header=F)

A_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_4h_A_hist_rpm.txt"),header=F)
B_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S2_4h_B_hist_rpm.txt"),header=F)

A_mock_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_4h_A_hist_rpm.txt"),header=F)
B_mock_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_4h_B_hist_rpm.txt"),header=F)

A_mock_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_24h_A_hist_rpm.txt"),header=F)
B_mock_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_24h_B_hist_rpm.txt"),header=F)

A_untr_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_untr_4h_A_hist_rpm.txt"),header=F)
B_untr_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_untr_4h_B_hist_rpm.txt"),header=F)


colnames(A_12h) = c("chr","position","count")
colnames(B_12h) = c("chr","position","count")
colnames(A_24h) = c("chr","position","count")
colnames(B_24h) = c("chr","position","count")
colnames(A_4h) = c("chr","position","count")
colnames(B_4h) = c("chr","position","count")
colnames(A_mock_4h) = c("chr","position","count")
colnames(B_mock_4h) = c("chr","position","count")
colnames(A_mock_24h) = c("chr","position","count")
colnames(B_mock_24h) = c("chr","position","count")
colnames(A_untr_4h) = c("chr","position","count")
colnames(B_untr_4h) = c("chr","position","count")

# SCV1 data

S1_A_12h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_12h_A_hist_rpm.txt"),header=F)
S1_B_12h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_12h_B_hist_rpm.txt"),header=F)

S1_A_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_24h_A_hist_rpm.txt"),header=F)
S1_B_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_24h_B_hist_rpm.txt"),header=F)

S1_A_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_4h_A_hist_rpm.txt"),header=F)
S1_B_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_S1_4h_B_hist_rpm.txt"),header=F)

S1_A_mock_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_S1_4h_A_hist_rpm.txt"),header=F)
S1_B_mock_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_S1_4h_B_hist_rpm.txt"),header=F)

S1_A_mock_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_S1_24h_A_hist_rpm.txt"),header=F)
S1_B_mock_24h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_mock_S1_24h_B_hist_rpm.txt"),header=F)

S1_A_untr_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_untr_S1_4h_A_hist_rpm.txt"),header=F)
S1_B_untr_4h = read.table(here("data/small_RNAseq_counts/calu3_smallRNA_untr_S1_4h_B_hist_rpm.txt"),header=F)

colnames(S1_A_12h) = c("chr","position","count")
colnames(S1_B_12h) = c("chr","position","count")
colnames(S1_A_24h) = c("chr","position","count")
colnames(S1_B_24h) = c("chr","position","count")
colnames(S1_A_4h) = c("chr","position","count")
colnames(S1_B_4h) = c("chr","position","count")
colnames(S1_A_mock_4h) = c("chr","position","count")
colnames(S1_B_mock_4h) = c("chr","position","count")
colnames(S1_A_mock_24h) = c("chr","position","count")
colnames(S1_B_mock_24h) = c("chr","position","count")
colnames(S1_A_untr_4h) = c("chr","position","count")
colnames(S1_B_untr_4h) = c("chr","position","count")

#library(ggplot2)

A_12h_plot=ggplot(A_12h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()

S1_A_12h_plot =ggplot(S1_A_12h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()

# want to present mean of A and B, both 12h and 24h, plus highlights in the side-by-side

#S1- create combined data table for dataviz
merge_S1_12h = merge(S1_A_12h, S1_B_12h, by = "position")
col_order = c("chr.x", "position", "count.x","count.y")
merge_S1_12h_ = merge_S1_12h[, col_order]
colnames(merge_S1_12h_)[1] = "chr" 
merge_S1_12h_$mean.count = apply(merge_S1_12h_[,3:4],1,mean)

merge_S1_24h = merge(S1_A_24h, S1_B_24h, by = "position")
col_order = c("chr.x", "position", "count.x","count.y")
merge_S1_24h_ = merge_S1_24h[, col_order]
colnames(merge_S1_24h_)[1] = "chr" 
merge_S1_24h_$mean.count = apply(merge_S1_24h_[,3:4],1,mean)


#library(dplyr)

merge_S1_24h_a =merge_S1_24h_ %>%
  mutate(Timepoint = "24 hpi")

merge_S1_12h_a =merge_S1_12h_ %>%
  mutate(Timepoint = "12 hpi")

S1_24h_12h_mean = rbind(merge_S1_24h_a, merge_S1_12h_a)


#S2 - create the data table for data viz
merge_S2_12h = merge(A_12h, B_12h, by = "position")
col_order = c("chr.x", "position", "count.x","count.y")
merge_S2_12h_ = merge_S2_12h[, col_order]
colnames(merge_S2_12h_)[1] = "chr" 
merge_S2_12h_$mean.count = apply(merge_S2_12h_[,3:4],1,mean)

merge_S2_24h = merge(A_24h, B_24h, by = "position")
col_order = c("chr.x", "position", "count.x","count.y")
merge_S2_24h_ = merge_S2_24h[, col_order]
colnames(merge_S2_24h_)[1] = "chr" 
merge_S2_24h_$mean.count = apply(merge_S2_24h_[,3:4],1,mean)

merge_S2_24h_a =merge_S2_24h_ %>%
  mutate(Timepoint = "24 hpi")

merge_S2_12h_a =merge_S2_12h_ %>%
  mutate(Timepoint = "12 hpi")

S2_24h_12h_mean = rbind(merge_S2_24h_a, merge_S2_12h_a)

# add highlights # https://stackoverflow.com/questions/32543176/highlight-areas-within-certain-x-range-in-ggplot2

S1_rects = read.delim(here("data/20240217_SCV1_highlights.txt")) 
S2_rects = read.delim(here("data/20240217_SCV2_highlights.txt")) 

#library(scales)

S1_A_12h_plot = S1_A_12h_plot + 
geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 4000, fill = Type), alpha=0.5) +
  scale_fill_manual (values = c("miRNA" = "violetred1","svRNA_TFO" = "slateblue1","HIS_svRNA" = "orange")) +
  labs(title="SARS-CoV-1 A_12h", y = "RPM", size = 10) +
  theme(text = element_text(size=10)) + 
  theme(legend.position="right") +
  theme(legend.text=element_text(size=8)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = .1))


colors_type = c("miRNA" = "violetred1","svRNA_TFO" = "slateblue1","HIS_svRNA" = "orange", "24 hpi" = "darkgrey", "12 hpi" = "blue")

S1_avg_plot = ggplot() +
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 6500, fill = Type), alpha=0.5) +
  geom_area(data = S1_24h_12h_mean, aes(x=position,y=mean.count, fill = Timepoint), position = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c(colors_type)) +
  labs(title="SARS-CoV-1", y = "RPM", size = 10) +
  ylim(0,6500) +
  theme_bw() +
  theme(text = element_text(size=10)) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = .1))


#A_12h_plot = A_12h_plot + 
geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 4000, fill = Type), alpha=0.5) +
  scale_fill_manual (values = c("miRNA" = "violetred1","svRNA_TFO" = "slateblue1","HIS_svRNA" = "orange")) +
  labs(title="SARS-CoV-2 A_12h", y = "RPM", size = 10) +
  theme(text = element_text(size=10)) + 
  theme(legend.position="right") +
  theme(legend.text=element_text(size=8)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = .1))


S2_avg_plot = ggplot() +
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 6500, fill = Type), alpha=0.5) +
  geom_area(data = S2_24h_12h_mean, aes(x=position,y=mean.count, fill = Timepoint), position = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c(colors_type)) +
  labs(title="SARS-CoV-2", y = "RPM", size = 10) +
  ylim(0,6500) +
  theme_bw() +
  theme(text = element_text(size=10)) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = .1))


#library(scales)


#B_12h plot
B_12h_plot =ggplot(B_12h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()


S1_B_12h_plot =ggplot(S1_B_12h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()


##24h mock A and B
#A_mock_24h plot
A_mock_24h_plot =ggplot(A_mock_24h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()

S1_A_mock_24h_plot =ggplot(S1_A_mock_24h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()

S1_A_mock_24h_plot=S1_A_mock_24h_plot + geom_rect(data=S1_rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,
                                                                                     ymax=10000, group=group), color="transparent", fill="orange", alpha=0.3)+
  labs(title="SARS-CoV-1 mock_A_24h", y = "RPM")

#B_mock_24h plot
B_mock_24h_plot =ggplot(B_mock_24h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()

S1_B_mock_24h_plot =ggplot(S1_B_mock_24h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()

S1_B_mock_24h_plot=S1_B_mock_24h_plot + geom_rect(data=S1_rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,
                                                                                     ymax=10000, group=group), color="transparent", fill="orange", alpha=0.3)+
  labs(title="SARS-CoV-1 mock_B_24h", y = "RPM")

###4h A and B
#A_4h plot
A_4h_plot =ggplot(A_4h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()

S1_A_4h_plot =ggplot(S1_A_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()

A_4h_plot + ylim(0,100)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 100, fill = Type), alpha=0.5)+
  scale_fill_manual (values = c("miRNA" = "violetred1","svRNA_TFO" = "slateblue1","HIS_svRNA" = "orange"))+
  labs(title="SARS-CoV-2 A_4h", y = "RPM")

S1_A_4h_plot + ylim(0,100)+
  geom_rect(data = S1_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 100, fill = Type), alpha=0.5)+
  scale_fill_manual (values = c("miRNA" = "violetred1","svRNA_TFO" = "slateblue1","HIS_svRNA" = "orange"))+
  labs(title="SARS-CoV-1 A_4h", y = "RPM")


#B_4h plot
B_4h_plot =ggplot(B_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()

S1_B_4h_plot =ggplot(S1_B_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()



##-------------------------------------------------------------------------
##                     Multiple plot function                             -
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



# save the multiplot

pdf("20240219_12_24_mean_multi.pdf", height=2.9173227, width=7.20472)
multiplot(S2_avg_plot, S1_avg_plot, cols=1)
dev.off()   

###4h mock A and B
#A_mock_4h plot
A_mock_4h_plot =ggplot(A_mock_4h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()


#B_mock_4h plot
B_mock_4h_plot =ggplot(B_mock_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()


S1_B_mock_4h_plot =ggplot(S1_B_mock_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()

S1_B_mock_4h_plot=S1_B_mock_4h_plot + geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,
                                                                                   ymax=10000, group=group), color="transparent", fill="orange", alpha=0.3)+
  labs(title="SARS-CoV-1 mock B_4h", y = "RPM")


###4h untr A and B
#A_untr_4h plot
A_untr_4h_plot =ggplot(A_untr_4h,aes(x=position,y=count))+
  geom_area(fill="blue")+
  theme_bw()


S1_A_untr_4h_plot =ggplot(S1_A_untr_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()

S1_A_untr_4h_plot=S1_A_untr_4h_plot + geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,
                                                                                   ymax=10000, group=group), color="transparent", fill="orange", alpha=0.3)+
  labs(title="SARS-CoV-1 untr_A_4h", y = "RPM")


#B_mock_4h plot
B_untr_4h_plot =ggplot(B_untr_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()


S1_B_untr_4h_plot =ggplot(S1_B_untr_4h,aes(x=position,y=count))+
  geom_area(fill="darkslategrey")+
  theme_bw()

S1_B_untr_4h_plot=S1_B_untr_4h_plot + geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,
                                                                                   ymax=10000, group=group), color="transparent", fill="orange", alpha=0.3)+
  labs(title="SARS-CoV-1 untr_B_4h", y = "RPM")


# load the multiplot function (above) before rendering the plot and saving it as a pdf file.

pdf("20240121_4h_multi.pdf", width=7.5, height=10)
multiplot(S1_A_4h_plot,S1_A_mock_4h_plot,S1_A_untr_4h_plot, S1_B_4h_plot, S1_B_mock_4h_plot,S1_B_untr_4h_plot, cols=2)
dev.off()


##-------------------------------------------------------------------------
##                      6. Zoomed in A vs B coverage                      -
##-------------------------------------------------------------------------

# combine A and B for zoom plotting

#library(dplyr)

# SCV2
B_24ha =B_24h %>%
  mutate(replicate = "B")

A_24ha =A_24h %>%
  mutate(replicate = "A")

AB_24h = rbind(A_24ha, B_24ha)


A_12ha =A_12h %>%
  mutate(replicate = "A")

B_12ha =B_12h %>%
  mutate(replicate = "B")

AB_12h = rbind(A_12ha, B_12ha)

A_4ha =A_4h %>%
  mutate(replicate = "A")

B_4ha =B_4h %>%
  mutate(replicate = "B")

AB_4h = rbind(A_4ha, B_4ha)

A_mock_4ha =A_mock_4h %>%
  mutate(replicate = "A")

B_mock_4ha =B_mock_4h %>%
  mutate(replicate = "B")

AB_mock_4h = rbind(A_mock_4ha, B_mock_4ha)


A_mock_24ha =A_mock_24h %>%
  mutate(replicate = "A")

B_mock_24ha =B_mock_24h %>%
  mutate(replicate = "B")

AB_mock_24h = rbind(A_mock_24ha, B_mock_24ha)

A_untr_4ha = A_untr_4h %>%
  mutate(replicate = "A")

B_untr_4ha = B_untr_4h %>%
  mutate(replicate = "B")

AB_untr_4h = rbind(A_untr_4ha, B_untr_4ha)

# SCV1
S1_B_24ha =S1_B_24h %>%
  mutate(replicate = "B")

S1_A_24ha =S1_A_24h %>%
  mutate(replicate = "A")

S1_AB_24h = rbind(S1_A_24ha, S1_B_24ha)


S1_A_12ha =S1_A_12h %>%
  mutate(replicate = "A")

S1_B_12ha =S1_B_12h %>%
  mutate(replicate = "B")

S1_AB_12h = rbind(S1_A_12ha, S1_B_12ha)

S1_A_4ha =S1_A_4h %>%
  mutate(replicate = "A")

S1_B_4ha =S1_B_4h %>%
  mutate(replicate = "B")

S1_AB_4h = rbind(S1_A_4ha, S1_B_4ha)

S1_A_mock_4ha =S1_A_mock_4h %>%
  mutate(replicate = "A")

S1_B_mock_4ha =S1_B_mock_4h %>%
  mutate(replicate = "B")

S1_AB_mock_4h = rbind(S1_A_mock_4ha, S1_B_mock_4ha)


S1_A_mock_24ha =S1_A_mock_24h %>%
  mutate(replicate = "A")

S1_B_mock_24ha =S1_B_mock_24h %>%
  mutate(replicate = "B")

S1_AB_mock_24h = rbind(S1_A_mock_24ha, S1_B_mock_24ha)


#zoom in to specific area

# N-ORF and 3'UTR
S1_N3UTR_24h = subset(S1_AB_24h, position>28119 & position<29741)
S1_N3UTR_24h_noAAA = subset(S1_AB_24h, position>28119 & position<29728)

#highlight the miRNA region from Morales et al. 2017 located from 28461 - 28481
sCV1_N_start <- c(28461,29717)
SCV1_N_end <- c(28481,29740)
morales_li_rects <- data.frame(start=SCV1_N_end, end=SCV1_N_end, group=c(1,2))

S1_N3UTR_24h_plot = ggplot(S1_N3UTR_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data=morales_li_rects, inherit.aes=FALSE, aes(xmin=sCV1_N_start, xmax=SCV1_N_end, ymin=0, ymax=10000, group=group), color="transparent", fill="orange", alpha=0.3)+
  labs(title="SARS-CoV-1 24h N-ORF_3'UTR", y = "RPM")+
  theme_bw()

pdf("2024_N_3UTR_S1_Moralies_Li.pdf", width=7.5, height=6)
S1_N3UTR_24h_plot
dev.off()

S1_N3UTR_24h_noAAA_plot = ggplot(S1_N3UTR_24h_noAAA,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  theme_bw()

# set highlight colors - distinguishing between svRNA-TFOs (pink) and TFOs (grey) etc.
colors_type_new = c("miRNA" = "slateblue1","svRNA_TFO" = "deeppink1", "TFO" = "darkviolet", "HIS_svRNA" = "orange")

# regions to scan that have TFOs, the 120 nt region is based on either a predicted pre-miRNA structure or Rangan CS region (see Materials and MEthods of Manuscript), or else the region was arbitrarily centred on the the TFO.


#TFO_nsp2: 1076 - 1195 
#svRNA_nsp2: 1558 - 1677
#TFO_nsp2: 2018 - 2137
#svRNA_nsp3: 2997 - 3116
#TFO_nsp3: 5867 - 5986
#TFO_nsp4: 9141 - 9260
#TFO_nsp6: 11104 - 11223
#TFO_nsp6: 11730 - 11849
#TFO_nsp12: 15645 - 15764
#TFO_nsp15: 19862 - 19981
#TFO_S: 21701 - 21820
#TFO_S: 21926 - 22045
#svRNA_S.1: 23304 - 23423
#svRNA_S.2: 23924 - 24043
#TFO_ORF3a: 25618 - 25737 (both TFOs are captured) #####
#svRNA_E: 26255 - 26374
#svRNA_N: 29343 - 29492


#nsp2.1

nsp2.1_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1076,1195)+
  labs(title="nsp2.1_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp2.1_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1076,1195)+
  labs(title="nsp2.2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.1_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1076,1195)+
  labs(title="nsp2.1_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.1_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1076,1195)+
  labs(title="nsp2.1_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.1_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1076,1195)+
  labs(title="nsp2.1_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.1_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1076,1195)+
  labs(title="nsp2.1_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(nsp2.1_24h_plot_notsub, nsp2.1_12h_plot_notsub, nsp2.1_4h_plot_notsub,nsp2.1_mock_24h_plot_notsub, nsp2.1_mock_4h_plot_notsub,nsp2.1_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp2.2

nsp2.2_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1558,1677)+
  labs(title="svRNA_TFO_nsp2_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))


nsp2.2_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1558,1677)+
  labs(title="svRNA_TFO_nsp2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.2_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1558,1677)+
  labs(title="svRNA_TFO_nsp2_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.2_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1558,1677)+
  labs(title="svRNA_TFO_nsp2_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.2_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1558,1677)+
  labs(title="svRNA_TFO_nsp2_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.2_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(1558,1677)+
  labs(title="svRNA_TFO_nsp2_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(nsp2.2_24h_plot_notsub, nsp2.2_12h_plot_notsub, nsp2.2_4h_plot_notsub,nsp2.2_mock_24h_plot_notsub, nsp2.2_mock_4h_plot_notsub,nsp2.2_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp2.3
nsp2.3_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2018,2137)+
  labs(title="nsp2.2_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))


nsp2.3_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2018,2137)+
  labs(title="nsp2.2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.3_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2018,2137)+
  labs(title="nsp2.2_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.3_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2018,2137)+
  labs(title="nsp2.2_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.3_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2018,2137)+
  labs(title="nsp2.2_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp2.3_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2018,2137)+
  labs(title="nsp2.2_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(nsp2.3_24h_plot_notsub, nsp2.3_12h_plot_notsub, nsp2.3_4h_plot_notsub,nsp2.3_mock_24h_plot_notsub, nsp2.3_mock_4h_plot_notsub,nsp2.3_untr_4h_plot_notsub, cols=2)
#saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement


# svRNA_TFO_nsp3
nsp3.1_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2997,3116)+
  labs(title="svRNA_TFO_nsp3_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp3.1_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2997,3116)+
  labs(title="svRNA_TFO_nsp3_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.1_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2997,3116)+
  labs(title="svRNA_TFO_nsp3_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.1_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2997,3116)+
  labs(title="svRNA_TFO_nsp3_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.1_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2997,3116)+
  labs(title="svRNA_TFO_nsp3_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.1_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(2997,3116)+
  labs(title="svRNA_TFO_nsp3_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(nsp3.1_24h_plot_notsub, nsp3.1_12h_plot_notsub, nsp3.1_4h_plot_notsub,nsp3.1_mock_24h_plot_notsub, nsp3.1_mock_4h_plot_notsub,nsp3.1_untr_4h_plot_notsub, cols=2)
#saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp3.2
nsp3.2_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(5867,5986)+
  labs(title="nsp3_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp3.2_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(5867,5986)+
  labs(title="nsp3_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.2_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(5867,5986)+
  labs(title="nsp3_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.2_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(5867,5986)+
  labs(title="nsp3_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.2_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(5867,5986)+
  labs(title="nsp3_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp3.2_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(5867,5986)+
  labs(title="nsp3_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(nsp3.2_24h_plot_notsub, nsp3.2_12h_plot_notsub, nsp3.2_4h_plot_notsub,nsp3.2_mock_24h_plot_notsub, nsp3.2_mock_4h_plot_notsub,nsp3.2_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp4
nsp4_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(9141,9260)+
  labs(title="nsp4_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp4_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(9141,9260)+
  labs(title="nsp4_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp4_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(9141,9260)+
  labs(title="nsp4_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp4_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(9141,9260)+
  labs(title="nsp4_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp4_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(9141,9260)+
  labs(title="nsp4_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp4_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(9141,9260)+
  labs(title="nsp4_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(nsp4_24h_plot_notsub, nsp4_12h_plot_notsub, nsp4_4h_plot_notsub,nsp4_mock_24h_plot_notsub, nsp4_mock_4h_plot_notsub,nsp4_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp6.1
nsp6.1_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11104,11223)+
  labs(title="nsp6.1_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp6.1_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11104,11223)+
  labs(title="nsp6.1_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.1_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11104,11223)+
  labs(title="nsp6.1_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.1_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11104,11223)+
  labs(title="nsp6.1_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.1_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11104,11223)+
  labs(title="nsp6.1_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.1_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11104,11223)+
  labs(title="nsp6.1_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(nsp6.1_24h_plot_notsub, nsp6.1_12h_plot_notsub, nsp6.1_4h_plot_notsub,nsp6.1_mock_24h_plot_notsub, nsp6.1_mock_4h_plot_notsub,nsp6.1_untr_4h_plot_notsub, cols=2)
#saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp6.2
nsp6.2_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11730,11849)+
  labs(title="nsp6.2_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp6.2_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11730,11849)+
  labs(title="nsp6.2_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.2_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11730,11849)+
  labs(title="nsp6.2_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.2_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11730,11849)+
  labs(title="nsp6.2_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.2_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11730,11849)+
  labs(title="nsp6.2_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp6.2_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(11730,11849)+
  labs(title="nsp6.2_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(nsp6.2_24h_plot_notsub, nsp6.2_12h_plot_notsub, nsp6.2_4h_plot_notsub,nsp6.2_mock_24h_plot_notsub, nsp6.2_mock_4h_plot_notsub,nsp6.2_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement


# nsp12.1
nsp12.1_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(15645,15764)+
  labs(title="nsp12.1_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp12.1_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(15645,15764)+
  labs(title="nsp12.1_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp12.1_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(15645,15764)+
  labs(title="nsp12.1_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp12.1_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(15645,15764)+
  labs(title="nsp12.1_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp12.1_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(15645,15764)+
  labs(title="nsp12.1_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp12.1_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(15645,15764)+
  labs(title="nsp12.1_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(nsp12.1_24h_plot_notsub, nsp12.1_12h_plot_notsub, nsp12.1_4h_plot_notsub,nsp12.1_mock_24h_plot_notsub, nsp12.1_mock_4h_plot_notsub,nsp12.1_untr_4h_plot_notsub, cols=2)
#saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# nsp15
nsp15_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(19862,19981)+
  labs(title="nsp15_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

nsp15_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(19862,19981)+
  labs(title="nsp15_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp15_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(19862,19981)+
  labs(title="nsp15_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp15_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(19862,19981)+
  labs(title="nsp15_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp15_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(19862,19981)+
  labs(title="nsp15_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

nsp15_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(19862,19981)+
  labs(title="nsp15_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(nsp15_24h_plot_notsub, nsp15_12h_plot_notsub, nsp15_4h_plot_notsub,nsp15_mock_24h_plot_notsub, nsp15_mock_4h_plot_notsub,nsp15_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement


# S.1 TFO
S.1_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21701,21820)+
  labs(title="S.1_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

S.1_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21701,21820)+
  labs(title="S.1_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.1_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21701,21820)+
  labs(title="S.1_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.1_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21701,21820)+
  labs(title="S.1_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.1_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21701,21820)+
  labs(title="S.1_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.1_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21701,21820)+
  labs(title="S.1_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(S.1_24h_plot_notsub, S.1_12h_plot_notsub, S.1_4h_plot_notsub,S.1_mock_24h_plot_notsub, S.1_mock_4h_plot_notsub,S.1_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# S.2 TFO
S.2_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21926,22045)+
  labs(title="S.2_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

S.2_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21926,22045)+
  labs(title="S.2_TFO_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.2_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21926,22045)+
  labs(title="S.2_TFO_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.2_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21926,22045)+
  labs(title="S.2_TFO_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.2_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21926,22045)+
  labs(title="S.2_TFO_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.2_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(21926,22045)+
  labs(title="S.2_TFO_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(S.2_24h_plot_notsub, S.2_12h_plot_notsub, S.2_4h_plot_notsub,S.2_mock_24h_plot_notsub, S.2_mock_4h_plot_notsub,S.2_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement


# svRNA_TFO_S.1
S.3_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23304,23423)+
  labs(title="svRNA_TFO_S.1_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

S.3_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23304,23423)+
  labs(title="svRNA_TFO_S.1_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.3_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23304,23423)+
  labs(title="svRNA_TFO_S.1_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.3_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23304,23423)+
  labs(title="svRNA_TFO_S.1_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.3_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23304,23423)+
  labs(title="svRNA_TFO_S.1_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.3_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23304,23423)+
  labs(title="svRNA_TFO_S.1_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(S.3_24h_plot_notsub, S.3_12h_plot_notsub, S.3_4h_plot_notsub,S.3_mock_24h_plot_notsub, S.3_mock_4h_plot_notsub,S.3_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# svRNA_TFO_S.2
S.4_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23924,24043)+
  labs(title="svRNA_TFO_S.2_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

S.4_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23924,24043)+
  labs(title="svRNA_TFO_S.2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.4_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23924,24043)+
  labs(title="svRNA_TFO_S.2_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.4_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23924,24043)+
  labs(title="svRNA_TFO_S.2_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.4_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23924,24043)+
  labs(title="svRNA_TFO_S.2_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

S.4_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(23924,24043)+
  labs(title="svRNA_TFO_S.2_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(S.4_24h_plot_notsub, S.4_12h_plot_notsub, S.4_4h_plot_notsub,S.4_mock_24h_plot_notsub, S.4_mock_4h_plot_notsub,S.4_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement


# ORF3a.1 and ORF3a.2 (both captured in the zoomed plot)
ORF3a_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(25618,25737)+
  labs(title="ORF3a_TFOs_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

ORF3a_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(25618,25737)+
  labs(title="ORF3a_TFOs_S.2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

ORF3a_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(25618,25737)+
  labs(title="ORF3a_TFOs_S.2_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

ORF3a_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(25618,25737)+
  labs(title="ORF3a_TFOs_S.2_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

ORF3a_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(25618,25737)+
  labs(title="ORF3a_TFOs_S.2_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

ORF3a_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(25618,25737)+
  labs(title="ORF3a_TFOs_S.2_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

multiplot(ORF3a_24h_plot_notsub, ORF3a_12h_plot_notsub, ORF3a_4h_plot_notsub,ORF3a_mock_24h_plot_notsub, ORF3a_mock_4h_plot_notsub,ORF3a_untr_4h_plot_notsub, cols=2)
# saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# svRNA_TFO_E
E_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(26255,26374)+
  labs(title="E_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

E_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(26255,26374)+
  labs(title="E_TFO_S.2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

E_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(26255,26374)+
  labs(title="E_TFO_S.2_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

E_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(26255,26374)+
  labs(title="E_TFO_S.2_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

E_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(26255,26374)+
  labs(title="E_TFO_S.2_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

E_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 500, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(26255,26374)+
  labs(title="E_TFO_S.2_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(E_24h_plot_notsub, E_12h_plot_notsub, E_4h_plot_notsub,E_mock_24h_plot_notsub, E_mock_4h_plot_notsub,E_untr_4h_plot_notsub, cols=2)
#saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# N

N_24h_plot_notsub = ggplot(AB_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(29343,29492)+
  labs(title="N_TFO_24h", y = "RPM")+
  theme_bw()+
  theme(legend.key.size = unit(0.5,'cm'))+
  theme(legend.title = element_text(size=10))+
  theme(legend.text = element_text(size=10))+
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

N_12h_plot_notsub = ggplot(AB_12h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(29343,29492)+
  labs(title="N_TFO_S.2_12h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

N_4h_plot_notsub = ggplot(AB_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(29343,29492)+
  labs(title="N_TFO_S.2_4h", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

N_mock_4h_plot_notsub = ggplot(AB_mock_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(29343,29492)+
  labs(title="N_TFO_S.2_4h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

N_mock_24h_plot_notsub = ggplot(AB_mock_24h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(29343,29492)+
  labs(title="N_TFO_S.2_24h_mock", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))

N_untr_4h_plot_notsub = ggplot(AB_untr_4h,aes(x=position,y=count))+
  geom_area(aes(colour = replicate, fill= replicate), position = 'identity', alpha=0.1)+
  geom_rect(data = S2_rects, inherit.aes=FALSE, aes(xmin = start,xmax = end,  ymin = 0, ymax = 5000, fill = Type), alpha=0.5)+
  scale_fill_manual(values = c(colors_type_new)) +
  xlim(29343,29492)+
  labs(title="N_TFO_S.2_4h_untreated", y = "RPM")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = 0.1, r = 3.5, b = 0.1, l = 0.1, unit = "cm"))


multiplot(N_24h_plot_notsub, N_12h_plot_notsub, N_4h_plot_notsub,N_mock_24h_plot_notsub, N_mock_4h_plot_notsub,N_untr_4h_plot_notsub, cols=2)
#saved as PDF 8.5 x 8.5 inches is viewable for including in the supplement

# write out the AB_24h table to examine RPMs

write.table(AB_24h, file=here("data/AB24h.txt"),sep = ",", row.names=FALSE, col.names = TRUE)


##-------------------------------------------------------------------------
##    7. Optional codes for calculating coverage over a window            -
##-------------------------------------------------------------------------

# https://www.biostars.org/p/392698/ 

# https://www.biostars.org/p/70577/

# the Kim et al (2021) paper partitioned the genome into 500 bins of 60 nt length, maximum coverage for each, then some were plotted at single nucleotide resolution

### https://www.biostars.org/p/9522228/  some instructions how to do the binning with bedtools

bedtools makewindows -n 500 -b chrom.sizes2.bed > genome.500windows.bed

bedtools map -a genome.500windows.bed -b bam2bed_sort.bed -o sum > out_map_sum

data2 = read.table("out_map_sum",sep="\t",header=F)

colnames(data2) = c("chr","start","end","count")
meanx=mean(data2$count)

head(data2)
sum(data2$count)
# 5,062,344

# End of script

#sessionInfo()
#R version 4.2.1 (2022-06-23 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 26100)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8
#[4] LC_NUMERIC=C                    LC_TIME=English_Canada.utf8    

#attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scales_1.2.1  dplyr_1.1.2   ggplot2_3.4.4 here_1.0.1    MASS_7.3-57  
