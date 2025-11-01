<p align="center">
<img src="https://github.com/damselflywingz/EpiVIRTEX/blob/master/EpiVIRTEX_logo.svg" alt="EpiVIRTEX logo" width="200"/>
</p>

# Epigenetic Viral Interference through RNA Triplex Exploration

## *Introduction*
This README provides scripts and other resources used for data analysis and visualization from: <i>Paulson A.R., Montoya V., and Joy, J.B. 2025. Functional  genomic analysis reveals mechanisms of epigenetic interference in SARS-CoV-1 and SARS-CoV-2. (submitted)</i>. 

For more information about the EpiVIRTEX pipeline please see the materials and methods of the manuscript.

## Built With 

* [Decipher](http://bioconductor.org/packages/release/bioc/html/DECIPHER.html) was used to create multiple sequence alignments
* [Synplot2](https://github.com/AndrewFirth12/synplot2) was used for synonymous-site conservation analysis
* [short-BLASTn](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/) was used for host-homology searches
* [IQ-TREE 2](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/) was used for phylogenetic reconstructions
* [ggplot2](https://ggplot2.tidyverse.org/) was used to visualize the results of our functional genomics analysis

## Getting Started

To get started:

1. Retrieve the entire https://github.com/damselflywingz/EpiVIRTEX directory, and move it and its contents to your desired location;
2. Open the Rproject "Paulson_AR_SARS_ncRNA-analysis.Rproj" in RStudio;</li>
3. Follow the scripts provided (see: https://github.com/damselflywingz/EpiVIRTEX/tree/main/Scripts)

## Small RNA-sequencing analysis

The analysis involves:

* Calculating counts per million (CPM) coverage of viral small RNA-sequencing count data downloaded from the [GSE148729](https://www.ncbi.nlm.nih.gov/search/all/?term=GSE148729);
* To avoid processing of raw reads the count data from bedtools genomecov for each of the [Wyler et al. 2021](https://www.cell.com/iscience/fulltext/S2589-0042(21)00119-X?elqTrackId=a5fb8ddc3ea24283adca366859fed54e) small RNA-seq libraries have also been provided in the EpiVIRTEX/data/ folder.
* See the manuscript and refer to [Updated20240121__2022_15_15_script_smallRNA_calu_3_proc_align_quant_S1_and_S2.R](https://github.com/damselflywingz/EpiVIRTEX/tree/main/Scripts/Updated20240121__2022_15_15_script_smallRNA_calu_3_proc_align_quant_S1_and_S2.R) for more information.
 
## Recombination breakpoint and hotspot analysis

This analysis involves:

* Using the R Stats Package to run permutation for testing non-random association between TFOs of concern and known genomic recombination breakpoints and hotspots reported in [Lytras et al. 2022](https://academic.oup.com/gbe/article/14/2/evac018/6524630);
* See the manuscript and refer to [EpiVIRTEX/Scripts/20240308_random_breakpoint_script.R](https://github.com/damselflywingz/EpiVIRTEX/tree/main/Scripts/20240308_random_breakpoint_script.R) for more information.

## Host-homology searches and visualization

This analysis involves:

* Running short-BLASTn searches for matches between the genomes of SARS-CoV-2 (Wuhan-Hu-1) or SARS-CoV-1 (Tor2) and the reference genomes of several host animals.
* See the manuscript and refer to [EpiVIRTEX/Scripts/20240308_R_script_blast_results_heatmap_PCA.R](https://github.com/damselflywingz/EpiVIRTEX/tree/main/Scripts/20240308_R_script_blast_results_heatmap_PCA.R) for more information.

## Synonymous-site conservation analysis and combined genomic visualization

This analysis involves:

* Multiple sequence alignment, phylogenetic reconstruction and visualization, tip label extraction, and synonymous-site conservation analysis.
* See the manuscript and refer to [EpiVIRTEX/Scripts/20240905_SynPlot2_Align_tree.R](https://github.com/damselflywingz/EpiVIRTEX/tree/main/Scripts/20240905_SynPlot2_Align_tree.R) for more information.

## Venn diagrams for potential TFO targets

This analysis involves:

* Visualization of gene-target lists in Venn Diagrams combining three different SARS-CoV-2 transcriptome studies and enrichment analysis.
* See the manuscript and refer to [EpiVIRTEX/Scripts/20251101_DEG_Venn_Enrichment_PaulsonAR_2025.R](https://github.com/damselflywingz/EpiVIRTEX/tree/main/Scripts/20251101_DEG_Venn_Enrichment_PaulsonAR_2025.R) for more information.

## Contact

Amber Paulson - <a href="https://twitter.com/dragonflywingz" target="_blank">@Dragonflywingz </a> - Amber [dot] Rose [dot] Paulson [at] gmail [dot] com

Project Link: [https://github.com/damselflywingz/EpiVIRTEX/](https://github.com/damselflywingz/EpiVIRTEX)

## Licensing and copyright

Licensed under Creative Commons Attribution 4.0 Internatinal (CC-BY-NC 4.0)
[https://creativecommons.org/licenses/by-nc/4.0/](https://creativecommons.org/licenses/by-nc/4.0/)
