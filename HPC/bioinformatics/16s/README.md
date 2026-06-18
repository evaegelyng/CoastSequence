Package versions: cutadapt 4.1, dada2 1.22.0, ggplot2 3.3.6, phyloseq 1.38.0, plyr 1.8.7, RColorBrewer 1.1_3, reshape2 1.4.4, scales 1.2.1, sickle-trim 1.33, stringr 1.4.1, vegan 2.6_2

Overview of folders/files and their contents: The folder "Basic processing" contains the gwf workflow file (workflow.py) and a folder with scripts for the basic processing of raw amplicon sequencing data, from demultiplexing of fastq files to the generation of an ASV table with a corresponding file of ASV sequences. The folder "Final_scripts" contains scripts for further processing, including taxonomic assignment, metadata formatting, ASV filtering based on negative controls and prevalence, and normalization by rarefying.

Instructions: The scripts for basic processing are all run from the workflow.py file. For the further processing, the scripts were run in the following order:

1. workflow.py (calling the script dada2_tax.r)
2. scripts/make_metadata_EES.R
3. scripts/clean_up_ASV_wise_EES.R (through scripts/clean_up_ASV.sh)
4. scripts/no_sing_ASV_wise_EES.R
5. scripts/clean_tax_ASV_silva_EES.R
6. scripts/normalize.r