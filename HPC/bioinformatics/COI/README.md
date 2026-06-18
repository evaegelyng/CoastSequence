Package versions: 

- ape 5.6_2
- biostrings 2.66.0
- blast 2.12.0
- data.table 1.14.4
- decipher 2.26.0
- dnoise 1.1
- dplyr 1.0.10
- ggplot2 3.3.6
- lattice 0.20_45
- lulu 0.1.0
- mjolnir 0.1.0
- obitools 1.2.13
- optparse 1.7.3
- phyloseq 1.38.0
- plyr 1.8.7
- python 3.8.13
- RColorBrewer 1.1_3
- reshape2 1.4.4
- scales 1.2.1
- seqinr 4.2_16
- seqkit 2.1.0
- stringr 1.4.0
- swarm 3.1.0
- taxizedb 0.3.0
- tibble 3.1.8
- tidyr 1.1.4
- vegan 2.6_2
- withr 2.4.2

Overview of folders/files and their contents: There are two gwf workflow files (workflow_part1.py and workflow.py) and a folder containing scripts that are either called from the workflow files or run separately. The folder conda_envs contains yml files describing the conda environments used to run the scripts.

Instructions: The scripts should be run in the following order:

1. fastqc.sh # Run FastQC on raw reads
2. workflow_part1.py # The first part of the workflow that includes demultiplexing, merging reads, quality filtering, dereplication, chimera removal, clustering into OTUs and taxonomic assignment
3. barcode_gaps.R # Estimate barcode gap to determine an appropriate minimum similarity threshold for BLAST, and if needed, rerun workflow_part1.py 
4. make_metadata.r # Produce complete metadata file
5. clean_up_OTU_wise.r # Clean the dataset based on blank controls
6. combine_cont_lists.sh and combine_cont_lists_detailed.sh # Combine contaminant lists across negative controls for an overview
7. no_sing_OTU_wise.r # Remove sequences found in a single PCR replicate
8. merge_otu_table_w_classified.r # The cleaned OTU table is then merged with the taxonomy from the output file of MJOLNIR (COSQ_final_dataset)
9. normalize_250303.r # Rarefy reads to normalize sequencing depth across samples

The complete dataset was then downloaded and the remaining analysis was run locally (see "HMSC_and_local" repository). The subset of sequences with hits of 97% similarity or more, was manually curated and then normalized

To get ASV level results, OTUs that could be confidently identified to marine species, and were found in at least 2 clusters were selected using select_OTUs_250505.R. Then, continue with workflow.py, followed by:

1. asv.sh # The denoised MOTU files were combined  
2. Script2.1_DnoisE_230310.r # Further filter OTUs and remove NUMTs
3. no_sing_selected.R # Remove sequences found in a single PCR replicate
4. normalize_ASV.r # Rarefy reads to normalize sequencing depth across samples