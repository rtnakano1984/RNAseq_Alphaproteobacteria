#!/bin/bash

# bash script to run AME and DREME using the MEME suite for given gene clusters
# Originally by Ryohei Thomas Nakano; nakano@mpipz.mpg.de
# 17 Feb 2021

# usage:
# 1) Copy config file to your local data directory and edit it as necessary
# 2) Log in to an HPC cluster node
# 3) Activate LSF system by lsfenv
# 4) Initiate the pipeline by:
#    /biodata/dep_psl/grp_psl/ThomasN/scripts/MEME_custom/ame.sh /path/to/your/config/ame_config.sh
# 5) Once all done, process the result files by:
#    /biodata/dep_psl/grp_psl/ThomasN/scripts/MEME_custom/ame_post.sh /path/to/your/config/ame_config.sh



# path to the cluster file (two columns delimiated by deliminated by tabs, with header "Gene\tCluster" in the first line)
export cluster_file="/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/k_means-zero_centered-9.sorted_hclust.txt"

# path to the directory where all intermediate files will be saved (no need to be already present)
export out_dir="/netscratch/dep_psl/grp_psl/ThomasN/MEME_suite/project3598/k_means-zero_centered-9"

# path to the directory you want to store the final results
export results="/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/AME-k_means-zero_centered-9"



# =============== Usually common for all analyses, no need to change =============== #
# path to the directory with header files
export header="/biodata/dep_psl/grp_psl/ThomasN/MEME_suite"

# path to the directory you store the ame_ind.sh script
export script_dir="/biodata/dep_psl/grp_psl/ThomasN/scripts"

# path to Rscript
export R_PATH="/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin"

# path to the motif database
export motif_dir="/netscratch/dep_psl/grp_psl/ThomasN/resources/motif_databases/ARABD"

# path to the genome data
export source="/netscratch/dep_psl/grp_psl/ThomasN/resources"

