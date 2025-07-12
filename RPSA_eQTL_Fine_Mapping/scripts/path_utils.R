#path_utils.R
#Centralizes common paths for RPSA eQTL fine mapping 

library(here)

# Define the base workflow directory (relative to project root)
RPSA_eQTL_Fine_Mapping_dir <- here("RPSA_eQTL_Fine_Mapping")

# Define subdirectories to use
data_dir    <- file.path(RPSA_eQTL_Fine_Mapping_dir, "data")
scripts_dir <- file.path(RPSA_eQTL_Fine_Mapping_dir, "scripts")
results_dir <- file.path(RPSA_eQTL_Fine_Mapping_dir, "Results")
pip_plot_dir <- file.path(results_dir, "SuSIE_PIP_Plots")
Coloc_dir <- file.path(results_dir, "SuSIE_Coloc_Results")

