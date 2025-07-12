This directory contains scripts, data, and supporting files to reproduce the fine mapping and allele fold change calculations in this manuscript

Some required input files (VCFs and gene-dosage files) contain private genetic data and are not uploaded to GitHub. See below for details on requesting access.

---
Directory Structure

data/
Contains raw input files and intermediate data:
- RPSA_Expression/: Fully processed, filtered and normalized gene expression matrices from GTEx v8, subset to just RPSA expression in the listed tissues (https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar), left empty due to GitHub size restrictions
- RPSA_Covariate_Adjusted_Expression/: Intermediate covariate-adjusted expression files produced by script 03_Covariate_Regression.py
- GTEx_VCF_files/: Private genotype VCF data from GTEx v8 (not uploaded). See notes below
- GTEx_Raw_Tissue_RNA_Counts/: Raw RNA counts from GTEx v8 (https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)
- GTEx_PLINK_files_for_RPSA/: Private gene-dosage files obtained from GTEx v8 vcf files (not uploaded). See notes below
- GTEx_Normalized_Tissue_RNA_counts/: Intermediate normalized RNA count files produced by script 02_Normalize_Log_Transform_Counts.R, left empty due to GitHub size restrictions
- GTEx_Analysis_v8_eQTL_covariates/: Covariate files from GTEx v8 (https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_covariates.tar.gz), left empty due to GitHub size restrictions
-Tissue_COV_PLINK: helper text file for 01_Fine_Map_Coloc_eQTLs.R


scripts/
Scripts to run the workflow. All scripts are written in Python or R and numbered in execution order:
- 01_Fine_Map_Coloc_eQTLs.R: Fine-mapping using SuSIE, pairwise colocalization, and creation of union credible sets  
- 02_Normalize_Log_Transform_Counts.R: DESeq normalization and log transformation of raw RNAseq counts  
- 03_Covariate_Regression.py**: Covariate regression on normalized counts  
- path_utils.R: Defines project root and data paths for R scripts
- path_utils.py: Defines project root and data paths for Python scripts 

results/
- SuSIE_PIP_Plots/: PIP plots from fine-mapping (01_Fine_Map_Coloc_eQTLs.R)
- SuSIE_Coloc_Results/: Union of credible sets that co-localized across tissues (01_Fine_Map_Coloc_eQTLs.R)
- aFCn_output/: Contains allelic fold change calculations  
  > allelic fold change was calculated using aFCn.py (https://github.com/PejLab/aFCn) on one thread, for the representative positions of each union set using the covariate adjusted expression files in RPSA_Covariate_Adjusted_Expression/,\ and the VCF files in GTEx_VCF_files/


---

Dependencies

Scripts written in R were run using R v4.3.3 and require following libraries (specific versions used in this analysis)
susieR (0.12.4)
Rfast (2.1.2)
Matrix (1.6.5)
coloc (5.2.3)
tidyverse (2.0.0)
gdata (3.0.1)
dplyr (1.1.4)
igraph ()
stringr () 
DESeq2 (1.42.1)
here (1.0.1)

Scripts are written in Python were run using python v3.11.9 and require the following libraries (specific versions used in this analysis):

numpy==1.26.0
pandas==2.1.1
statsmodels==0.14.4

---
Notes
- Raw genotype data (GTEx VCF) cannot be shared directly. Please obtained from dbGaP accession number phs000424v.8.p2
- The VCF files in GTEx_VCF_files/ were produced by subsetting the multi-sample VCF file provided in phs000424v.8.p2 to just the donors for the tissue, and to chr3:38406720-40412542 with bcftools.
- PLINK files in GTEx_PLINK_files_for_RPSA/ were produced by subsetting the multi-sample VCF file provided in phs000424v.8.p2 to just the donors for the tissue, and to chr3:38406720-40412542 with bcftools. The subset multi-sample VCF file was converted to a gene-dosage file using using PLINK1.9 with the commands --keep-allele-order and --recode A-transpose
- 01_Fine_Map_Coloc_eQTLs.R depends on the private gene-dosage files in GTEx_PLINK_files_for_RPSA/
- Calculating allelic fold change depends on the private vcf files in GTEx_VCF_files/
- Publicly available GTEx data was not uploaded due to size restrictions, and can be obtained from the provided URLs


