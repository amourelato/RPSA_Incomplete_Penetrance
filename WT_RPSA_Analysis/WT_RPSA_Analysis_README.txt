
This directory contains scripts, data, and supporting files to reproduce WT RPSA haplotype analysis of the manuscript.  


Some required input files (phased and unphased VCFs) contain private genetic data and are not uploaded to GitHub. See below for details on requesting access.

---
Directory structure

data/
Contains raw input files and intermediate data:
- SHAPEIT5_Phased/ Private phased genotype data from the sequenced individuals in the ICA-RPSA cohort (not uploaded). See notes below                  
- Unphased_chr3_VCF/ Private phased genotype data from the sequenced members of kindreds ICA-AX and ICA-AS in the ICA-RPSA cohort (not uploaded). See notes below    
- Multiplex_Structure_and_Denovo_Haplotype_Assignment/ helper text files for 01_Identify_WT_Haplotype_Export_Distance_WT_Seq.py to indicate samples that belong to the same kindred, and haplotypes corresponding to the WT RPSA carrying sequence for the individuals with denovo mutations.  
- Individual_Blocks_Group_Phenotype/ helper text files for 07_Block_Bootstrap_Resampling_Coding.py, 08_Block_Bootstrap_Resampling_5UTR.py, 09_Block_Bootstrap_Resampling_Odds_Ratio_Coding.py, 10_Block_Bootstrap_Resampling_Odds_Ratio_5UTR.py to indicate blocks of samples 
- GnomAD_v3_chr3_MAF/ chromosome 3 minor allele frequencies obtained from GnomAD v3, required for 05_ICA_AX_Sibling_Discordant.py and 06_ICA_AS_Sibling_Discordant.py, not uploaded due to GitHub size restrictions
- Distance_Haplotype_Intermediate_Files/ Intermediate hamming distance and WT RPSA haplotype sequence files produced by by script 01_Identify_WT_Haplotype_Export_Distance_WT_Seq.py (WT RPSA haplotype sequences are not uploaded as they contain private genotype data)




Scripts/
Scripts to run the workflow. All scripts are written in Python numbered in execution order:
- 01_Identify_WT_Haplotype_Export_Distance_WT_Seq: Identifies shared WT RPSA carrying haplotype and calculates WT RPSA hamming distance matrix, depends on phased vcf file in SHAPEIT5_Phased/ which is not uploaded as it contains private genotype data which 
- 02_Sibship_Likelihood.py: Computes maximum likelihood binomial p-value using the sibships in the ICA-RPSA cohort
- 03_Grouping_WT_at_rep_pos.py: Groups individuals by specific haplotype positions, depends on WT RPSA haplotype sequences in Distance_Haplotype_Intermediate_Files/ which are not uploaded as they contain private genotype data
- 04_Plotting_Base_Distance_WT_G3.py: Plots base distances between G3 WT RPSA hamming distances, depends on hamming distance file in Distance_Haplotype_Intermediate_Files/
- 05_ICA_AX_Sibling_Discordant.py : Finds non-identical by state (IBS) segments on chr3 between siblings in ICA-AX, depends on vcf file in Unphased_chr3_VCF/, which is not uploaded as it contains private genotype data
- 06_IBS_Segments_B.py : Finds non-identical by state (IBS) segments on chr3 between siblings in ICA-AS, depends on vcf file in Unphased_chr3_VCF/, which is not uploaded as it contains private genotype data
- 07_Block_Bootstrap_Resampling_Coding.py : Block-bootstrap replicates of chi-squared statistic from individuals with coding RPSA mutations, depends on blocks in Individual_Blocks_Group_Phenotype/
- 08_Block_Bootstrap_Resampling_5UTR.py : Block-bootstrap replicates of chi-squared statistic from individuals with 5'UTR RPSA mutations, depends on blocks in Individual_Blocks_Group_Phenotype/
- 09_Block_Bootstrap_Resampling_Odds_Ratio_Coding.py: Block-bootstrap replicates of the odds ratio from individuals with coding RPSA mutations, depends on blocks in Individual_Blocks_Group_Phenotype/
- 10_Block_Bootstrap_Resampling_Odds_Ratio_5UTR.py: Block-bootstrap replicates of the odds ratio from individuals with coding RPSA mutations, depends on blocks in Individual_Blocks_Group_Phenotype/
- path_utils.py: Defines project root and data paths for Python scripts 


Results/
- WT_RPSA_Grouping/: WT RPSA group assignment based on the representative positions of each credible set (03_Grouping_Individuals.py)
- WT_RPSA_G3_Base_Distance/: Base distance between the WT RPSA haplotypes of the members of G3 (04_Plotting_Base_Distance.py)
- ChiSquared_Bootstrap/: Histogram of bootstrap replicates of the chi-squared test statistic (07_ChiSquared_Bootstrap.py and 08_ChiSquared_Figure.py)


---

Dependencies

Scripts are written in Python were run using python v3.11.9 and require the following libraries (specific versions used in this analysis):

* numpy==1.26.0
* pandas==2.1.1
* scipy==1.11.3
* matplotlib==3.8.4
* seaborn==0.13.2


---
Notes

- Some input data (phased and unphased VCF files) contain private genotype data and are not uploaded. Please contact the corresponding authors (casanova@rockefeller.edu; bebo283@rockefeller.edu) to request access to private genetic data
- Raw WGS reads were processed and aligned to hg38 using the Illumina DRAGEN germline pipeline v4.0.3, and a multi-sample vcf file for all sequenced individuals in the ICA-RPSA cohort was produced using GATK v3.6 GenotypeGVCFs. Quality control (QC) was performed on the multisample VCF to remove low-quality variants with Plink 1.9 and VCFtools: Genotypes with a ratio of reads for the least covered allele (reference or variant allele) to the total number of reads covering the position (minor read ratio, MRR) of <20% for heterozygous calls and/or a GATK genotype quality (GQ) < 20% and/or a depth of coverage (DP) < 8x were considered missing. We excluded sites: i) with a call rate < 95%, (ii) with a p-value for departure from Hardy-Weinberg equilibrium < 1x10-6, (iii) located in low-complexity or decoy regions, (iv) multiallelic with more than four alleles, (v) spanning more than 25 nucleotides.
- The quality-controlled multi-sample VCF file of chromosome 3 was phased using SHAPEIT5 (v5.1.1) (https://github.com/odelaneau/shapeit5) phase_common with a genetic map of chr3 (https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b38.tar.gz), pedigree information provided by a custom text file reflecting the familial relationships in our cohort, and chr3 reference haplotypes obtained from the 1000 Genomes project (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel.vcf.gz). SHAPEIT5 phase_common was run on one thread with seed = 123456, mcmc-iterations = 10b,1p,1b,1p,1b,1p,1b,1p,1m, and pbwt-depth = 8. All other parameters were default. After performing quality control and haplotype phasing, the phased multi-sample VCF for all sequenced individuals in the cohort was subset to a 1MB window (chr3:38,909,631-39,909,631; hg38) surrounding the RPSA locus with bcftools. This phased and subset vcf file was used as input for 01_Identify_WT_Haplotype_Export_Distance_WT_Seq.py 
-MAF file for chr3 was not uploaded due to GTEx size restrictions.

