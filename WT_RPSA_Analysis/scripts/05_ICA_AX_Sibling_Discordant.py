"""
Identifies regions that are not identical by state (IBS) on chr3 carrying the WT copy for PRSA for siblings in ICA-AX

Requirements:
- pandas==2.1.1

"""

import pandas as pd
import itertools
from path_utils import get_workflow_dir

def read_vcf(file_path):
    """Reads a VCF file into a Pandas DataFrame.

    Args:
        file_path (str): The path to the VCF file.

    Returns:
        pd.DataFrame: A Pandas DataFrame representing the VCF file content.
    """
    with open(file_path, 'r') as f:
        num_header_lines = sum(1 for line in f if line.startswith('##'))

    df = pd.read_csv(file_path, sep='\t', skiprows=num_header_lines)
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)  # Rename the CHROM column

    return df

def subset_genotype_vcf(vcf_unphased):
    """
    Extract genotype calls from VCF sample columns in a DataFrame.

    Args:
        vcf_unphased (pd.DataFrame): DataFrame read from a VCF file, containing 
                                     sample columns and other metadata columns.

    Returns:
        pd.DataFrame: Modified DataFrame where each sample column contains only 
                      unphased genotype calls (e.g., '0/1', '1/1').
    """
    for column in vcf_unphased.columns:
        if column.startswith("JM"):
            #Remove everything but the genotype information
            vcf_unphased[column] = vcf_unphased[column].str.replace("|", "/").str.split(":").str[0]
    
    return vcf_unphased

def merge_subset_vcf_MAF(vcf_unphased, MAF_CSV_File):
    """
    Merge unphased VCF data with allele frequency data and filter by MAF threshold (MAF > 0.03).

    Args:
        vcf_unphased (pd.DataFrame): DataFrame containing VCF variant data and sample genotypes.
        MAF_CSV_File (pd.DataFrame): DataFrame containing allele frequency data with at least
                                     the columns 'CHROM', 'POS', 'REF', 'ALT', and 'AF_total'.

    Returns:
        pd.DataFrame: Filtered merged DataFrame containing only variants with MAF > 0.03.
    """
    merged_df = pd.merge(vcf_unphased, MAF_CSV_File, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')
    merged_df_above_threshold = merged_df[merged_df['AF_total'] >= 0.03]

    return merged_df_above_threshold

def mendelian_phase(row, parent_col, child1_col, child2_col):
    """
    Perform Mendelian phasing for two children based on the genotype of a single parent.

    Args:
        row (pd.Series): A row from the dataframe containing genotypes for parent and children.
        parent_col (str): Column name containing the parent's genotype (e.g., '0/0', '1/1', '0/1').
        child1_col (str): Column name containing the first child's genotype.
        child2_col (str): Column name containing the second child's genotype.

    Returns:
        pd.Series: A Series with two elements:
            - child1_other (str or None): The allele carried by child1 on the haplotype
              not inherited from the parent (e.g., '0', '1', or None if phasing not possible).
            - child2_other (str or None): The allele carried by child2 on the haplotype
              not inherited from the parent.
    """

    parent_genotype = row[parent_col]
    child1_genotype = row[child1_col]
    child2_genotype = row[child2_col]
    
    #Initialize WT RPSA haplotype sequence
    child1_other = None
    child2_other = None

    #Only proceed if parent genotype is valid
    if parent_genotype == '0/0':
        
        if child1_genotype in ['0/0']:
            child1_other = '0'
        elif child1_genotype in ['0/1', '1/0']:
            child1_other = '1'
        
        if child2_genotype in ['0/0']:
            child2_other = '0'
        elif child2_genotype in ['0/1', '1/0']:
            child2_other = '1'

    elif parent_genotype == '1/1':
     
        if child1_genotype in ['1/1']:
            child1_other = '1'
        elif child1_genotype in ['1/0', '0/1']:
            child1_other = '0'
        
        if child2_genotype in ['1/1']:
            child2_other = '1'
        elif child2_genotype in ['1/0', '0/1']:
            child2_other = '0'

    # For parent == '0/1' or missing, leave as None as we cannot phase 
    return pd.Series(
        [child1_other, child2_other],
        index=['child1_WT', 'child2_WT']   
    )

def main():
    """
    Perform Mendelian phasing on selected samples and identify shared/discordant haplotype regions.

    1. Read in an unphased VCF file.
    2. Remove all columns except 'CHROM', 'POS', 'REF', 'ALT', and three specific sample columns.
    3. Subset sample columns to contain only genotype information (e.g., '0/1', '1/1').
    4. Read in a CSV file of population allele frequencies and rename relevant columns.
    5. Merge VCF data with MAF data and filter to variants with MAF ≥ 0.03.
    6. Further filter to retain only positions where the parent sample is homozygous (0/0 or 1/1).
    7. Perform Mendelian phasing to assign haplotypes to siblings based on the parent.
    8. Subset phased variants to positions that differ between the siblings on the haplotype inherited 
       from the non-phasing parent.
    9. Identify consecutive discordant positions and group those separated by less than 1,000,000 bases.
    10. For each such group, print the smallest and largest genomic positions.


    Returns:
        None
    """
        
    #set the workflow directory
    workflow_dir = get_workflow_dir()

    #Read in unphased vcf_file for all disc sib-pairs and their parent that passed on the Mutant RPSA copy for mendelian phasing
    vcf_unphased_path = workflow_dir / 'data' / 'unphased_Chr3_VCF' / 'chr3_DisSib_WGS_25nt_filtered.bcf' 
    vcf_unphased = read_vcf(vcf_unphased_path)
    #Remove all columns except the chromosome, position, ref, alt, and ICA-AX sample
    vcf_unphased = vcf_unphased.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'JM9401','JM9402','JM9403'], axis=1) 
    #Subset to just genotype information
    vcf_genotype = subset_genotype_vcf(vcf_unphased)
    
    #Read in MAF csv file and name the columns
    MAF_CSV_File_path = workflow_dir/ 'data' / 'GnomAD_v3_chr3_MAF' / 'chr03_GnomAD_v3.csv'
    MAF_CSV_File = pd.read_csv(MAF_CSV_File_path, sep=",")
    MAF_CSV_File.columns = ['CHROM', 'POS', 'REF',	'ALT','AC_total', 'AN_total', 'AF_total']
    
    #Merge the MAF and vcf files, and subset to just positions in the VCF with MAF > 0.03
    df = merge_subset_vcf_MAF(vcf_genotype, MAF_CSV_File)
    
    #Filter to keep only homozygous positions in the parent (0/0 or 1/1)
    df_homo_parent = df[df['JM9145'].str.contains('0/0|1/1')]
    
    #Perform the mendelian phasing
    df[['child1_WT', 'child2_WT']] = df.apply(
        mendelian_phase, axis=1,
        parent_col='JM9145', child1_col='JM9375', child2_col='JM9146'
    )
    
    #Subset to just positions that are different between the siblings on their WT RPSA allele
    diff_df = df[
        df['child1_WT'].notnull() &
        df['child2_WT'].notnull() &
        (df['child1_WT'] != df['child2_WT'])
    ]
    
    #Identify IBD segments
    diff_df = diff_df.copy()
    diff_df['diff'] = diff_df['POS'].diff().fillna(0)
    diff_df['group'] = (diff_df['diff'] >= 1000000).cumsum()
    
    IBD_boundaries = diff_df.groupby('group')['POS'].agg(start='min', end='max').reset_index()
    
    print(IBD_boundaries)

if __name__ == "__main__":
    main()





