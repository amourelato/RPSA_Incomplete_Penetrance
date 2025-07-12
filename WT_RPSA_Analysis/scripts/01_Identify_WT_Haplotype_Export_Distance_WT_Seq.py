"""
Splits multi-sample phased VCF file into an A and B haplotype for each sample
Identifies the haplotype carrying the RPSA mutation for members of multiplex families
Uses manually phased long-read sequencing to identify the haplotype carrying the RPSA mutation for individuals with de-novo mutations
Exports the hamming distance between every possible pairwise combination of WT RPSA carrying haplotypes

Requirements:
- numpy==1.26.0
- pandas==2.1.1
- scipy==1.11.3

"""

import numpy as np
from itertools import product
from scipy.spatial.distance import pdist, squareform, hamming
import pandas as pd
from pathlib import Path
from path_utils import get_workflow_dir

def read_vcf(file_path):
    """Reads a VCF file into a Pandas DataFrame.

    Args:
        file_path (str): Path to the VCF file.

    Returns:
        pd.DataFrame: A pandas DataFrame representing the VCF file content.
    """
    with open(file_path, 'r') as f:
        num_header_lines = sum(1 for line in f if line.startswith('##'))

    df = pd.read_csv(file_path, sep='\t', skiprows=num_header_lines)
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)  # Rename the CHROM column

    return df

def SplitHaplotype(ID, vcf):
    """
    Splits a sample's phased genotype column (formatted as 'A|B') into two separate columns in a pandas DataFrame.

    Args:
        ID (str): Sample ID corresponding to a column in the VCF DataFrame.
        vcf (pd.DataFrame): DataFrame representing the VCF file content.

    Returns:
        pd.DataFrame: A DataFrame with two columns containing the individual phased alleles for the specified sample.
    """
    output = pd.DataFrame()
    output[['A', 'B']] = vcf[ID].str.split('|', expand=True)
    return output

def total_pairwise_hamming(haplotypes):
    """
    Computes the total pairwise Hamming distance among a set of haplotypes.

    Args:
        haplotypes (list of list of int): A list where each element is a haplotype represented as a list of alleles (0 and 1).

    Returns:
        float: The sum of Hamming distances for all unique haplotype pairs.
    """

    hap_array = np.array(haplotypes).astype(int)
    dist = pdist(hap_array, metric='hamming')
    return dist.sum()

def find_min_hamming_for_family(family, vcf):
    """
    Finds the sequence of the haplotype carrying the RPSA mutation (A or B) for all individuals in a family 
    by minimizing the total pairwise Hamming distance.

    Args:
        family (list of str): List of individual sample IDs belonging to the same family.
        vcf (pd.DataFrame): DataFrame representing the VCF file content.

    Returns:
        dict: A dictionary containing:
            - 'family': The input list of family member IDs.
            - 'min_hamming': The minimum total pairwise Hamming distance found.
            - 'assignment': A dict mapping each individual to 'A' or 'B' indicating the haplotype carrying the RPSA mutation.
    """

    #Split phased VCF A and B for each individual, and place A and B in a dictionary
    hap_dict = {}
    for indiv in family:
        split = SplitHaplotype(indiv, vcf)
        hap_dict[indiv] = (split['A'].tolist(), split['B'].tolist())
    
    #Generate all possible combinations of haplotype assignments (A or B) for all members of the same family
    best_combo = None
    min_dist = float('inf')
    individuals = list(hap_dict.keys())
    choices = list(product([0, 1], repeat=len(individuals)))  # 0 = A, 1 = B
    
    #Find the combination of haplotypes within a family that minimizes the total pairwise Hamming distance
    for combo in choices:
        haplotypes = []
        for i, choice in enumerate(combo):
            indiv = individuals[i]
            hap = hap_dict[indiv][choice]
            haplotypes.append(hap)
        
        dist = total_pairwise_hamming(haplotypes)
        
        if dist < min_dist:
            min_dist = dist
            best_combo = combo
    
    #Define which of A or B is the haplotype that carries the RPSA mutation for each member of a multiplex family
    assignment_dict = {indiv: ('A' if bit == 0 else 'B') for indiv, bit in zip(individuals, best_combo)}
    
    return {
        'family': family,
        'min_hamming': min_dist,
        'assignment': assignment_dict
    }

def compute_haplotype_distance_matrices(vcf, assignment, assignment_WT):
    """
    Computes pairwise Hamming distances between the haplotypes carrying the RPSA mutation, and between the haplotypes carrying the WT copy of RPSA
    for the samples in the given VCF DataFrame.

    Args:
        vcf (pd.DataFrame): DataFrame representing the VCF file content.
        assignment (dict): Dictionary mapping each sample ID to the sequence ('A' or 'B') carrying the mutation.
        assignment_WT (dict): Dictionary mapping each sample ID to the sequence ('A' or 'B') carrying the wild-type haplotype.

    Returns:
        dist_mutant_df (pd.DataFrame): DataFrame containing pairwise Hamming distances between mutant haplotypes.
        dist_WT_df (pd.DataFrame): DataFrame containing pairwise Hamming distances between wild-type haplotypes.
    """

    individuals = list(assignment.keys())
    
    #Create a list of RPSA mutation carrying haplotypes
    hap_mutant = []
    for indiv in individuals:
        A, B = vcf[indiv].str.split('|', expand=True).values.T
        hap = A if assignment[indiv] == 'A' else B
        hap_mutant.append(hap.astype(int))
    
    #Create a list of WT RPSA carrying haplotypes
    hap_WT = []
    for indiv in individuals:
        A, B = vcf[indiv].str.split('|', expand=True).values.T
        hap = A if assignment_WT[indiv] == 'A' else B
        hap_WT.append(hap.astype(int))

    #Convert both lists to arrays
    hap_mutant = np.array(hap_mutant)
    hap_WT = np.array(hap_WT)

    #Compute pairwise Hamming distances
    dist_mutant = squareform(pdist(hap_mutant, metric='hamming'))
    dist_WT = squareform(pdist(hap_WT, metric='hamming'))

    #Return pairwise Hamming distances as DataFrames with individuals as index/columns
    dist_mutant_df = pd.DataFrame(dist_mutant, index=individuals, columns=individuals)
    dist_WT_df = pd.DataFrame(dist_WT, index=individuals, columns=individuals)

    return dist_mutant_df, dist_WT_df

def ReturnHaplotypes(Mutant_Haplotype_Dictionary, vcf):
    """
    Splits a VCF DataFrame containing phased genotypes into two separate DataFrames:
    one for haplotypes carrying the mutant RPSA allele, and one for those carrying the wild-type RPSA allele.

    Args:
        vcf (pd.DataFrame): DataFrame representing the VCF file content.
        Mutant_Haplotype_Dictionary (dict): Dictionary mapping each sample ID to the sequence ('A' or 'B') carrying the mutation.

    Returns:
        NM_Haplotypes_All (pd.DataFrame): DataFrame containing haplotypes carrying the mutant allele.
        M_Haplotypes_All (pd.DataFrame): DataFrame containing haplotypes carrying the wild-type allele.
    """

    NM_Haplotypes_All = pd.DataFrame()
    M_Haplotypes_All = pd.DataFrame()

    for individual in Mutant_Haplotype_Dictionary.keys():

        Split_Haplotype = SplitHaplotype(individual, vcf) 
        
        #decide which haplotype to append to NM_Haplotypes_Dictionary based on which haplotype is the mutant haplotype
        NM_letter_to_append = ''
        M_letter_to_append = ''
        mutant_haplotype = Mutant_Haplotype_Dictionary[individual]
        if mutant_haplotype == "A":
            M_letter_to_append = "A"
            NM_letter_to_append = "B"

        elif mutant_haplotype == "B":
            M_letter_to_append = "B"
            NM_letter_to_append = "A"

        NM_Haplotypes_All[individual] = Split_Haplotype[NM_letter_to_append]    
        M_Haplotypes_All[individual] = Split_Haplotype[M_letter_to_append]

    #make add in the position column from the VCF file    
    NM_Haplotypes_All.insert(0, 'POS', vcf['POS'])
    M_Haplotypes_All.insert(0, 'POS', vcf['POS'])
    
    return NM_Haplotypes_All, M_Haplotypes_All

def main():
    """
    Main function

    This function:
    - Loads phased genotype data from a VCF file.
    - Loads family grouping information for multiplex kindreds
    - Loads the mutation carrying haplotype for individuals with denovo RPSA variants
    - Splits phased haplotypes for each individual.
    - Minimizes within-family Hamming distances to identify the RPSA WT and mutation carrying haplotypes
    - Outputs pairwise distances between WT carrying haplotypes and the sequences of the WT carrying haplotypes for all samples in the VCF file

    Returns:
        None
    """
    #set the workflow directory
    workflow_dir = get_workflow_dir()

    #read in pedigree file and create a list of family lists
    pedigree_path = workflow_dir / 'data' / 'Multiplex_Structure_and_Denovo_Haplotype_Assignment' / 'Individuals_Families.txt'
    pedigree_df = pd.read_csv(pedigree_path, sep="\t", dtype=str)
    family_lists = pedigree_df.groupby('Kindred')['Individual'].apply(list).tolist()
    
    #read in shapeit5 phased vcf, the actual phased VCF contains private genetic data and cannot be uploaded
    
    vcf_path = workflow_dir / 'data' / 'SHAPEIT5_Phased' / 'SHAPEIT5_pedigree_ref_phased_combined_wgs_chr3_filtered_25nt_all_indv_chr_3_38909631_39909631.bcf'
    SHAPEIT5 = read_vcf(vcf_path)
    
    #read in dataframe containing SHAPEIT5 WT and mutant haplotypes for individuals with denovo variants phased with read-backed phasing
    denovo_path = workflow_dir / 'data' / 'Multiplex_Structure_and_Denovo_Haplotype_Assignment' / 'Denovo_WT_Mutant.txt'
    denovo_df = pd.read_csv(denovo_path, sep='\t', dtype=str)
    
    #build the mutant haplotype dictionary for the individuals in multiplex families
    final_assignment_mutant_haplotype = {}
    final_assignment_WT_haplotype = {}

    for family in family_lists:
        if len(family) > 1:
            res = find_min_hamming_for_family(family, SHAPEIT5)
            final_assignment_mutant_haplotype.update(res['assignment'])
            
    #build the WT haplotype dictionary for the individuals in multiplex families
    for indiv, hap in final_assignment_mutant_haplotype.items():
        final_assignment_WT_haplotype[indiv] = 'B' if hap == 'A' else 'A'
    
    #add the denovo individuals with variants phased with read-backed phasing to the mutant and WT haplotype dictionaries
    for row in denovo_df.index.values.tolist():
        final_assignment_mutant_haplotype[denovo_df['Individual'][row]] = denovo_df['Mutant'][row]
        final_assignment_WT_haplotype[denovo_df['Individual'][row]] = denovo_df['WT'][row]
    
    dist_A, dist_B = compute_haplotype_distance_matrices(SHAPEIT5, final_assignment_mutant_haplotype, final_assignment_WT_haplotype)
    
    WT_Hap, Mut_Hap = ReturnHaplotypes(final_assignment_mutant_haplotype, SHAPEIT5)
    
    output_dir = workflow_dir / 'data' / 'Distance_Haplotype_Intermediate_Files'

    dist_B.to_csv(output_dir / 'Hamming_Distance_WT_Haplotypes.csv', sep=',')
    WT_Hap.to_csv(output_dir / 'WT_Haplotypes.csv', sep=',', index=False) #the WT haplotypes contain private genetic data and cannot be uploaded

if __name__ == "__main__":
    main()

   

