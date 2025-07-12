"""
Performs block bootstrap resampling to place a 95% CI on the Odds-Ratio

Requirements:
- numpy==1.26.0
- pandas==2.1.1

"""

import random
import pandas as pd
import numpy as np
from pathlib import Path
from path_utils import get_workflow_dir

#define a function that creates a dictionary of families from a pedigree file
def CreateFamilyDict(Individual_DF):
    """
    Creates a dictionary mapping blocks to lists of individuals

    Args:
        Individual_DF (pd.DataFrame): A pandas DataFrame with the individual identifiers and the block to which that individual belongs
                                      
    Returns:
        Family_Dict (dict): A dictionary where keys are blocks names and values are lists of
              individual identifiers belonging to the block
    """

    Family_Dict = {}
 
    for row in Individual_DF.index.values.tolist():
       
        key = Individual_DF["Block"][row]
        value = row
 
        if key in Family_Dict.keys():
            Family_Dict[key].append(value)
        else:
            Family_Dict[key] = [value]
 
    return Family_Dict

def BlockBootstrap(Fam_Dict, Individual_DF, iterations):
    """
    Re-samples blocks of individuals with replacement and calculates an odds ratio after grouping individuals in G1 and G2 together  
    (applying a Haldane-Anscombe correction if the table has a 0 count) from each re-sampling iteration.

    Args:
        Fam_Dict (dict): A dictionary mapping block names to lists of individual identifiers
                         within each block.
        Individual_DF (pd.DataFrame): A pandas DataFrame containing individual identifiers 
                                      and the blocks they belong to.
        iterations (int): Number of bootstrap iterations to perform.

    Returns:
        list: A list of odds-ratios (one per bootstrap iteration).
    """

    random.seed(12345)
    Odds_Ratio = []
 
   
    while len(Odds_Ratio) < iterations:
    
        Indivdauls_from_blocks = []
    
        for _ in range(0, len(Fam_Dict.keys())):
            Block = random.choice(list(Fam_Dict.values()))
            Indivdauls_from_blocks.extend(Block)
    
        a = 0
        b = 0
        c = 0
        d = 0
        e = 0
        f = 0
 
        for indv in Indivdauls_from_blocks:
    
            Hap = Individual_DF["Haplotype"][indv]
            Phen = Individual_DF["Phenotype"][indv]
           
            if int(Phen) == 0:
                if Hap == "G1":
                    b += 1
                elif Hap == "G2":
                    d += 1
                elif Hap == "G3":
                    f += 1
    
            elif int(Phen) == 1:
                if Hap == "G1":
                    a += 1
                elif Hap == "G2":
                    c += 1
                elif Hap == "G3":
                    e += 1
    
        #define a new table based on the bootstrap replicate, grouping G1 and G2 together
        #the table will be a 2x2 Rows are Group (G1+G2/G3), Columns are Spleen Phenotype (ICA or Spleen)
        
        table = np.array([[(b+d), (a+c)],[f, e]])

        #If the table has a 0 count, Apply Haldane-Anscombe correction and calculate an odds ratio based on the corrected table
        if np.any(table == 0):
            HA_table = table + 0.5
            HA_a, HA_b = HA_table[0]
            HA_c, HA_d = HA_table[1]
            HA_or = (HA_a * HA_d) / (HA_b * HA_c)
            Odds_Ratio.append(HA_or)
        #If the table does not have a 0 count, directly calculate an odds ratio
        else:
            Sum_a, Sum_b = table[0]
            Sum_c, Sum_d = table[1]
            Sum_or = (Sum_a * Sum_d) / (Sum_b * Sum_c)
            Odds_Ratio.append(Sum_or)


    return Odds_Ratio

def main():
    """
    Main function

    This function:
    - Loads the individual identifiers of the the sequenced individuals with coding RPSA variants, their blocks, their phenotype, and their WT RPSA Group
    - Creates a dictionary mapping blocks to the individuals in the block
    - Obtains 100,000 bootstrap replicates of the odds ratio
    - Prints the mean, and 2.5 and 97.5 percentiles of the replicates of the odds ratio

    Returns:
        None
    """

    #set the workflow directory
    workflow_dir = get_workflow_dir()

    #read in Individual Dataframe
    Individual_DF_path = workflow_dir / 'data' / 'Individual_Blocks_Group_Phenotype' / 'Individual_Block_Group_Phenotype_Coding.txt'
    Individual_DF = pd.read_csv(Individual_DF_path, sep="\t", index_col = "Individual")
     
    #Create the dictionary that maps blocks to the individuals in the blocks
    Fam_Dict = CreateFamilyDict(Individual_DF)
    
    # Obtain bootstrap replicates of the chi-squared test statistic
    Bootstrapped_Odds_Ratio = BlockBootstrap(Fam_Dict, Individual_DF, 10**5)
    

    # Calculate the percentile 95% CI
    obs_p = np.sort(Bootstrapped_Odds_Ratio)
    mean = np.mean(Bootstrapped_Odds_Ratio)
    percentile_low = np.percentile(Bootstrapped_Odds_Ratio, 2.5)
    percentile_high = np.percentile(Bootstrapped_Odds_Ratio, 97.5)
    
    # Show the mean and 95% CI
    print(f"Mean:\t{mean}")
    print(f"0.025 Percentile:\t{percentile_low}")
    print(f"0.975 Percentile:\t{percentile_high}")
    


if __name__ == "__main__":
    main()