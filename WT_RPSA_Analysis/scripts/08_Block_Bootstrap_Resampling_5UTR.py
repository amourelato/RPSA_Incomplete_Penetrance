"""
Performs block bootstrap resampling to place a 95% CI on the chi-squared test statistic

Requirements:
- numpy==1.26.0
- pandas==2.1.1
- scipy==1.11.3
- matplotlib==3.8.4

"""

import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
from pathlib import Path
from path_utils import get_workflow_dir


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
    Re-samples blocks of individuals with replacement and calculates a chi-squared 
    test statistic from each re-sampling iteration.

    Args:
        Fam_Dict (dict): A dictionary mapping block names to lists of individual identifiers
                         within each block.
        Individual_DF (pd.DataFrame): A pandas DataFrame containing individual identifiers 
                                      and the blocks they belong to.
        iterations (int): Number of bootstrap iterations to perform.

    Returns:
        list: A list of chi-squared test statistics (one per bootstrap iteration).
    """

    random.seed(12345)
    Chi_Squared_Stat = []
 
    
    while len(Chi_Squared_Stat) < iterations:
    
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
        g = 0
        h = 0
    
        for indv in Indivdauls_from_blocks:
    
            Hap = Individual_DF["Haplotype"][indv]
            Phen = Individual_DF["Phenotype"][indv]
           
            if int(Phen) == 0:
                if Hap == "G0":
                    b += 1
                elif Hap == "G1":
                    d += 1
                elif Hap == "G2":
                    f += 1
                elif Hap == "G3":
                    h += 1
    
            elif int(Phen) == 1:
                if Hap == "G0":
                    a += 1
                elif Hap == "G1":
                    c += 1
                elif Hap == "G2":
                    e += 1
                elif Hap == "G3":
                    g += 1
    
        if a==0 and b==0:
            continue
    
        elif c==0 and d==0:
            continue
    
        elif e==0 and f==0:
            continue
    
        elif g==0 and h==0:
            continue
    

        else:
            table = table = np.array([[a, c, e, g],[b, d, f, h]])
            chi2, p, dof, expected = chi2_contingency(table)
            Chi_Squared_Stat.append(chi2)


    return Chi_Squared_Stat


def FreedmanDiaconisBins(List):
    """
    Calculate the number of histogram bins using the Freedman-Diaconis rule.

    Args:
        List (list): A list or array of numeric values.

    Returns:
        int: The recommended number of histogram bins.
    """

    q25, q75 = np.percentile(List, [25, 75])
    bin_width = 2 * (q75 - q25) * len(List) ** (-1/3)
    bins = round((max(List) - min(List)) / bin_width)
  
    return bins

def main():
    """
    Main function

    This function:
    - Loads the individual identifiers of the the sequenced individuals with 5UTR RPSA variants, their blocks, their phenotype, and their WT RPSA Group
    - Creates a dictionary mapping blocks to the individuals in the block
    - Obtains 100,000 bootstrap replicates of the chi-squared test statistic
    - Prints the mean, and 2.5 and 97.5 percentiles of the replicates of the chi-squared test statistic
    - Exports a histogram of the replicates of the chi-squared test statistic

    Returns:
        None
    """


    #set the workflow directory
    workflow_dir = get_workflow_dir()

    #read in Individual Datasheet
    Individual_DF_path = workflow_dir / 'data' / 'Individual_Blocks_Group_Phenotype' / 'Individual_Block_Group_Phenotype_5UTR.txt'
    Individual_DF = pd.read_csv(Individual_DF_path, sep="\t", index_col = "Individual")
     
    #Create the dictionary that maps blocks to the individuals in the blocks
    Fam_Dict = CreateFamilyDict(Individual_DF)
    
    # Obtain bootstrap replicates of the chi-squared test statistic
    Bootstrapped_Test_Statistic = BlockBootstrap(Fam_Dict, Individual_DF, 10**5)
    

    # Calculate the percentile 95% CI
    obs_p = np.sort(Bootstrapped_Test_Statistic)
    mean = np.mean(Bootstrapped_Test_Statistic)
    percentile_low = np.percentile(Bootstrapped_Test_Statistic, 2.5)
    percentile_high = np.percentile(Bootstrapped_Test_Statistic, 97.5)
    
    # Show the mean and 95% CI
    print(f"Mean:\t{mean}")
    print(f"0.025 Percentile:\t{percentile_low}")
    print(f"0.975 Percentile:\t{percentile_high}")
    
    
    # Plot the Bootstraped Chi-squared test statistics
    bins = FreedmanDiaconisBins(np.sort(Bootstrapped_Test_Statistic))
    plt.hist(Bootstrapped_Test_Statistic, density=True, bins=bins, color='#1338BE', edgecolor='black', linewidth=0.5)
    plt.axvline(x=percentile_low, color='black', linestyle='--', linewidth=1.5)
    plt.axvline(x=mean, color='black', linestyle='--', linewidth=1.5)
    plt.axvline(x=percentile_high, color='black', linestyle='--', linewidth=1.5)
    plt.axvline(x=5.991, color='red', linestyle='--', linewidth=1.5)
    
    # Remove y-axis
    plt.gca().axes.get_yaxis().set_visible(False)
    
    #export chi-squared histogram
    output_dir = workflow_dir / 'Results' / 'Block_Bootstrap_Resampling'
    plt.savefig(output_dir / "Block_Bootstrap_Chi_Squared_Histogram_5UTR.pdf", format='pdf', dpi=600, bbox_inches='tight')
 
if __name__ == "__main__":
    main()