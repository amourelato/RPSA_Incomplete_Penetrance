"""
Corrects for significant linear effects of significant covariates

Requirements:
- numpy==1.26.0
- pandas==2.1.1
- statsmodels==0.14.4
"""

import statsmodels.api as sm
import numpy as np
import pandas as pd
from pathlib import Path
from path_utils import get_workflow_dir

#define a function that reads in the expression file and writes out pandas df with re-formatted donor names subset to just RPSA expression
def ReadExpressionFile(file_path):
    """
    Reads in normalized and log transformed expression for the tissue, and subsets to RPSA expression

    Args:
        file_path : A path to the normalized and covariate adjusted csv file for the tissue
        Expression_RPSA (pd.Dataframe): A pandas dataframe containing normalized and log transformed RPSA expression for the tissue

    Returns:
        pd.DataFrame: A pandas DataFrame with normalized and log transformed RPSA expression
    """

    #read in expression file and subset to just RPSA
    Expression = pd.read_csv(file_path, sep=",")
    Expression.columns.values[0] = 'gene_id'
    Expression_RPSA = Expression[Expression['gene_id'].str.startswith("ENSG00000168028")]

    #re-format donor names to match covariate file donor names
    new_Expression_colnames = [column.split(".")[0] +"-"+ column.split(".")[1] for column in Expression_RPSA if column.startswith("GTEX")]
    new_Expression_colnames.insert(0, 'gene_id')
    Expression_RPSA.columns = new_Expression_colnames
    
    return Expression_RPSA

def CovCorrection(Covariates, Expression_RPSA):
    """Takes covariates and normalized RPSA expression, and perfomrs covariate correction

    Args:
        Covariates (pd.DataFrame): A pandas dataframe containing covariates for the tissue
        Expression_RPSA (pd.Dataframe): A pandas dataframe containing normalized and log transformed RPSA expression for the tissue

    Returns:
        pd.DataFrame: A pandas DataFrame containing the covariate correctred RPSA expression
    """


    #get a list of donors that overlap between covariate and expression 
    donors = list(set([column for column in Expression_RPSA if column.startswith("GTEX")]) & set([column for column in Covariates if column.startswith("GTEX")]))

    
    RPSA_Expression_Array = Expression_RPSA[donors].iloc[0].to_numpy().astype(float)
    Cov_array = Covariates[donors].to_numpy().T
    Cov_array_with_intercept = sm.add_constant(Cov_array)

    #solve regression coefficient with ordinary least squares
    model = sm.OLS(RPSA_Expression_Array, Cov_array_with_intercept).fit()

    #now keep just the covariates with a p val < 0.01, we will only regress significant covariates
    mask = model.pvalues < 0.01
    Cov_array_subset = Cov_array_with_intercept.T[mask]

    #add intercept of one to the subset covariate array 
    Cov_array_subset_with_intercept = sm.add_constant(Cov_array_subset.T)

    #solve regression coefficient using just significant covariates
    model_subset = sm.OLS(RPSA_Expression_Array, Cov_array_subset_with_intercept).fit()

    #predict expression from regression coefficients
    y_pred = model_subset.predict(Cov_array_subset_with_intercept)

    #find the residuals
    residuals = RPSA_Expression_Array - y_pred

    #turn the residuals into a pandas dataframe and add the gene name for aFCn
    residuals_reshaped = residuals.reshape(1, -1)
    RPSA_Residual_Exp = pd.DataFrame(residuals_reshaped)
    RPSA_Residual_Exp.insert(0, 'Name', "ENSG00000168028")

    #add the donor names as the colnames 
    donors.insert(0, "Name")
    RPSA_Residual_Exp.columns = donors

    #return the residuals for aFCn
    return RPSA_Residual_Exp

def main():
    """
    Main function

    This function:
    - Loads covariates GTEx identified for each tissue
    - Loads normalized and log transformed RPSA expression for each tissue
    - Performs covariate correction for significant linear effects of significant covariates
    - Outputs covariate corrected RPSA expression

    Returns:
        None
    """
    #set the workflow directory, and define the covariate and normalized
    workflow_dir = get_workflow_dir()
    cov_dir = workflow_dir / 'data' / 'GTEx_Analysis_v8_eQTL_covariates'
    norm_exp_dir = workflow_dir / 'data' / 'GTEx_Normalized_Tissue_RNA_counts'

    #Initialize an empty dictionary to store dataframes and read in covariates as dataframes
    covariate_dfs = {}
    for file_path in cov_dir.glob("*.txt"):
        name = file_path.stem.split(".v8")[0]
        df = pd.read_csv(file_path, sep="\t")
        covariate_dfs[name] = df

    #Initialize an empty dictionary to store dataframes and read in normalized expression as dataframes
    normalized_expression_dfs = {}
    for file_path in norm_exp_dir.glob("*.csv"):
        name = file_path.stem.split(".csv")[0]
        df = ReadExpressionFile(file_path)
        normalized_expression_dfs[name] = df

    #Initialize an empty dictionary to store covariate corrected expression, and correct expression
    corrected_dfs = {}
    for key in covariate_dfs.keys():
        df = CovCorrection(covariate_dfs[key], normalized_expression_dfs[key])
        corrected_dfs[key] = df

    #Define the output directory and export the covariate adjusted RPSA expression for each tissue
    output_dir = workflow_dir / 'data' / 'RPSA_Covariate_Adjusted_Expression'
    for key, df in corrected_dfs.items():
        output_path = output_dir / f"{key}_RPSA_Residual_Exp.csv"
        df.to_csv(output_path, sep=",", index=False)

if __name__ == "__main__":
    main()