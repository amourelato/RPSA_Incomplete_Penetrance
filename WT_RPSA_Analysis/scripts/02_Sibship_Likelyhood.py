"""
Uses a maximum likelihood binomial method to test the null hypothesis that 1Mb sequences flanking the WT copy of RPSA and spleen phenotype were independent in siblings

Requirements:
- numpy==1.26.0
- scipy==1.11.3

"""

import numpy as np
from scipy.stats import norm

def likelihood(alpha, n, S, x, y):
    """
    Binomial likelihood function

    Args:
        alpha (int): Parameter alpha (probability parameter).
        n (int): number of siblings in the sibship.
        S (int): number of siblings in a sibship with ICA.
        x (int): number of siblings with ICA who inherited A.
        y (int): number of siblings without ICA who inherited A.

    Returns:
        float: The computed likelihood value.
    """

    term1 = alpha**(x + n - S - y) * (1 - alpha)**(S - x + y)
    term2 = (1 - alpha)**(x + n - S - y) * alpha**(S - x + y)
    return 0.5 * (term1 + term2)


def total_log_likelihood(alpha, families):
    """
    Compute the total log-likelihood across multiple sibships.
    Evaluates the likelihood function for each sibship and sums the logarithms to compute the overall log-likelihood.

    Args:
        alpha (float): Probability parameter at which to evaluate the likelihood function.
                       Must be between 0 and 1 (exclusive).
        families (array-like): An array or list where each element represents a sibship
                               and contains the fields 'n', 'S', 'x', and 'y' as integers.

    Returns:
        float: The total log-likelihood across all sibships.
               Returns negative infinity (-np.inf) if alpha is outside (0,1)
               or if any individual likelihood is non-positive.
    """

    if alpha <= 0 or alpha >= 1:
        return -np.inf
    log_sum = 0.0
    for fam in families:
        n, S, x, y = fam['n'], fam['S'], fam['x'], fam['y']
        L = likelihood(alpha, n, S, x, y)
        if L <= 0:
            return -np.inf
        log_sum += np.log(L)
    return log_sum


def total_likelihood(alpha, families):
    """
    Compute the total likelihood across multiple sibships.

    Evaluates the likelihood function for each sibship and returns the product
    of the individual likelihoods.

    Args:
        alpha (float): Probability parameter at which to evaluate the likelihood function.
                       Must be between 0 and 1 (exclusive).
        families (array-like): An array or list where each element represents a sibship
                               and contains the fields 'n', 'S', 'x', and 'y' as integers.

    Returns:
        float: The total likelihood across all sibships.
               Returns 0 if alpha is outside (0,1) to avoid invalid computation.
    """
    
    if alpha <= 0 or alpha >= 1:
        return 0  # To avoid log(0) or division by zero
    product = 1.0
    for fam in families:
        n, S, x, y = fam['n'], fam['S'], fam['x'], fam['y']
        L = likelihood(alpha, n, S, x, y)
        product *= L
    
    return product

def estimate_alpha_and_pvalue():
    """
    Perform a grid search to estimate the alpha parameter that maximizes 
    the total log-likelihood across sibships, and compute a Z-score and p-value.

    The function:
    - Defines an array of sibship data with fields 'n', 'S', 'x', and 'y'. This definition is based on the hamming distances calculated between the siblings WT RPSA haplotypes
    - Evaluates the total log-likelihood over alpha values in [0.5, 1] with step 0.001.
    - Identifies the alpha that maximizes the log-likelihood (MLE).
    - Computes the Z-score comparing the likelihood at MLE alpha vs alpha=0.5.
    - Computes the corresponding one-sided p-value from the normal distribution.

    Returns:
        tuple: A tuple containing:
            - float: MLE estimate of alpha.
            - float: Maximum log-likelihood value.
            - float: Computed Z-score.
            - float: Computed p-value.
    """
    # List of family arrays
    families = [
    {'n': 2, 'S': 2, 'x': 2, 'y': 0}, #ICA-A
    {'n': 2, 'S': 2, 'x': 2, 'y': 0}, #ICA-H
    {'n': 3, 'S': 3, 'x': 3, 'y': 0}, #ICA-BV
    {'n': 3, 'S': 0, 'x': 0, 'y': 0}, #ICA-AG
    {'n': 2, 'S': 2, 'x': 2, 'y': 0}, #ICA-G
    {'n': 2, 'S': 1, 'x': 1, 'y': 0}, #ICA-AV
    {'n': 2, 'S': 1, 'x': 1, 'y': 0}, #ICA-BT
    {'n': 2, 'S': 1, 'x': 0, 'y': 0}, #ICA-AX
    {'n': 3, 'S': 1, 'x': 1, 'y': 0}, #ICA-AG
    {'n': 2, 'S': 1, 'x': 1, 'y': 0}, #ICA-AO
    {'n': 2, 'S': 1, 'x': 1, 'y': 0}, #ICA-AZ
    {'n': 2, 'S': 1, 'x': 1, 'y': 0}, #ICA-AZ
    {'n': 2, 'S': 1, 'x': 0, 'y': 0}, #ICA-AS
    ]

    # Range of alpha values to evaluate
    alphas = np.arange(0.5, 1.001, 0.001)
    log_likelihoods = [total_log_likelihood(a, families) for a in alphas]

    # Find index of max log-likelihood
    max_index = np.argmax(log_likelihoods)
    mle_alpha = alphas[max_index]
    max_logL = log_likelihoods[max_index]

    # Compute Z-score and p-value
    likelihood_mle = total_likelihood(mle_alpha, families)
    likelihood_null = total_likelihood(0.5, families)
    Z_score = np.sqrt(2 * np.log(likelihood_mle / likelihood_null))
    pval = norm.sf(Z_score)

    return mle_alpha, max_logL, Z_score, pval

if __name__ == "__main__":
    mle_alpha, max_logL, Z_score, pval = estimate_alpha_and_pvalue()
    print(pval)