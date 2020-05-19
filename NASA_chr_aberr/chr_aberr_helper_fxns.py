# pandas/numpy for handling data
import pandas as pd
import numpy as np

# seaborn/matplotlib for graphing
import matplotlib.pyplot as plt
import seaborn as sns

# statistics
from statistics import mean 
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats
import scikit_posthocs as sp
from statsmodels.stats.anova import AnovaRM

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.libqsturng import psturng


def combine_midflight(row):
    if 'mid-flight 1' in row or 'mid-flight 2' in row:
        row = 'mid-flight'
        return row
    else:
        return row
    

def scipy_anova_post_hoc_tests(df=None, flight_status_col='flight status new',
                               sig_test=stats.f_oneway):
    """
    df should be melted by aberration type
    """
    # make list of aberrations
    aberrations = list(df['aberration type'].unique())
    
    # loop through aberrations & perform anovas between pre/mid/post
    for aberr in aberrations:
        
        g_1 = df[(df[flight_status_col] == 'Pre-Flight') & (df['aberration type'] == aberr)]['count per cell']
        g_2 = df[(df[flight_status_col] == 'Mid-Flight') & (df['aberration type'] == aberr)]['count per cell']
        g_3 = df[(df[flight_status_col] == 'Post-Flight') & (df['aberration type'] == aberr)]['count per cell']
        statistic, p_value = sig_test(g_1, g_2, g_3)
        print(aberr, p_value)

        # if anova detects sig diff, perform post-hoc tests            
        if p_value <= 0.05:
            mc = MultiComparison(df[df['aberration type'] == aberr]['count per cell'], 
                                 df[df['aberration type'] == aberr][flight_status_col])
            mc_results = mc.tukeyhsd()
            print(mc_results)
            res = mc_results
            print(f'pvalues: {list(psturng(np.abs(res.meandiffs / res.std_pairs), len(res.groupsunique), res.df_total))}')
            print('\n')