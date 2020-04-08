# enables access to directories/files
import os

# for handling data
import numpy as np
from numpy import array
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile

# graphing
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import seaborn as sns

# statistics
from statsmodels.graphics.gofplots import qqplot
from scipy import stats
import scikit_posthocs as sp
from scipy.stats import zscore
from scipy.stats import ks_2samp
from statistics import mean 
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats
from scipy.stats import zscore
from scipy.stats import ks_2samp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
import scikit_posthocs as sp
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.libqsturng import psturng

import re
from ast import literal_eval
import more_itertools
import math

from matplotlib import lines
from matplotlib.offsetbox import AnchoredText
import imgkit

# machine learning
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from xgboost import XGBRegressor
from sklearn.metrics import explained_variance_score, r2_score
from sklearn.metrics import median_absolute_error
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import cross_val_score
from sklearn import preprocessing
import scipy.cluster.hierarchy as hac
import matplotlib.gridspec as gridspec
import random
import six

from sklearn.preprocessing import LabelEncoder


def generate_dictionary_for_telomere_length_data(patharg):
  
    """
    USAGE:
    telomere_data_dict = generate_dictionary_for_telomere_length_data(directory)
    Where the directory contains images of files containing telomere length data in
    a predefined format. This function is written specifically for the Excel file templates
    that I use, and will provide in this repository, but could be altered for any format.
   
    The individual telomere lengths column is extracted, cleansed of missing values & DAPI-intensity 
    values; outliers (3 std devs from mean of column) are removed; and the telomere length values are 
    standardized to each other by use of fluorescent beads which calibrate according to inherent 
    differences between microscope imaging sessions. The individual's ID & timepoint (from filename) (KEY) 
    is associated with its respective individual telomere length data (VALUE) as a KEY:VALUE pair 
    in the dictionary. The dictionary can then be looped over to initialize all timepoint data
    for that individual for analysis, i.e visualizations, statistics, etc.
    """
    
    # initialize dictionary to hold our data
    dict_astro_individ_telos_dfs = {}

    # loop through directory to grab files
    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
            print(f'{file.name} telomere data acquisition in progress..')
        
            try:
                df = pd.read_excel(file)

            except:
                print(f'{file.name} File not found..')
                return -1

            df.rename(columns={'Unnamed: 3':'Individ Telos'}, inplace=True)
            
            # these numbers correspond to rows containing information about the DAPI counterstain, NOT telomeres, so we drop
            DAPI_values_to_drop=[5, 192, 379, 566, 753, 940, 1127, 1314, 1501, 1688, 1875, 2062,
                    2249, 2436, 2623, 2810, 2997, 3184, 3371, 3558, 3745, 3932, 4119, 4306, 4493, 
                    4680, 4867, 5054, 5241, 5428]

            # grabbing individual telomere length data from the file & dropping DAPI info
            individual_telos_lengths = (df['Individ Telos'])
            individual_telos_lengths = individual_telos_lengths.drop(labels=DAPI_values_to_drop)
            
            # first pass at generating synthetic data for github exposition; to initialize actual
            # data, comment out the line below, and uncomment the .iloc[] line
#             individual_telos_lengths = individual_telos_lengths.sample(2500, random_state=1)
            individual_telos_lengths = individual_telos_lengths.iloc[7:5611]

            # ensure the telomere measurements are a numeric data type, drop any missing values, 
            # make data into a dataframe
            telos_str_toNaN = pd.to_numeric(individual_telos_lengths, errors='coerce')
            individual_telos_cleaned = telos_str_toNaN.dropna(axis=0, how='any')
            telos_df = individual_telos_cleaned.to_frame(name=None)
            
            # remove any telomere measurements that lie beyond 3 standard deviations of the mean
            # the data is relatively normal in shape, & this process removes about ~10-20 telos from ~5520
            # modest loss, acceptable to help standardize
            telos_individ_df = telos_df[(np.abs(stats.zscore(telos_df)) < 3).all(axis=1)]
            
            # logic clauses for recognizing which astronaut ID is in the sample name
            # different astronauts were imaging at different times and thus associated with 
            # different Cy3 calibrations for the microscope, thus data is standardized according to Cy3
            
            if ('5163' in file.name) or ('1536' in file.name):
                telos_individ_df_cy3Cal = telos_individ_df.div(59.86)

            elif '2171' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(80.5)

            elif '7673' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(2.11)

            elif '2479' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(2.18)

            elif '1261' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(2.16)

            else:
                telos_individ_df_cy3Cal = telos_individ_df
            
            #for publications & nasa data request 
            #average of all cy3 calibrated control telo measurements (11 age matched controls)
#             telos_individ_df_cy3Cal = telos_individ_df_cy3Cal.div(116.1848153)

            file_name_trimmed = file.name.replace('.xlsx', '')
            dict_astro_individ_telos_dfs[file_name_trimmed] = telos_individ_df_cy3Cal

    print('Done collecting all astronaut telomere length excel files')
    return dict_astro_individ_telos_dfs



def astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone, axsNUMtwo):

        astroDF = astroDF.to_numpy()
        astroquartile = astroquartile.to_numpy()

        N, bins, patches = axs[axsNUMone,axsNUMtwo].hist(astroDF, bins=n_bins, range=(0, 400), edgecolor='black')

        for a in range(len(patches)):
            if bins[a] <= np.quantile(astroquartile, 0.25):
                patches[a].set_facecolor('#fdff38')

            elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
                patches[a].set_facecolor('#d0fefe')

            elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
                patches[a].set_facecolor('#d0fefe')

            elif bins[a] > np.quantile(astroquartile, 0.75): 
                patches[a].set_facecolor('#ffbacd')


        axs[axsNUMone,axsNUMtwo].set_title(f"Histogram of {astroname}'s Telomeres")
        axs[axsNUMone,axsNUMtwo].set_xlabel('Bins of Individ. Telomeres')
        axs[axsNUMone,axsNUMtwo].set_ylabel('Freqs of Individ. Telomeres')
        axs[axsNUMone,axsNUMtwo].xaxis.set_major_locator(plt.MaxNLocator(12))
        

        
def astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone):

    astroDF = astroDF.to_numpy()
    astroquartile = astroquartile.to_numpy()


    N, bins, patches = axs[axsNUMone].hist(astroDF, bins=n_bins, range=(0, 400), edgecolor='black')

    for a in range(len(patches)):
        if bins[a] <= np.quantile(astroquartile, 0.25):
            patches[a].set_facecolor('#fdff38')

        elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
            patches[a].set_facecolor('#d0fefe')

        elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
            patches[a].set_facecolor('#d0fefe')

        elif bins[a] > np.quantile(astroquartile, 0.75): 
            patches[a].set_facecolor('#ffbacd')


    axs[axsNUMone].set_title(f'Histogram of Individual Telomeres for {astroname}')
    axs[axsNUMone].set_xlabel('Bins of Individ. Telomeres')
    axs[axsNUMone].set_ylabel('Freqs of Individ. Telomeres')
    axs[axsNUMone].xaxis.set_major_locator(plt.MaxNLocator(19))
    
    
def gen_missing_values_andimpute_or_randomsampledown(n_cells, telosPercell, astro_df, option=None):

    if astro_df.size > 5520:
        astro_dfsampled = astro_df.sample(5520)
        return astro_dfsampled

    if astro_df.size > 25 and astro_df.size <= 2760:
        missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
        rsampled = astro_df.sample(missing_data_difference, replace=True, random_state=28)
        concat_ed = pd.concat([rsampled, astro_df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed

    if astro_df.size > 25 and astro_df.size < 5520:
        missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )

        if option == 'rsamp':
            rsampled = astro_df.sample(missing_data_difference, random_state=28)
            concat_ed = pd.concat([rsampled, astro_df], sort=False)

            np.random.shuffle(concat_ed.to_numpy())
            concat_ed.reset_index(drop=True, inplace=True)
            return concat_ed
        else:
            return astro_df

    else:
        return astro_df
    
    
def statistics_between_timepoints(astro_pre, astro_mid1, astro_mid2, astro_post, 
    astro_prename, astro_mid1name, astro_mid2name, astro_postname, test):

    print(  astro_prename + '  vs  ' + astro_mid1name,
            test(astro_pre, astro_mid1), '\n',

            astro_prename + '  vs  ' + astro_mid2name,
            test(astro_pre, astro_mid2),'\n', 
            
            astro_mid1name + '  vs  ' + astro_postname,
            test(astro_mid1, astro_post),'\n', 

            astro_mid1name + '  vs  ' + astro_mid2name,
            test(astro_mid1, astro_mid2),'\n', 

            astro_mid2name + '  vs  ' + astro_postname,
            test(astro_mid2, astro_post),'\n', 

            astro_prename + '  vs  ' + astro_postname,
            test(astro_pre, astro_post),'\n', )


    
def statistics_between_timepoints_prepost_only(astro_pre, astro_post, astro_prename, astro_postname):

    print(astro_prename + '  compared vs  ' + astro_postname,
            stats.mannwhitneyu(astro_pre, astro_post),'\n', )
    
    
    
def get_astro_number_from_id(astro_id):
    astro_num = ''
    
    if astro_id == '5163':
        astro_num = 1
        synth = 'synthetic 1'
        
    elif astro_id == '1536':
        astro_num = 2
        synth = 'synthetic 2'
        
    elif astro_id == '7673':
        astro_num = 3
        synth = 'synthetic 3'
        
    elif astro_id == '2479':
        astro_num = 4
        synth = 'synthetic 4'
        
    elif astro_id == '2171':
        astro_num = 5
        synth = 'synthetic 5'
    
    elif astro_id == '1261':
        astro_num = 7
        synth = 'synthetic 7'
    
    elif astro_id == '3228':
        astro_num = 8
        synth = 'synthetic 8'
        
    elif astro_id == '2381':
        astro_num = 9 
        synth = 'synthetic 9'
        
    elif astro_id == '4819':
        astro_num = 10
        synth = 'synthetic 10'
        
    elif astro_id == '1062':
        astro_num = 11
        synth = 'synthetic 11'
        
    elif astro_id == '2494':
        astro_num = 12
        synth = 'synthetic 12'
        
    return astro_num, synth



def relative_flight_timepoint(name_key):
    if 'L' in name_key:
        flight_status = 'Pre-Flight'
    elif 'FD' in name_key:
        flight_status = 'Mid-Flight'
    elif 'R' in name_key:
        flight_status = 'Post-Flight'
        
    return flight_status


def quartile_cts_rel_to_df1(df1, df2):
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    
    quartile_1 = df2[df2 <= df1.quantile(0.25)].count()
    
    quartile_2_3 = df2[(df2 > df1.quantile(0.25)) & (df2 < df1.quantile(0.75))].count()

    quartile_4 = df2[df2 >= df1.quantile(0.75)].count()
    
    return quartile_1.values, quartile_2_3.values, quartile_4.values


def get_timepoint(name_key):
    timepoint_5_char = ['L-270', 'L-180', 'FD140', 'FD260', 'R+105', 'R+180', 'R+270']
    timepoint_4_char = ['L-60', 'FD45', 'FD90', 'R+60']
    timepoint_3_char = ['R+5', 'R+7'] 
    
    for timepoint in timepoint_5_char:
        if timepoint in name_key:
            timepoint = name_key[-5:]
            return timepoint.strip()
    
    for timepoint in timepoint_4_char:
        if timepoint in name_key:
            timepoint = name_key[-4:]
            return timepoint.strip()
            
    for timepoint in timepoint_3_char:
        if timepoint in name_key:
            timepoint = name_key[-3:]
            return timepoint.strip()
        
        
def make_quartiles_columns(astro_df):
    
    pos_1, pos_2, pos_3 = 6, 7, 8
    astro_id, timepoint, flight, telo_data = 1, 2, 3, 4

    for i, row in astro_df.iterrows():
        
        astro_id_4digit = row[astro_id]
        
        if row[flight] == 'Pre-Flight' and row[timepoint] == 'L-270':
            preFlight_telos = row[telo_data]
            astro_df.iat[i, pos_1], astro_df.iat[i, pos_2], astro_df.iat[i, pos_3] = (quartile_cts_rel_to_df1(preFlight_telos, preFlight_telos))
            
        elif row[flight] == 'Pre-Flight' and row[timepoint] == 'L-180':
            if 'L-270' in list(astro_df[astro_df['astro id'] == astro_id_4digit]['timepoint']):
#                 print(f'L-270 present for {row[astro_id]}.. continuing')
                astro_df.iat[i, pos_1], astro_df.iat[i, pos_2], astro_df.iat[i, pos_3] = (quartile_cts_rel_to_df1(preFlight_telos, row[telo_data]))
                
            elif 'L-270' not in list(astro_df[astro_df['astro id'] == astro_id_4digit]['timepoint']):
#                 print(f'L-270 is NOT present for {row[astro_id]}.. assigning L180')

                preFlight_telos = row[telo_data]
                astro_df.iat[i, pos_1], astro_df.iat[i, pos_2], astro_df.iat[i, pos_3] = (quartile_cts_rel_to_df1(preFlight_telos, preFlight_telos))
            
        elif row[flight] == 'Pre-Flight':
            astro_df.iat[i, pos_1], astro_df.iat[i, pos_2], astro_df.iat[i, pos_3] = (quartile_cts_rel_to_df1(preFlight_telos, row[telo_data]))

        elif row[flight] == 'Mid-Flight':
            astro_df.iat[i, pos_1], astro_df.iat[i, pos_2], astro_df.iat[i, pos_3] = (quartile_cts_rel_to_df1(preFlight_telos, row[telo_data]))

        elif row[flight] == 'Post-Flight':
            astro_df.iat[i, pos_1], astro_df.iat[i, pos_2], astro_df.iat[i, pos_3] = (quartile_cts_rel_to_df1(preFlight_telos, row[telo_data]))

        else:
            print('unknown label in row[1] of the all patients df.. please check patient timepoint names')
    
    return astro_df



def graphing_statistics_telomere_data(dict_astro_individ_telos_dfs):  
    
    astro_list_of_IDs = ['5163', '2171', '1536', '7673', '4819', '3228', 
                         '2494', '2479', '2381', '1261', '1062']
    
    timepoint_series = ['L-270', 'L-180', 'L-60', 'FD45', 'FD90', 'FD140', 
                       'FD260', 'R+5', 'R+7', 'R+60', 'R+180', 'R+270']

    n=0
    
    for idNO in astro_list_of_IDs:
        
        n+=1
    #   #initialize blank list of timepoints
        data = [[1, 0, 0, 0], [0]]
        
        emptydata = pd.DataFrame(data)
        astro_L270 = pd.DataFrame(data)
        astro_L180 = pd.DataFrame(data)
        astro_L60 = pd.DataFrame(data)
        astro_Mid1 = pd.DataFrame(data)
        astro_Mid2 = pd.DataFrame(data)
        astro_R7 = pd.DataFrame(data)
        astro_R60 = pd.DataFrame(data)
        astro_R180 = pd.DataFrame(data)
        astro_R270 = pd.DataFrame(data)
        

        astro_L270name = ''
        astro_L180name = ''
        astro_L60name = ''
        astro_Mid1name = ''
        astro_Mid2name = ''
        astro_R7name = ''
        astro_R60name = ''
        astro_R180name = ''
        astro_R270name = ''

        for j in timepoint_series:
            for i in dict_astro_individ_telos_dfs.keys():
                
                
                if (idNO in i) and j == 'L-270' and ('L-270' in i):
                    astro_L270 = dict_astro_individ_telos_dfs[i]
                    astro_L270name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

                elif (idNO in i) and j == 'L-180' and ('L-180' in i):
                    astro_L180 = dict_astro_individ_telos_dfs[i]
                    astro_L180name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

                elif (idNO in i) and j == 'L-60' and ('L-60' in i):
                    astro_L60 = dict_astro_individ_telos_dfs[i]
                    astro_L60name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')
    
                elif (idNO in i) and (j == 'FD45' or j == 'FD90') and (j in i):
                    astro_Mid1 = dict_astro_individ_telos_dfs[i]                 
                    astro_Mid1name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')
                    
                elif (idNO in i) and (j == 'FD140' or j == 'FD260') and (j in i):
                    astro_Mid2 = dict_astro_individ_telos_dfs[i]
                    astro_Mid2name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

                elif (idNO in i) and j == 'R+7' and (j in i):
                    astro_R7 = dict_astro_individ_telos_dfs[i]               
                    astro_R7name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

                elif (idNO in i) and j == 'R+60' and (j in i):
                    astro_R60 = dict_astro_individ_telos_dfs[i]                      
                    astro_R60name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

                elif (idNO in i) and j == 'R+180' and (j in i):
                    astro_R180 = dict_astro_individ_telos_dfs[i]                 
                    astro_R180name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

                elif (idNO in i) and j == 'R+270' and (j in i):
                    astro_R270 = dict_astro_individ_telos_dfs[i]           
                    astro_R270name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')
                    
                else:
                    continue

        if idNO == '5163' or idNO == '2171' or idNO == '1536':

            if (astro_L270.size > 25 or astro_L180.size > 25) and (astro_Mid1.size > 25 and astro_Mid2.size > 25 ) and (astro_R180.size > 25 or astro_R270.size > 25):
                
                n_cells = 30
                
                astro_L270 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_L270, 'rsamp')
                astro_L180 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_L180, 'rsamp')
                astro_Mid1 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_Mid1, 'rsamp')
                astro_Mid2 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_Mid2, 'rsamp')
                astro_R180 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_R180, 'rsamp')
                astro_R270 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_R270, 'rsamp')
    
                n_bins = 30
                fig, axs = plt.subplots(2,2, sharey=True, tight_layout=False, figsize = (16, 12))

                if astro_L270name != '': 
                    if astro_R270name != '':
                        
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L270, astro_L270, astro_L270name, 0, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L270, astro_Mid1name, 0, 1)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L270, astro_Mid2name, 1, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R270, astro_L270, astro_R270name, 1, 1)
#                         print('stats')
#                         statistics_between_timepoints(astro_L270, astro_Mid1, astro_Mid2, astro_R270, 
#                                       astro_L270name, astro_Mid1name, astro_Mid2name, astro_R270name)

                    elif astro_R270name == '':
        
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L270, astro_L270, astro_L270name, 0, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L270, astro_Mid1name, 0, 1)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L270, astro_Mid2name, 1, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R180, astro_L270, astro_R180name, 1, 1)
#                         print('stats')
#                         statistics_between_timepoints(astro_L270, astro_Mid1, astro_Mid2, astro_R180, 
#                                       astro_L270name, astro_Mid1name, astro_Mid2name, astro_R180name)

                elif astro_L270name == '':
                    if astro_R270name == '':
                
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L180, astro_L180, astro_L180name, 0, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L180, astro_Mid1name, 0, 1)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L180, astro_Mid2name, 1, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R180, astro_L180, astro_R180name, 1, 1)
#                         print('randomly sampled stats')
#                         statistics_between_timepoints(astro_L180, astro_Mid1, astro_Mid2, astro_R180, 
#                                       astro_L180name, astro_Mid1name, astro_Mid2name, astro_R180name)

                    elif astro_R270name != '':
        

                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L180, astro_L180, astro_L180name, 0, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L180, astro_Mid1name, 0, 1)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L180, astro_Mid2name, 1, 0)
                        astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R270, astro_L180, astro_R270name, 1, 1)
#                         print('randomly sampled stats')
#                         statistics_between_timepoints(astro_L180, astro_Mid1, astro_Mid2, astro_R270, 
#                                       astro_L180name, astro_Mid1name, astro_Mid2name, astro_R270name)

                else:
                    continue

                # plt.savefig('Final telomere histogram random sampling dso'+idNO+'.pdf')
                plt.show()



        if idNO in ['7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']:
            if (astro_L270.size > 25) and (astro_R270.size > 25):
                
                n_cells = 30
                
#                 astro_L270name = f'synthetic astronaut {n} L+270'
#                 astro_R270name = f'synthetic astronaut {n} R+270'
                astro_L270 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_L270, 'rsamp')
                astro_R270 = gen_missing_values_andimpute_or_randomsampledown(n_cells, 184, astro_R270, 'rsamp')

                n_bins = 30
                fig, axs = plt.subplots(2, sharey=True, tight_layout=False, figsize = (12, 14))
                astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astro_L270, astro_L270, astro_L270name, 0)
                astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astro_R270, astro_L270, astro_R270name, 1)

#                 statistics_between_timepoints_prepost_only(astro_L270, astro_R270, astro_L270name, astro_R270name)

            else:
                continue

#             plt.savefig('Resampled telomere histogram dso'+idNO+'.pdf')
            plt.show()
            
        
        
        
def grab_control_values_generate_dictionary(patharg):

    """

    """

    dict_mean_individ_telos_dfs = {}

    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
            print(file.name, 'telomere data acquisition in progress..')
        
            try:
                df = pd.read_excel(file)

            except:
                print('File not found..')
                return -1

            df.rename(columns={'Unnamed: 3':'Mean Individ Telos'}, inplace=True)

            mean_values_of_individual_telomere_lengths = (df['Mean Individ Telos'])
            mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.drop(labels=[5, 192, 379, 566, 753, 940, 1127, 1314,
                    1501, 1688, 1875, 2062, 2249, 2436, 2623, 2810, 2997, 3184, 3371, 3558, 3745, 3932, 4119, 4306, 4493, 4680, 4867, 5054, 5241, 5428])
            mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.iloc[7:5611]
            meantelos_str_toNaN = pd.to_numeric(mean_values_of_individual_telomere_lengths, errors='coerce')
            mean_individual_telos_cleaned = meantelos_str_toNaN.dropna(axis=0, how='any')
            mean_individ_df = mean_individual_telos_cleaned.to_frame(name=None)
            mean_individ_df = mean_individ_df[(np.abs(stats.zscore(mean_individ_df)) < 3).all(axis=1)]

            if '0397' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(2.285)

            elif '3907' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(2.179)

            elif '1826' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(2.143)

            elif '0100' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(59.86)

            elif '0912' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(80.5)

            elif '0646' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(80.5)

            else:
                mean_individ_df_cy3Cal = mean_individ_df

            file_name_trimmed = file.name.replace('.xlsx', '')
            
            
            mean_individ_df_cy3Cal = gen_missing_values_andimpute_or_randomsampledown(30, 184, mean_individ_df_cy3Cal, 'rsamp')
            
            dict_mean_individ_telos_dfs[file_name_trimmed] = mean_individ_df_cy3Cal

    print('data collection complete')
    return dict_mean_individ_telos_dfs


def grab_control_telo_values_per_cell_generate_dictionary(patharg):

    """

    """

    dict_mean_individ_telos_dfs = {}

    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
            print(file.name, 'telomere data acquisition in progress..')
        
            try:
                df = pd.read_excel(file, skiprows=3)
                df = df.iloc[0:30, 12].to_frame()

            except:
                print('File not found..')
                return -1

            mean_individ_df = df.dropna(axis=0, how='any')

            if '0397' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(2.285)

            elif '3907' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(2.179)

            elif '1826' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(2.143)

            elif '0100' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(59.86)

            elif '0912' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(80.5)

            elif '0646' in file.name:
                mean_individ_df_cy3Cal = mean_individ_df.div(80.5)

            else:
                mean_individ_df_cy3Cal = mean_individ_df

            file_name_trimmed = file.name.replace('.xlsx', '')
            
            mean_individ_df_cy3Cal = mean_individ_df_cy3Cal.div(116.1848153)
            dict_mean_individ_telos_dfs[file_name_trimmed] = mean_individ_df_cy3Cal

    print('data collection complete')
    return dict_mean_individ_telos_dfs


def grab_astro_telo_values_per_cell_generate_dictionary(patharg):

    dict_astro_individ_telos_dfs = {}

    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
            print(f'{file.name} telomere data acquisition in progress..')
        
            try:
                df = pd.read_excel(file, skiprows=3)
                df = df.iloc[0:30, 12].to_frame()

            except:
                print(f'{file.name} File not found..')
                return -1
            
            telos_individ_df = df.dropna(axis=0, how='any')
            
            if ('5163' in file.name) or ('1536' in file.name):
                telos_individ_df_cy3Cal = telos_individ_df.div(59.86)

            elif '2171' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(80.5)

            elif '7673' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(2.11)

            elif '2479' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(2.18)

            elif '1261' in file.name:
                telos_individ_df_cy3Cal = telos_individ_df.div(2.16)

            else:
                telos_individ_df_cy3Cal = telos_individ_df
            
            telos_individ_df_cy3Cal = telos_individ_df_cy3Cal.div(116.1848153)

            file_name_trimmed = file.name.replace('.xlsx', '')
            dict_astro_individ_telos_dfs[file_name_trimmed] = telos_individ_df_cy3Cal

    print('Done collecting all astronaut telomere length excel files')
    return dict_astro_individ_telos_dfs


def raincloud_plot_astros_groups(x=None, y=None, data=None, 
                                 groupby=None, iterable=None):
    
    group_df = data.groupby(groupby)
    
    for item in iterable:
        plot_df = group_df.get_group(item)
        if x == 'timepoint':
        #this line only needed for timepoint
            plot_df[x].cat.remove_unused_categories(inplace=True)
    
        ax = sns.set(font_scale=1)
        #bw = sigma
        ax = pt.RainCloud(x = x, y = y, data = plot_df, palette = "Set2", bw = .20, 
                     width_viol = .8, figsize = (8,6), move=0.21, orient = "h")
        plt.title(f'{item} telos', fontsize=16)        
    
        
def make_astronaut_dataframe(dict_astro_individ_telos_dfs):
    data = []
    
    for name_key, telo_value in dict_astro_individ_telos_dfs.items():
        astro_id = name_key[3:7]
        astro_num, synth = get_astro_number_from_id(astro_id)
        time_point = get_timepoint(name_key)
        flight_status = relative_flight_timepoint(name_key)
        telo_value = gen_missing_values_andimpute_or_randomsampledown(30, 184, pd.Series(telo_value.values.reshape(-1,)), 'rsamp')

        data.append([astro_num, astro_id, time_point, flight_status, telo_value, np.mean(telo_value.values)])

    astro_df = pd.DataFrame(data, columns = ['astro number', 'astro id', 'timepoint', 'flight status', 'telo data', 'telo means'])

    sorter = ['L-270', 'L-180', 'L-60', 'FD45', 'FD90', 'FD140', 'FD260', 'R+5', 'R+7', 'R+60', 'R+105', 'R+180', 'R+270']
    astro_df['timepoint'] = astro_df['timepoint'].astype('category')
    astro_df['timepoint'].cat.set_categories(sorter, inplace=True)

    astro_df['Q1'] = 'telos preF Q1 <0.25'
    astro_df['Q2-3'] = 'telos preF Q2-3 >0.25 & <0.75'
    astro_df['Q4'] = 'telos preF Q4 >0.75'

    astro_df = astro_df.sort_values(['astro number', 'timepoint']).reset_index(drop=True)
    
    return astro_df


def make_astronaut_cell_data_dataframe(dict_astro_individ_telos_dfs):
    data = []
    
    for name_key, telo_value in dict_astro_individ_telos_dfs.items():
        astro_id = name_key[3:7]
        astro_num, synth = get_astro_number_from_id(astro_id)
        time_point = get_timepoint(name_key)
        flight_status = relative_flight_timepoint(name_key)
        telo_value = pd.Series(telo_value.values.reshape(-1,))

        data.append([astro_num, astro_id, time_point, flight_status, telo_value, np.mean(telo_value.values)])

    astro_df = pd.DataFrame(data, columns = ['astro number', 'astro id', 'timepoint', 'flight status', 'telo data per cell', 'telo means'])

    sorter = ['L-270', 'L-180', 'L-60', 'FD45', 'FD90', 'FD140', 'FD260', 'R+5', 'R+7', 'R+60', 'R+105', 'R+180', 'R+270']
    astro_df['timepoint'] = astro_df['timepoint'].astype('category')
    astro_df['timepoint'].cat.set_categories(sorter, inplace=True)

    astro_df['Q1'] = 'telos preF Q1 <0.25'
    astro_df['Q2-3'] = 'telos preF Q2-3 >0.25 & <0.75'
    astro_df['Q4'] = 'telos preF Q4 >0.75'

    astro_df = astro_df.sort_values(['astro number', 'timepoint']).reset_index(drop=True)
    
    return astro_df
        
    
def make_control_dataframe(dict_astro_individ_telos_dfs):
    data = []
    
    for name_key, telo_value in dict_astro_individ_telos_dfs.items():
        astro_id = name_key[3:7]
#         astro_num, synth = get_astro_number_from_id(astro_id)
        time_point = get_timepoint(name_key)
        flight_status = relative_flight_timepoint(name_key)
        telo_value = pd.Series(telo_value.values.reshape(-1,))

        data.append([astro_id, time_point, flight_status, telo_value, np.mean(telo_value.values)])

    astro_df = pd.DataFrame(data, columns = ['control id', 'timepoint', 'flight status controls', 'telo data', 'telo means'])

    sorter = ['L-270', 'L-180', 'L-60', 'FD45', 'FD90', 'FD140', 'FD260', 'R+5', 'R+7', 'R+60', 'R+105', 'R+180', 'R+270']
    astro_df['timepoint'] = astro_df['timepoint'].astype('category')
    astro_df['timepoint'].cat.set_categories(sorter, inplace=True)

    astro_df = astro_df.sort_values(['control id', 'timepoint']).reset_index(drop=True)
    
    return astro_df


def make_control_cell_data_dataframe(dict_astro_individ_telos_dfs):
    data = []
    
    for name_key, telo_value in dict_astro_individ_telos_dfs.items():
        astro_id = name_key[3:7]
#         astro_num, synth = get_astro_number_from_id(astro_id)
        time_point = get_timepoint(name_key)
        flight_status = relative_flight_timepoint(name_key)
        telo_value = pd.Series(telo_value.values.reshape(-1,))

        data.append([astro_id, time_point, flight_status, telo_value, np.mean(telo_value.values)])

    astro_df = pd.DataFrame(data, columns = ['control id', 'timepoint', 'flight status controls', 'telo data per cell', 'telo means'])

    sorter = ['L-270', 'L-180', 'L-60', 'FD45', 'FD90', 'FD140', 'FD260', 'R+5', 'R+7', 'R+60', 'R+105', 'R+180', 'R+270']
    astro_df['timepoint'] = astro_df['timepoint'].astype('category')
    astro_df['timepoint'].cat.set_categories(sorter, inplace=True)

    astro_df = astro_df.sort_values(['control id', 'timepoint']).reset_index(drop=True)
    
    return astro_df
        

def mid_split(row):
    if 'FD90' in row or 'FD45' in row:
        return 'Mid-Flight 1'
    elif 'FD140' in row or 'FD260' in row:
        return 'Mid-Flight 2'
    elif 'L' in row:
        return 'Pre-Flight'
    elif 'R' in row:
        return 'Post-Flight'
    

def histogram_plot_groups(x=None, data=None, 
                                 groupby=None, iterable=None):
    
    group_df = data.groupby(groupby)
    
    for item in iterable:
        plot_df = group_df.get_group(item)
        
        non_irrad = plot_df[plot_df['timepoint'] == '1 non irrad'][x]
        irrad_4_Gy = plot_df[plot_df['timepoint'] == '2 irrad @ 4 Gy'][x]
        three_B = plot_df[plot_df['timepoint'] == '3 B'][x]
        four_C = plot_df[plot_df['timepoint'] == '4 C'][x]

        n_bins = 70
        fig, axs = plt.subplots(2, 2, sharey=True, tight_layout=False, figsize=(20, 13))
        
        ax = sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
        ax = sns.set(font_scale=1)
        
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, non_irrad, non_irrad, f'patient #{item} 1 non rad', 0, 0)
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, irrad_4_Gy, non_irrad, f'patient #{item} 2 irrad @ 4 Gy', 0, 1)
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, three_B,  non_irrad, f'patient #{item} 3 B', 1, 0)
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, four_C,  non_irrad, f'patient #{item} 4 C', 1, 1)
        
            
        
            
def initialize_telo_data_1st_timepoint_variable(timepoint=None, df=None):

    if timepoint in list(df['timepoint'].unique()):
        variable = df[df['timepoint'] == str(timepoint)]['telo data exploded']
        return variable
    
    elif timepoint not in list(df['timepoint'].unique()):
        variable = pd.DataFrame([[0,1],[0,1]])
        return variable
    
    
def initialize_telo_data_timepoint_or_blank(timepoint, df):
    if timepoint in list(df['timepoint'].unique()):
        timepoint_telo_data = df[df['timepoint'] == str(timepoint)]['telo data exploded']
        
        name_id = str(df['astro id'].unique()[0])
        name_timepoint = f' {timepoint}'
        name_total = 'dso' + name_id + name_timepoint
        return name_total, timepoint_telo_data
        
    elif timepoint not in list(df['timepoint'].unique()):
        timepoint_telo_data = pd.DataFrame([0,1],[0,1])
        name = ''
        return name, timepoint_telo_data
    
########################################################################################################################
########################################################################################################################

# FUNCTIONS FOR GRAPHING INDIVIDUAL TELOMERES 

########################################################################################################################
########################################################################################################################
    
            
def graph_four_histograms(quartile_ref, n_bins, df1, df2, df3, df4,
                                                name1, name2, name3, name4):
    
    n_bins = n_bins
    fig, axs = plt.subplots(2,2, sharey=True, sharex=True, constrained_layout=True, figsize = (8, 6))
    sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
    
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    
    plt.rc('xtick',labelsize=10)
    plt.rc('ytick',labelsize=10)
    
#     csfont = {'fontname':'sans-serif'}
#     plt.suptitle(f"Individual Telomere Length Distributions at \nPre, Mid-1, Mid-2, and Post-Flight: {name1[0:8]}", 
#                  y=.95, fontsize=14, **csfont)

    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df1, quartile_ref, name1, 0, 0)
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df2, quartile_ref, name2, 0, 1)
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df3, quartile_ref, name3, 1, 0)
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df4, quartile_ref, name4, 1, 1)
    
    
    
def graph_two_histograms(quartile_ref, n_bins, df1, df2,
                                               name1, name2, controls=None):
    
    n_bins = n_bins
    fig, axs = plt.subplots(2, sharey=True, constrained_layout=True, figsize = (8, 6))
    sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
    
    for ax in axs.flat:
        ax.label_outer()
        
    plt.rc('xtick',labelsize=10)
    plt.rc('ytick',labelsize=10)
    
    astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, df1, quartile_ref, name1, 0)
    astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, df2, quartile_ref, name2, 1)
    
#     csfont = {'fontname':'sans-serif'}
#     plt.suptitle(f"Individual Telomere Length Distributions at Pre and Post-Flight: {name1[0:8]}", 
#                  y=.95, fontsize=14, **csfont)
    
#     if controls == True:
#         csfont = {'fontname':'sans-serif'}
#         plt.suptitle(f"Individual Telomere Length Distributions at Pre and Post-Flight: All Control Samples", 
#                      y=.95, fontsize=14, **csfont)

    
def astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone, axsNUMtwo):

    astroDF = astroDF.to_numpy()
    astroquartile = astroquartile.to_numpy()

    N, bins, patches = axs[axsNUMone,axsNUMtwo].hist(astroDF, bins=n_bins, range=(0, 400), edgecolor='black')

    for a in range(len(patches)):
        if bins[a] <= np.quantile(astroquartile, 0.25):
            patches[a].set_facecolor('#fdff38')

        elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
            patches[a].set_facecolor('#d0fefe')

        elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
            patches[a].set_facecolor('#d0fefe')

        elif bins[a] > np.quantile(astroquartile, 0.75): 
            patches[a].set_facecolor('#ffbacd')

    modified_astroname = astroname.replace('astro', '')
    axs[axsNUMone,axsNUMtwo].set_title(f"{modified_astroname}", fontsize=10,)
    
    font_axes=10

    if axsNUMone == 0 and axsNUMtwo == 0:
        axs[axsNUMone,axsNUMtwo].set_ylabel("Individual Telomere Counts", fontsize=font_axes)
        
    if axsNUMone == 1 and axsNUMtwo == 0:
        axs[axsNUMone,axsNUMtwo].set_ylabel("Individual Telomere Counts", fontsize=font_axes)
        axs[axsNUMone,axsNUMtwo].set_xlabel("Bins of Individual Telomeres (RFI)", fontsize=font_axes)
            
    if axsNUMone == 1 and axsNUMtwo == 1:
        axs[axsNUMone,axsNUMtwo].set_xlabel("Bins of Individual Telomeres (RFI)", fontsize=font_axes)
                  
    axs[axsNUMone,axsNUMtwo].xaxis.set_major_locator(plt.MaxNLocator(10))
        

        
def astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone):

    astroDF = astroDF.to_numpy()
    astroquartile = astroquartile.to_numpy()


    N, bins, patches = axs[axsNUMone].hist(astroDF, bins=n_bins, range=(0, 400), edgecolor='black')

    for a in range(len(patches)):
        if bins[a] <= np.quantile(astroquartile, 0.25):
            patches[a].set_facecolor('#fdff38')
        elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
            patches[a].set_facecolor('#d0fefe')
        elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
            patches[a].set_facecolor('#d0fefe')
        elif bins[a] > np.quantile(astroquartile, 0.75): 
            patches[a].set_facecolor('#ffbacd')
            
    axs[axsNUMone].set_title(f"{astroname}", fontsize=10,)
    
    font_axes=10

    if axsNUMone == 0 or axsNUMone == 1:
        axs[axsNUMone].set_ylabel("Individual Telomere Counts", fontsize=font_axes)
    if axsNUMone == 1:
        axs[axsNUMone].set_xlabel("Bins of Individual Telomeres (RFI)", fontsize=font_axes)
            
    axs[axsNUMone].xaxis.set_major_locator(plt.MaxNLocator(10))
    
    
    
def make_histograms_colored_by_quartile_for_astronauts(exploded_telos_df=None, astro_ids=None, nbins=45):

#     astro_ids = ['5163', '2171', '1536', '7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']
    
    grouped_data = exploded_telos_df.groupby('astro id')

    # by looping through astronaut ids, we'll pull out their respective dataframes
    # once we have the astronauts respective dfs, we'll figure out the quartile df & 
    for astro_id_num in astro_ids:
        
        if astro_id_num not in grouped_data.groups.keys():
            break
        plot_df = grouped_data.get_group(astro_id_num)
        for timepoint in ['L-270', 'L-180']:
            first_timepoint = initialize_telo_data_1st_timepoint_variable(timepoint=timepoint, df=plot_df)
            if first_timepoint.size > 30:
                break
        quartile_ref = first_timepoint

    #     okay, now we have the first timepoint as the reference for making quartile cutoffs! 
    #     now need to intialize other values!

        name_L270, astro_L270 = initialize_telo_data_timepoint_or_blank('L-270', plot_df)
        name_L180, astro_L180 = initialize_telo_data_timepoint_or_blank('L-180', plot_df)

        if '5163' == astro_id_num or '1536' == astro_id_num:
            name_Mid1, astro_Mid1 = initialize_telo_data_timepoint_or_blank('FD90', plot_df)
            name_Mid2, astro_Mid2 = initialize_telo_data_timepoint_or_blank('FD140', plot_df)
        if '2171' == astro_id_num:
            name_Mid1, astro_Mid1 = initialize_telo_data_timepoint_or_blank('FD45', plot_df)
            name_Mid2, astro_Mid2 = initialize_telo_data_timepoint_or_blank('FD260', plot_df)

        name_R180, astro_R180 = initialize_telo_data_timepoint_or_blank('R+180', plot_df)
        name_R270, astro_R270 = initialize_telo_data_timepoint_or_blank('R+270', plot_df)


        if ('5163' == astro_id_num) or ('2171' == astro_id_num) or ('1536' == astro_id_num):
            n_bins = n_bins

            if name_L270 != '': 
                        if name_R270 != '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L270, astro_Mid1, astro_Mid2, astro_R270,
                                                                    name_L270, name_Mid1, name_Mid2, name_R270)
                        elif name_R270 == '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L270, astro_Mid1, astro_Mid2, astro_R180,
                                                                    name_L270, name_Mid1, name_Mid2, name_R180)
            elif name_L270 == '':
                        if name_R270 != '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L180, astro_Mid1, astro_Mid2, astro_R270,
                                                                    name_L180, name_Mid1, name_Mid2, name_R270)
                        elif name_R270 == '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L180, astro_Mid1, astro_Mid2, astro_R180,
                                                                    name_L180, name_Mid1, name_Mid2, name_R180)

        elif astro_id_num in ['7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']:
            
            n_bins = 60
            graph_two_histograms(quartile_ref, n_bins, astro_L270, astro_R270, name_L270, name_R270)
            
        plt.savefig(f'../individual telomere length histogram distributions/png/dso{astro_id_num} histogram of individual telomere length distributions.png', dpi=600)
        
        plt.savefig(f'../individual telomere length histogram distributions/svg/dso{astro_id_num} histogram of individual telomere length distributions.svg', format='svg', dpi=1500)
    
    
def initialize_encoded_telo_data_timepoint_or_blank(timepoint, df):
    if timepoint in list(df['timepoint'].unique()):
        timepoint_telo_data = df[df['timepoint'] == str(timepoint)]['telo data exploded']
        
        name_id = str(df['encoded astro id'].unique()[0])
        name_timepoint = f' {timepoint}'
        name_total = 'astro ' + name_id + name_timepoint
        return name_total, timepoint_telo_data
        
    elif timepoint not in list(df['timepoint'].unique()):
        timepoint_telo_data = pd.DataFrame([0,1],[0,1])
        name = ''
        return name, timepoint_telo_data
    
    
def make_histograms_colored_by_quartile_for_encoded_astronauts(exploded_telos_df=None, astro_ids=None, n_bins=60):
    grouped_data = exploded_telos_df.groupby('encoded astro id')
    for astro_id_num in astro_ids:
        if astro_id_num not in grouped_data.groups.keys():
            break
        plot_df = grouped_data.get_group(astro_id_num)

        for timepoint in ['L-270', 'L-180']:
            first_timepoint = initialize_telo_data_1st_timepoint_variable(timepoint=timepoint, df=plot_df)
            if first_timepoint.size > 30:
                break
        quartile_ref = first_timepoint

        name_L270, astro_L270 = initialize_encoded_telo_data_timepoint_or_blank('L-270', plot_df)
        name_L180, astro_L180 = initialize_encoded_telo_data_timepoint_or_blank('L-180', plot_df)

        if 'B' == astro_id_num or 'C' == astro_id_num:
            name_Mid1, astro_Mid1 = initialize_encoded_telo_data_timepoint_or_blank('FD90', plot_df)
            name_Mid2, astro_Mid2 = initialize_encoded_telo_data_timepoint_or_blank('FD140', plot_df)
        if 'A' == astro_id_num:
            name_Mid1, astro_Mid1 = initialize_encoded_telo_data_timepoint_or_blank('FD45', plot_df)
            name_Mid2, astro_Mid2 = initialize_encoded_telo_data_timepoint_or_blank('FD260', plot_df)
        name_R180, astro_R180 = initialize_encoded_telo_data_timepoint_or_blank('R+180', plot_df)
        name_R270, astro_R270 = initialize_encoded_telo_data_timepoint_or_blank('R+270', plot_df)


        if ('B' == astro_id_num) or ('A' == astro_id_num) or ('C' == astro_id_num):
            n_bins = n_bins
            if name_L270 != '': 
                        if name_R270 != '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L270, astro_Mid1, astro_Mid2, astro_R270,
                                                                    name_L270, name_Mid1, name_Mid2, name_R270)
                        elif name_R270 == '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L270, astro_Mid1, astro_Mid2, astro_R180,
                                                                    name_L270, name_Mid1, name_Mid2, name_R180)
            elif name_L270 == '':
                        if name_R270 != '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L180, astro_Mid1, astro_Mid2, astro_R270,
                                                                    name_L180, name_Mid1, name_Mid2, name_R270)
                        elif name_R270 == '':
                            graph_four_histograms(quartile_ref, n_bins, astro_L180, astro_Mid1, astro_Mid2, astro_R180,
                                                                    name_L180, name_Mid1, name_Mid2, name_R180)

#         elif astro_id_num in ['7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']:
#             n_bins = 60
#             graph_two_histograms(quartile_ref, n_bins, astro_L270, astro_R270, name_L270, name_R270)
            
        plt.savefig(f'../individual telomere length histogram distributions/png/dso{astro_id_num} histogram of individual telomere length distributions.png', dpi=600)
        
        plt.savefig(f'../individual telomere length histogram distributions/svg/dso{astro_id_num} histogram of individual telomere length distributions.svg', format='svg', dpi=1500)


########################################################################################################################
########################################################################################################################

# FUNCTIONS FOR CORRELATING TELOMERES WITH ANALYTE DATA

########################################################################################################################
########################################################################################################################

def select_astros_of_interest(analyte_df, telomere_df, astro_ids_of_interest, target):
    
    telomere_df['astro id'] = telomere_df['astro id'].astype('str')
    
    if 'astro id' in analyte_df.columns:
        analyte_df['astro id'] = analyte_df['astro id'].astype('str')
    
    if 'sample type' in analyte_df.columns:
        analyte_df.drop('sample type', axis=1, inplace=True)
    
    # dropping unnecessary cols from telo df
    trim_astro_df = telomere_df.drop(['astro number', 'timepoint'], axis=1)
    
    if astro_ids_of_interest == 'all astros':
        
        # i.e as of 10/7/19 I only have n=4 (contains astro id col) & n=11 (no astro id) dataframes for analytes
        # I think when I received n=3 astros.. just type astro ids for astro_ids_of_interest, it will work properly
        # or.. if i receive n=11 dataframe with labeled astros.. 
        # just rewrite this area to accept n=11 df w/ astro id col
        if 'astro id' in analyte_df.columns:
            (print("Possible error.. the astro id column is present.. all astros were requested but this df potentially" 
                  "contains less than all 11 astros.. drop astro id col and retry"))
            return
        else: 
            # retain all astro ids
            selected_astros = trim_astro_df
            id_values = ['flight status']

    elif astro_ids_of_interest != 'all astros':
        # subset astro ids of interest 
        selected_astros = trim_astro_df[trim_astro_df['astro id'].isin(astro_ids_of_interest)].reset_index(drop=True)
        id_values = ['astro id', 'flight status']
        
    return analyte_df, selected_astros, id_values


def merge_analyte_telomere_data(analyte_df, selected_astros, id_values, telos_percent_change, target):
    
    # take mean telomere length values of all astronauts or per astros of interest & merge with analytes 
    mean_selected_astros = selected_astros.groupby(id_values).agg('mean').reset_index()
    
    if telos_percent_change == 'yes':
        mean_selected_astros[target] = (mean_selected_astros[target]
                                                      .apply(lambda row: make_telos_percent_change(row)))
        
    merge_analyte_df = analyte_df.merge(mean_selected_astros, on=id_values)
    
    # prepare to drop any columns w/ missing data
    indexer=['timepoint', target]
    for id_value in id_values:
        indexer.append(id_value)
        
    return merge_analyte_df, indexer


def how_drop_missing_values(merge_analyte_df, how_drop_missing, indexer):
    # drop every analyte (columns) with missing data
    if how_drop_missing == 'by column':
        pivot_merge = (merge_analyte_df.pivot_table(index=indexer, columns='biochemistry analyte', 
                                                    values='measured analyte').reset_index())
        pivot_merge.dropna(axis=1, inplace=True)
        cleaned_data = pivot_merge.melt(id_vars=indexer, var_name='biochemistry analyte', 
                                        value_name='measured analyte').reset_index(drop=True)
    
    # drop missing data on per analyte/timepoint/astro (row) basis 
    elif how_drop_missing == 'by melted row':
        cleaned_data = merge_analyte_df.dropna(axis=0)
        
    return cleaned_data


def retain_flight_status(cleaned_data, retain_what_flight_status):
    # retaining analytes for which flight status
    if retain_what_flight_status == 'any':
        retained_data = cleaned_data
        
    elif bool(set(retain_what_flight_status) & set(['Pre-Flight', 'Mid-Flight', 'Post-Flight'])) == True:
        retained_data = cleaned_data[cleaned_data['flight status'].isin(retain_what_flight_status)].copy()
        
    elif retain_what_flight_status == 'require at least one per status':
        total_analytes = list(cleaned_data['biochemistry analyte'].unique())
        analytes_3_unique_flight = []
        groupby_analyte = cleaned_data.groupby('biochemistry analyte')
        
        for analyte in total_analytes:
            # make groups by analyte
            get_group_by_analyte = groupby_analyte.get_group(analyte)
            
            # look at unique flight status values per analyte
            g_f_s_t = list(get_group_by_analyte['flight status'].unique())
            
            # if pre, mid, and post flight values in unique value list per analyte, then add this analyte to a list
            if 'Pre-Flight' in g_f_s_t and 'Mid-Flight' in g_f_s_t and 'Post-Flight' in g_f_s_t:
                analytes_3_unique_flight.append(analyte)

        # retain only analytes with at least one measurement per flight status 
        analytes_only_3_unique_df = cleaned_data[cleaned_data['biochemistry analyte'].isin(analytes_3_unique_flight)].copy()

        
        return analytes_only_3_unique_df
        
    return retained_data


def make_telos_percent_change(row):
    percent_chg_telos = ((row - 0.938117) / 0.938117) * 100
    return percent_chg_telos


def correlate_astro_analytes_telomeres_pipeline(analyte_df=None, telomere_df=None, target=None,
                                                astro_ids_of_interest=None,
                                                how_drop_missing=None, retain_what_flight_status=None,
                                                telos_percent_change='no'):
    """
    High level fxn description 
    
    Args:
        analyte_df (pandas dataframe): Contains either n=4 or n=11 biochemical analyte data in tidy data format.
        
        telomere_df (pandas dataframe): Must contain complete telomere length data in tidy data format.
        
        astro_ids_of_interest (str or list of str): Accepts either 'all astros' as str, whereby all astronaut data is 
        used for correlating telo/analyte data, or a list of astro ids to subset data for analysis.
        
        how_drop_missing (str): Accepts either 'by column', which drops any analyte containing at least one missing value,
        or 'by melted row', which drops only single instances of missing values.
        
        retain_what_flight_status (str or list of tring): decides how to subset individual analytes based on what 
        flight status labels they have 
        
            Accepts: 'any', whereby no subselection is placed on analytes based on flight status, 
            or: subset data by flight status (list of str) for all analytes as a GROUP i.e ['Pre-Flight'] or ['Pre-Flight', 'Post-Flight']
            or: 'require at least one per status', where EACH analytes must have at least one measurement per flight status 

    Returns:
        retained_data (pandas dataframe): Data subject to the processing steps described above.
    """
    
    # selecting astros of interest & capturing id values for handling merges 
    analyte_df, selected_astros, id_values = select_astros_of_interest(analyte_df, telomere_df, astro_ids_of_interest, target)

    # merging analyte & telomere data, capturing indexer for handling missing data
    merge_analyte_df, indexer = merge_analyte_telomere_data(analyte_df, selected_astros, id_values, telos_percent_change, target)
    
    # dropping missing values based on input
    cleaned_data = how_drop_missing_values(merge_analyte_df, how_drop_missing, indexer)

    # subsetting values based on flight status labels 
    retained_data = retain_flight_status(cleaned_data, retain_what_flight_status)
    
    return retained_data


def find_high_correlates_analytes_mean_telos(merged_analyte_blood_tidy_df, corr_cutoff, corr_loc=0, astro_ids=False, target=None):
    
    if astro_ids == False:
        corr_value_tests = []

        grouped_by_analyte = merged_analyte_blood_tidy_df.groupby('biochemistry analyte')

        for group in list(merged_analyte_blood_tidy_df['biochemistry analyte'].unique()):
            corr_value = grouped_by_analyte.get_group(group).corr()[target][corr_loc]

            if abs(corr_value) > corr_cutoff:
                corr_value_tests.append([group, corr_value])
#                 print(f"{group}: {corr_value:.4f}")

        return corr_value_tests
    
    elif astro_ids == True:
        
        corr_value_requested = input('Please state index for correlation value in corr().. 0 or 1')
        corr_value_tests = []
        astro_ids = list(merged_analyte_blood_tidy_df['astro id'].unique())
        astro_id_group = merged_analyte_blood_tidy_df.groupby('astro id')

        for astro in astro_ids:
            individ_astro_df = astro_id_group.get_group(astro)
            analyte_grouped_by_individ = individ_astro_df.groupby('biochemistry analyte')
            analytes = list(individ_astro_df['biochemistry analyte'].unique())

            for analyte in analytes:
                corr_value = analyte_grouped_by_individ.get_group(analyte).corr()[target][int(corr_value_requested)]
                corr_value_tests.append([astro, analyte, corr_value])
                
#                 if abs(corr_value) > corr_cutoff:
#                     print(f"{astro} - {analyte}: {corr_value:.4f}")
                    
        return corr_value_tests
    
    
def plot_diverging_correlations(list_correlates=None, target_name=None, figsize=(11,7), dpi=300):
    df = list_correlates
    x = df['correlation value']
    df['colors'] = ['black' if x < 0 else 'green' for x in df['correlation value']]
    df.sort_values('correlation value', inplace=True)
    df.reset_index(inplace=True, drop=True)

    plt.figure(figsize=figsize, dpi=dpi)
    plt.hlines(y=df.index, xmin=0, xmax=df['correlation value'], color=df['colors'], alpha=0.6, linewidth=7)

    # Decorations
    plt.yticks(df.index, df['biochemistry analyte'], fontsize=12)
    plt.xticks(fontsize=14)
    plt.xlabel(target_name, fontsize=16)
    plt.ylabel('Blood biochemistry analytes', fontsize=16)
#     plt.title(f'Correlation between {target_name} and Analytes', x=0.4, fontdict={'size':18})

    plt.grid(linestyle='-', alpha=.2, color='black')
    plt.tight_layout()
    my_xticks = np.array([-1, -.5, 0, .5, 1])
    plt.xticks(my_xticks[::1])
#     plt.savefig(f'diverging bars {target_name} n=11.png')
    
    
def scipy_anova_post_hoc_tests(df=None, flight_status_col='flight status', target='telo data per cell',
                               sig_test=stats.f_oneway, post_hoc=sp.posthoc_ttest):

    g_1 = df[df[flight_status_col] == 'Pre-Flight'][target]
    g_2 = df[df[flight_status_col] == 'Mid-Flight'][target]
    g_3 = df[df[flight_status_col] == 'Post-Flight'][target]
    statistic, p_value = sig_test(g_1, g_2, g_3)
    print(f'ONE WAY ANOVA for telomere length: {p_value}')
        
    # if anova detects sig diff, perform post-hoc tests
    if p_value <= 0.05:
        print('bonferroni')
        display(sp.posthoc_ttest(df, val_col=target, group_col=flight_status_col, equal_var=True,
                                 p_adjust=None))
        
        
def telos_scipy_anova_post_hoc_tests(df0=None, time_col='flight status', target='individual telomeres',
                                     sig_test=stats.f_oneway, post_hoc=None, repeated_measures=False):
    df = df0.copy()
    df.rename({'telo data per cell': 'telo_data_per_cell',
               'flight status': 'flight_status',
               'Mean Telomere Length (qPCR)': 'Mean_Telomere_Length_(qPCR)',
               'astro id': 'astro_id'}, axis=1, inplace=True)
              
    if ' ' in time_col:
        time_col = time_col.replace(' ', '_')
    if ' ' in target:
        target = target.replace(' ', '_')
    
    if repeated_measures == False:
        g_1 = df[df[time_col] == 'Pre-Flight'][target]
        g_2 = df[df[time_col] == 'Mid-Flight'][target]
        g_3 = df[df[time_col] == 'Post-Flight'][target]
        statistic, p_value = sig_test(g_1, g_2, g_3)
        print(f'ONE WAY ANOVA for telomere length: {p_value}')
              
    elif repeated_measures:
        results = AnovaRM(df, target, 'astro_id', 
                          within=[time_col], aggregate_func='mean').fit()
        # pvalue
        p_value = results.anova_table['Pr > F'][0]
        print(f'REPEATED MEASURES ANOVA for telomere length: {p_value}')     
          
    # if anova detects sig diff, perform post-hoc tests
    if p_value <= 0.05:
        mc = MultiComparison(df[target], df[time_col])
        mc_results = mc.tukeyhsd()
        print(mc_results)
        res = mc_results
        print(f'TukeyHSD pvalues: {list(psturng(np.abs(res.meandiffs / res.std_pairs), len(res.groupsunique), res.df_total))}')
        
#         print('\nbonferroni pvalues')
#         display(sp.posthoc_ttest(df, val_col=target, group_col=time_col, equal_var=False,
#                                  p_adjust='bonferroni'))
        

def id_encode_letters(row):
    if row == '2171':
        row = 'A'
    elif row == '5163':
        row = 'B'
    elif row == '1536':
        row = 'C'
    return row


def eval_make_test_comparisons(df=None, timepoints=None, test=None, test_name=None, 
                               target='individual telos'):
              
    timepoints = list(df['timepoint'].unique())
    timept_pairs = []
    row = []
    df_list = []
              
    for timept in timepoints:
        df_list.append(df[df['timepoint'] == timept][target])
              
    for iter1, df in zip(timepoints, df_list):
        for iter2, i in zip(timepoints, range(len(df_list))):
            pair1, pair2 = f"{iter1}:{iter2}", f"{iter2}:{iter1}"
            if iter1 != iter2 and pair1 not in timept_pairs and pair2 not in timept_pairs:
                stat, pvalue = test(df, df_list[i])
                print(f'{test_name} | {iter1} vs {iter2} {pvalue}')
                timept_pairs.append(pair1)
                timept_pairs.append(pair2)
                row.append([test_name, iter1, iter2, pvalue])
    return timept_pairs, row


def make_post_flight_df_and_merge(astro_df=None, exploded_telos=None, timepoint=None):
    """
    parse out mean telomere length & #s short/long telomeres from specific post-flight (R+7, R+60, ... R+270) timepoints
    and merge with exploded_telos dataframe for machine learning prep 
    """
    # parsing out post-flight data of interest
    timepoint_df = astro_df[astro_df['timepoint'] == timepoint].copy()
    for col in ['telo means', 'Q1', 'Q4']:
        timepoint_df.rename({col: f'{timepoint} {col}'}, axis=1, inplace=True)
        
    timepoint_df.drop(['astro number', 'timepoint', 'flight status'], axis=1, inplace=True)
    
    # extracting pre-flight individual telomere data only
    exploded_telos_pref = exploded_telos[exploded_telos['flight status'] == 'Pre-Flight'].copy()
    exploded_telos_pref.drop(['astro number', 'flight status'], axis=1, inplace=True)
    
    merge_df = exploded_telos_pref.merge(timepoint_df, on=['astro id'])
    return merge_df


class make_features(BaseEstimator, TransformerMixin):
    def __init__(self, make_log_individ_telos=False, make_log_target=False):
        self.make_log_individ_telos = make_log_individ_telos
        self.make_log_target = make_log_target
        
        
    def fit(self, X, y=None):
        return self
    
    
    def create_log_individ_telos(self, X, y=None):
        X['individual telos'] = np.log1p(X['individual telos'])
        return X
    
    
    def create_log_target(self, X, y=None):
        X['4 C telo means'] = np.log1p(X['4 C telo means'])
        return X
        
        
    def transform(self, X, y=None):
        if self.make_log_individ_telos:
            X = self.create_log_individ_telos(X)
            
        if self.make_log_target:
            X = self.create_log_target(X)
        return X
    
    
class make_dummies(BaseEstimator, TransformerMixin):
    def __init__(self, drop_first=True, cols_to_dummify=['timepoint'], how_dummify='encode'):
        self.drop_first = drop_first
        self.cols_to_dummify = cols_to_dummify
        self.how_dummify=how_dummify
        
    
    def fit(self, X, y=None):
        return self
    
    
    def transf_dummies(self, X, y=None):
        dummies = pd.get_dummies(X, drop_first=self.drop_first, columns=self.cols_to_dummify)
        return dummies
    
    def label_encode(self, X, y=None):
        label_encoder = preprocessing.LabelEncoder()
        X['encoded_timepoint'] = label_encoder.fit_transform(X[self.cols_to_dummify].values.ravel())
        X.drop(['timepoint'], axis=1, inplace=True)
        return X
    
    
    def transform(self, X, y=None):
        if self.how_dummify == 'get_dummies':
            X = self.transf_dummies(X)
        elif self.how_dummify == 'encode':
            X = self.label_encode(X)
        return X
    
    
class clean_data(BaseEstimator, TransformerMixin):
    def __init__(self, drop_astro_id=True, timepoint='R+7', target='telo means'):
        self.drop_astro_id = drop_astro_id
        self.timepoint_target = f'{timepoint} {target}'
        self.timepoint = timepoint
        self.target = target
    
    
    def fit(self, X, y=None):
        return self
    
    
    def transform(self, X, y=None):       

#         enforcing col types
        cols = list(X.columns)
        for col in cols:
            if 'individual telomeres' in col or 'telo means' in col:
                X[col] = X[col].astype('float64')
            else:
                X[col] = X[col].astype('int64')
        
        if self.drop_astro_id:
            X.drop(['astro id'], axis=1, inplace=True)
            
        X.reset_index(drop=True, inplace=True)
        
        target_cols = ['telo means', 'Q1', 'Q4']
        target_cols.remove(self.target)
        
        for item in target_cols:
            for col in X.columns:
                if item in col:
                    X.drop([col], axis=1, inplace=True)
        
#         if 'telo means' in self.target:
#             X.drop([f'{timepoint} Q1', f'{timepoint} Q4'], axis=1, inplace=True)
#         elif 'Q1' in self.target:
#             X.drop([f'{timepoint} Q1', f'{timepoint} Q4'], axis=1, inplace=True)
             
#         X = X[['encoded_timepoint', 'individual telomeres', self.timepoint_target]].copy()
        return X
    
    
def cv_score_fit_mae_test(train_set=None, test_set=None, target=None,
                          model=None, cv=5, scoring='neg_mean_absolute_error', verbose=True):
    random.seed(888)
    row = []
    features = [col for col in train_set.columns if col != target and col != 'astro id']
    
    X_train = train_set[features].copy()
    X_test = test_set[features].copy()
    
    y_train = train_set[target].copy()
    y_test = test_set[target].copy()
    
    # cv
    scores = -1 * cross_val_score(model, X_train, y_train, cv=5, scoring=scoring)
    if verbose:
        print(f'MAE per CV fold: \n{scores} \n')
        print(f'MEAN of MAE all folds: {scores.mean()}')
        print(f'STD of MAE all folds: {scores.std()}\n')

    # fitting the model
    model.fit(X_train, y_train)

    # predict y_test from X_test - this is using the train/test split w/o shuff;ing
    predict_y_test = model.predict(X_test)
    if verbose:
        print(f"MAE of predict_y_test & y_test: {mean_absolute_error(y_test, predict_y_test)}")
        print(f'R2 between predict_y_test & y_test: {r2_score(y_test, predict_y_test)}')
    
    row.append(['XGBoost', features, target, round(scores.mean(), 4),
                                             round(scores.std(), 4),
                                             round(mean_absolute_error(y_test, predict_y_test), 4), 
                                             round(r2_score(y_test, predict_y_test), 4)])
    
    return model, row 


def myMetric(x, y):
    r = stats.pearsonr(x, y)[0]
    return 1 - r 


def plot_dendogram(Z, target=None, indexer=None):
    with plt.style.context('fivethirtyeight' ): 
        plt.figure(figsize=(10, 2.5))
        plt.title(f'Dendrogram of clusters by {target}', fontsize=22, fontweight='bold')
        plt.xlabel('astro IDs', fontsize=22, fontweight='bold')
        plt.ylabel('distance', fontsize=22, fontweight='bold')
        hac.dendrogram(Z, labels=indexer, leaf_rotation=90.,    # rotates the x axis labels
                        leaf_font_size=15., ) # font size for the x axis labels
        plt.show()


def plot_results2(timeSeries, D, cut_off_level, y_size, x_size, verbose, time, target):
    result = pd.Series(hac.fcluster(D, cut_off_level, criterion='maxclust'))
    
    if verbose:
        clusters = result.unique() 
        fig = plt.subplots(figsize=(x_size, y_size))   
        mimg = math.ceil(cut_off_level/2.0)
        gs = gridspec.GridSpec(mimg,2, width_ratios=[1,1])
        cluster_indexed = pd.concat([result, timeSeries.reset_index(drop=True)], axis=1)
        
        columns = list(cluster_indexed.columns[1:])
        columns = ['clusters'] + columns
        cluster_indexed.columns = columns
        
        for ipic, c in enumerate(clusters):
            clustered = cluster_indexed[cluster_indexed['clusters'] == c].copy()
            print(ipic, "Cluster number %d has %d elements" % (c, len(clustered['astro id'])))
            melt = clustered.melt(id_vars=['astro id', 'clusters'], var_name=time,value_name=target)
            ax1 = plt.subplot(gs[ipic])
            sns.lineplot(x=time, y=target, hue='astro id', data=melt, legend=False, ax=ax1)
            ax1.set_title((f'Cluster number {c}'), fontsize=15, fontweight='bold')
        plt.tight_layout()
        
    return result
        
        
def cluster_data_return_df(df, target='inversions', time='timepoint', cut_off_n=4, 
                           metric=myMetric, method='single',
                           y_size=6, x_size=10, verbose=True):

    df = df.copy()
    
    label_enc = LabelEncoder()
    labels = list(df[time])
    encoded_labels = list(LabelEncoder().fit_transform(df[time]))
    cypher_dict = dict(zip(encoded_labels, labels))
    df[time] = LabelEncoder().fit_transform(df[time])
    
    df = df.pivot(index='astro id', columns=time, values=target).reset_index()
    
    # run the clustering    
    cluster_Z = hac.linkage(df, method=method, metric=metric)
    if verbose:
        plot_dendogram(cluster_Z, target=target, indexer=df.index)
    # return df bearing cluster groups
    indexed_clusters = plot_results2(df, cluster_Z, cut_off_n, y_size, x_size, verbose, time, target)
    
    # concat clusters to original df and return
    ready_concat = df.reset_index(drop=True)
    clustered_index_df = pd.concat([ready_concat, indexed_clusters], axis=1)
    clustered_index_df.columns = list(clustered_index_df.columns[:-1]) + [f'{target} cluster groups']

    melted = clustered_index_df.melt(id_vars=['astro id', f'{target} cluster groups'], var_name=time, value_name=target)
    melted[time] = melted[time].apply(lambda row: cypher_dict[row])
    
    return melted


def fish_assign_clustering(row):
    cluster_dict = {'5163': 1,
                    '2381': 1,
                    '2494': 1,
                    '1261': 2,
                    '1536': 2,
                    '7673': 2,
                    '2171': 2,
                    '4819': 2,
                    '3228': 3,
                    '1062': 3,
                    '2479': 3}
    return cluster_dict[row]


def qpcr_assign_cluster(row):
    cluster_grp_dict = {'2381': 1,
                        '4819': 2,
                        '5163': 2,
                        '2171': 2,
                        '3228': 2,
                        '1062': 2,
                        '2494': 2,
                        '7673': 2,
                        '2479': 3,
                        '1261': 3,
                        '1536': 3}
    return cluster_grp_dict[row]


def graph_cluster_groups(df, time=None, target=None, hue=None, colors='Set1', 
                         n_cols=3, y_label_name=None, figsize=(7,3.2),
                         fontsize=14, save=True, bbox_to_anchor=(0.5, 1.18),
                         y_lim=None, path_labels='11 astros'):
    
    colors = sns.color_palette(colors)
    
    plt.figure(figsize=figsize)
    ax = sns.lineplot(x=time, y=target, data=df, hue=hue, markers=True,
                      palette=sns.color_palette(colors[:len(df[hue].unique())]),
                      style=hue, **{'markersize':11, 'mec':'black', 'mew':1})

    plt.setp(ax.get_xticklabels(), 
#              rotation=45, 
             fontsize=fontsize)
    ax.set_ylabel(f'{y_label_name}', fontsize=fontsize)
        
    ax.set_xlabel('', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    
    legend = ax.legend()
    legend.texts[0].set_text('Cluster groups')
    
    if y_lim != None:
        ax.set_ylim(y_lim)
                         
    plt.legend(loc='upper center', bbox_to_anchor=bbox_to_anchor,
          ncol=n_cols, fancybox=True, fontsize=fontsize)
    if save:
        plt.savefig(f'../MANUSCRIPT 11 ASTROS/figures/{path_labels} lineplot {target} clustering.png', 
                dpi=600, bbox_inches = "tight")


def convert_mid_timepoint(row):
    if row == 'FD45' or row == 'FD90':
        return 'Mid-1'
    elif row == 'FD140' or row == 'FD260':
        return 'Mid-2'
    else:
        return row
    
    
def set_categories_sort(telomere_df=None, time='timepoint', sort_list=None):
    df = telomere_df.copy()
    if sort_list == None:
        sort_list = ['L-270', 'L-180', 'L-60', 'R+7', 'R+60', 'R+180', 'R+270']
        
    df[time] = df[time].astype('category')
    df[time].cat.set_categories(sort_list, inplace=True)
    return df


def ext_telo_data_longitudinal_clustering(telomere_df=None, 
                                          telomere_col_name='telo means',
                                          col_to_pivot='timepoint',
                                          timepts_of_interest=None):
    df = telomere_df.copy()
    if timepts_of_interest == None:
        timepts_of_interest = ['L-270', 'L-180', 'L-60', 'R+7', 'R+60', 'R+180', 'R+270']
    
    # parse cols of interest
    parsed_df = df[['astro id', col_to_pivot, telomere_col_name]].copy()
    parsed_df[col_to_pivot] =  parsed_df[col_to_pivot].astype('str')
    
    # pivot out timepoints
    pivot_df = parsed_df.pivot_table(index=['astro id'], columns=col_to_pivot, values=telomere_col_name).reset_index()
    pivot_df.set_index('astro id', inplace=True)
    cluster_ready_df = pivot_df[timepts_of_interest].copy()
    return cluster_ready_df


def rename_imputed_df(imputed_df=None, original_df=None):
    imputed_df.columns = original_df.columns
    imputed_df.index = original_df.index
    imputed_df.columns.name = ''
    return imputed_df


def clustermap_plot(df=None, method='single', metric='correlation', 
                    color_map='PRGn', col_cluster=False, fontsize=14, z_score=0,
                    y_label='Mean Telomere Length (Telo-FISH)', path_labels='11 astros',
                    save=True):

    g = sns.clustermap(df, method=method, metric=metric, z_score=z_score, figsize=(7,7), 
                       cmap=color_map, col_cluster=col_cluster) 

    # colorbar 
    g.cax.set_position([-0.05, .2, .03, .45])
    g.cax.set_ylabel(y_label, rotation=90, fontsize=fontsize)
    g.cax.tick_params(labelsize=12)

    # modifying y axis
    g.ax_heatmap.set_ylabel('Astronaut ID', fontsize=fontsize)
    g.ax_heatmap.set_xlabel('')
    labels = g.ax_heatmap.yaxis.get_majorticklabels()
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=fontsize)
    plt.setp(g.ax_heatmap.yaxis.get_minorticklabels(), fontsize=fontsize)
    g.ax_heatmap.set_yticklabels(labels, rotation=0, fontsize=fontsize, va="center")

    # modifying x axis
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, fontsize=fontsize)

    for a in g.ax_row_dendrogram.collections:
        a.set_linewidth(1)
    for a in g.ax_col_dendrogram.collections:
        a.set_linewidth(1)
    
    if save:
        plt.savefig(f'../MANUSCRIPT 11 ASTROS/figures/{path_labels} {y_label} cluster map.png', dpi=600, bbox_inches = "tight")
        
        
def flight_status(row):
    if 'FD90' in row or 'FD45' in row:
        return 'Mid-Flight'
    elif 'FD140' in row or 'FD260' in row:
        return 'Mid-Flight'
    elif 'L' in row:
        return 'Pre-Flight'
    elif 'R' in row:
        return 'Post-Flight'
    
    
# def encode_timepts(row):
#     encode_dict = {'L-270' : 1,
#                    'R+7': 2, 
#                    'R+270': 3}
#     return encode_dict[row]


# def myMetric(x, y):
#     r = stats.pearsonr(x, y)[0]
#     return 1 - r 


# def plot_dendogram(Z, target=None, indexer=None):
#     with plt.style.context('fivethirtyeight' ): 
#         plt.figure(figsize=(10, 2.5))
#         plt.title(f'Dendrogram of clusters by {target}', fontsize=22, fontweight='bold')
#         plt.xlabel('astro IDs', fontsize=22, fontweight='bold')
#         plt.ylabel('distance', fontsize=22, fontweight='bold')
#         hac.dendrogram(Z, labels=indexer, leaf_rotation=90.,    # rotates the x axis labels
#                         leaf_font_size=15., ) # font size for the x axis labels
#         plt.show()

        
# def plot_results(timeSeries, D, cut_off_level, y_size, x_size, verbose):
#     result = pd.Series(hac.fcluster(D, cut_off_level, criterion='maxclust'))
#     if verbose:
#         clusters = result.unique() 
#         fig = plt.subplots(figsize=(x_size, y_size))   
#         mimg = math.ceil(cut_off_level/2.0)
#         gs = gridspec.GridSpec(mimg,2, width_ratios=[1,1])
#         cluster_indexed = pd.concat([result, timeSeries.reset_index()], axis=1)
#         cluster_indexed.rename({0: 'clusters'}, axis=1, inplace=True)
        
#         for ipic, c in enumerate(clusters):
#             clustered = cluster_indexed[cluster_indexed['clusters'] == c].copy()
#             print(ipic, "Cluster number %d has %d elements" % (c, len(clustered['astro id'])))
#             clustered.drop(['index'], axis=1, inplace=True)
#             melt = clustered.melt(id_vars=['astro id', 'clusters'], var_name='timepoint',value_name='telo means')
#             melt = set_categories_sort(telomere_df=melt, sort_list=['L-270', 'R+7', 'R+270'])
#             ax1 = plt.subplot(gs[ipic])
#             melt
#             sns.lineplot(x='timepoint', y='telo means', hue='astro id', data=melt, legend=False, ax=ax1)
#             ax1.set_title((f'Cluster number {c}'), fontsize=15, fontweight='bold')
#         plt.tight_layout()
        
#     return result
        
        
# def cluster_telomere_data_return_df(df=None, target='telo means', cut_off_n=4, 
#                                     metric=myMetric, method='single',
#                                     y_size=6, x_size=10, verbose=True):
    
# #     astro_ids = df.index
# #     knn_telo_qpcr.reset_index(drop=True, inplace=True)
    
#     # run the clustering    
#     cluster_Z = hac.linkage(df, method=method, metric=metric)
    
#     if verbose:
#         plot_dendogram(cluster_Z, target=target, indexer=df.index)
#     # return df bearing cluster groups
#     df0 = df.copy().reset_index()
#     indexed_clusters = plot_results(df0, cluster_Z, cut_off_n, y_size=y_size, x_size=x_size, verbose=verbose)
    
#     # concat clusters to original df and return
#     ready_concat = df.reset_index()
#     clustered_index_df = pd.concat([ready_concat, indexed_clusters], axis=1)
#     clustered_index_df.rename(columns={clustered_index_df.columns[-1]: f'{target} cluster groups',
#                                        1: 'L-270',
#                                        2: 'R+7',
#                                        3: 'R+270'}, inplace=True)
#     melted = clustered_index_df.melt(id_vars=['astro id', f'{target} cluster groups'], 
#                                      var_name='timepoint', value_name=target)
#     return melted


# def graph_cluster_groups(df, target=None, hue=None, figsize=(7,3.2), ncol=3):
#     flatui = ["#9b59b6",  "#2ecc71", "#e74c3c", "#95a5a6", "#34495e",  "#3498db"]
    
#     plt.figure(figsize=figsize)
#     ax = sns.lineplot(x='timepoint', y=target, data=df, hue=hue,
#                       palette=sns.color_palette(flatui[:len(df[hue].unique())]),
#                       style=hue)

#     plt.setp(ax.get_xticklabels(), 
# #              rotation=45, 
#              fontsize=14)
#     if target == 'telo means':
#         ax.set_ylabel('Mean Telomere Length (Telo-FISH)', fontsize=14)
#     elif '(qPCR)' in target:
#         ax.set_ylabel('Mean Telomere Length (qPCR)', fontsize=14)
        
#     ax.set_xlabel('', fontsize=14)
#     ax.tick_params(labelsize=14)
    
#     legend = ax.legend()
#     legend.texts[0].set_text('Cluster groups')
    
#     plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18),
#           ncol=ncol, fancybox=True, fontsize=14)
    
#     plt.savefig(f'11 astronauts CLUSTERING {target}.png', 
#             dpi=600, bbox_inches = "tight")
    
    
def combine_midflight(row):
    if 'mid-flight 1' in row or 'mid-flight 2' in row:
        row = 'mid-flight'
        return row
    else:
        return row
    

def scipy_anova_post_hoc_tests(df=None, flight_status_col='flight status new',
                               sig_test=stats.f_oneway, post_hoc=sp.posthoc_ttest,
                               equal_var=False, pool_sd=False, repeated_measures=False):
    """
    df should be melted by aberration type
    """
    # make list of aberrations
    aberrations = list(df['aberration type'].unique())
    
    # loop through aberrations & perform anovas between pre/mid/post
    for aberr in aberrations:
    
        if repeated_measures == False:        
            g_1 = df[(df[flight_status_col] == 'Pre-Flight') & (df['aberration type'] == aberr)]['count per cell']
            g_2 = df[(df[flight_status_col] == 'Mid-Flight') & (df['aberration type'] == aberr)]['count per cell']
            g_3 = df[(df[flight_status_col] == 'Post-Flight') & (df['aberration type'] == aberr)]['count per cell']
            statistic, p_value = sig_test(g_1, g_2, g_3)
            print(aberr, p_value)

        elif repeated_measures:
            results = AnovaRM(df[df['aberration type'] == aberr].copy(), 'count per cell', 'astro id', 
                              within=[flight_status_col], aggregate_func='mean').fit()
            # pvalue
            p_value = results.anova_table['Pr > F'][0]


        # if anova detects sig diff, perform post-hoc tests
        if p_value <= 0.05:
            display(sp.posthoc_ttest(df[df['aberration type'] == aberr], val_col='count per cell', 
                                     group_col='flight status new', equal_var=equal_var, p_adjust='bonferroni',
                                     pool_sd=pool_sd))
            print('\n')
            
            
def rename_aberr(row):
    if row == 'sister chromatid exchanges':
        return 'classic SCEs'
    elif row == 'total inversions':
        return 'inversions'
    elif row == 'satellite associations':
        return 'sat. associations'
    else:
        return row
    
    
def rename_flights(row):
    if row == 'pre-flight':
        return 'Pre-Flight'
    elif row == 'mid-flight':
        return 'Mid-Flight'
    elif row == 'post-flight':
        return 'Post-Flight'
    
    
def pull_telofish_df():
    telof_df = pd.read_csv('../data/compiled and processed data/exploded_cells_astros_df.csv')
    telof_df_grouped = telof_df.groupby(by=['astro id', 'timepoint', 'flight status']).agg('mean').reset_index()
    telof_df_grouped['astro id'] = telof_df_grouped['astro id'].astype('int64')
    return telof_df_grouped


def pull_qpcr_df():
    # astronauts telomere qpcr df
    qpcr_df = pd.read_excel('../data/raw data/qpcr_telomere_astros.xlsx', usecols=[0, 1, 2])
    qpcr_df.dropna(axis=0, inplace=True)
    
    qpcr_df['astro id'] = qpcr_df['astro id'].astype('int64')    
    qpcr_df['flight status'] = qpcr_df['timepoint'].apply(lambda row: flight_status(row))
    
    qpcr_grouped = qpcr_df.groupby(by=['astro id', 'timepoint', 'flight status']).agg('mean').reset_index()
    qpcr_grouped['timepoint'] = qpcr_grouped['timepoint'].apply(lambda row: convert_mid_timepoint(row))
    return qpcr_grouped


def pull_aberr_df():
    melt_all_astro_chr_aberr = pd.read_csv('../data/compiled and processed data/All_astronauts_chromosome_aberration_data_tidy_data.csv')

    # reformatting (float -> int -> str)
    melt_all_astro_chr_aberr['astro id'] = melt_all_astro_chr_aberr['astro id'].astype('int')
    melt_all_astro_chr_aberr['astro id'] = melt_all_astro_chr_aberr['astro id'].astype('str')

    astro_chr_aberr = melt_all_astro_chr_aberr.copy()

    astro_chr_aberr['aberration type'] = astro_chr_aberr['aberration type'].apply(lambda row: rename_aberr(row))
    astro_chr_aberr['flight status'] = astro_chr_aberr['flight status'].apply(lambda row: combine_midflight(row))    
    astro_chr_aberr['flight status'] = astro_chr_aberr['flight status'].apply(lambda row: rename_flights(row))
    astro_chr_aberr['flight status'] = astro_chr_aberr['flight status'].astype('category')
    astro_chr_aberr['flight status'].cat.reorder_categories(['Pre-Flight', 'Mid-Flight', 'Post-Flight'], inplace=True)
    
    pivot_chr = astro_chr_aberr.pivot_table(index=['astro id', 'flight status'], 
                                            columns='aberration type', 
                                            values='count per cell')
    pivot_chr.reset_index(inplace=True)
    pivot_chr['astro id'] = pivot_chr['astro id'].astype('int64')
    return pivot_chr


def find_time_col(df1, df2):
    # combine col names into list to check fi col of interest is in both dfs
    # TO DO: check presence of cols of interest 
    dfs_columns = list(df1.columns) + list(df2.columns)
    if dfs_columns.count('flight status') == 2 and dfs_columns.count('timepoint') == 2:
        raise Exception(f'TWO TIMEPOINT REFERENTIAL COLUMNS IN merger dataframes.. please remove one')
    elif dfs_columns.count('flight status') == 2:
        return 'flight status'
    elif dfs_columns.count('timepoint') == 2:
        return 'timepoint'
    elif 'flight status' not in dfs_columns and 'timepoint' not in dfs_columns:
        raise Exception('both flight status AND timepoint cols not in merger dfs.. please check contents')
    else:
        raise Exception('one of merger dataframes lacks the desired merger col')
        
        
def pull_merge_all_data(merge_what_data=None, how_groupby_telo_data=None, what_flight_status=None,
                        what_aberrations=['dicentrics', 'inversions', 'translocations']):
    
    telofish_df = pull_telofish_df()
    qpcr_df = pull_qpcr_df()
    aberr_df = pull_aberr_df()
    
    if how_groupby_telo_data:
        telofish_df = telofish_df.groupby(how_groupby_telo_data).agg('mean').reset_index()
        qpcr_df = qpcr_df.groupby(how_groupby_telo_data).agg('mean').reset_index()

    if what_flight_status:
        telofish_df = telofish_df[telofish_df['flight status'].isin(what_telomere_flight_status)].copy()
        qpcr_df = qpcr_df[qpcr_df['flight status'].isin(what_telomere_flight_status)].copy()
        aberr_df = aberr_df[aberr_df['flight status'].isin(what_telomere_flight_status)].copy()
        
    if what_aberrations:
        aberr_df = aberr_df[['astro id', 'flight status'] + what_aberrations].copy()
        
    data_dict = {'telofish': telofish_df,
                 'qpcr': qpcr_df,
                 'aberr': aberr_df}
    
    if merge_what_data:
        # if 2 merge requests, pull dataframes w/ relevant data from dict & merge
        if len(merge_what_data) == 2:
            df0 = data_dict[merge_what_data[0]]
            df1 = data_dict[merge_what_data[1]]
            df_merged = df0.merge(df1, on=['astro id', find_time_col(df0, df1)])
            
        # if 3 merge requests, pull dataframes w/ relevant data from dict & merge
        elif len(merge_what_data) == 3:
            df0 = data_dict[merge_what_data[0]]
            df1 = data_dict[merge_what_data[1]]
            df2 = data_dict[merge_what_data[2]]
            df_temp = df0.merge(df1, on=['astro id', find_time_col(df0, df1)])
            df_merged = df_temp.merge(df2, on=['astro id', find_time_col(df_temp, df2)])
    return df_merged


