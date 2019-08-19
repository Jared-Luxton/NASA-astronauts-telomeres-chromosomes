#e nables access to directories/files
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
from ptitprince import PtitPrince as pt

# statistics
from scipy import stats
from statsmodels.graphics.gofplots import qqplot





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
    astro_prename, astro_mid1name, astro_mid2name, astro_postname):

    print(astro_prename + '  compared vs  ' + astro_mid1name,
            stats.mannwhitneyu(astro_pre, astro_mid1), '\n\n\n',

            astro_prename + '  compared vs  ' + astro_mid2name,
            stats.mannwhitneyu(astro_pre, astro_mid2),'\n\n\n', 
            
            astro_mid1name + '  compared vs  ' + astro_postname,
            stats.mannwhitneyu(astro_mid1, astro_post),'\n\n\n', 

            astro_mid1name + '  compared vs  ' + astro_mid2name,
            stats.mannwhitneyu(astro_mid1, astro_mid2),'\n\n\n', 

            astro_mid2name + '  compared vs  ' + astro_postname,
            stats.mannwhitneyu(astro_mid2, astro_post),'\n\n\n', 

            astro_prename + '  compared vs  ' + astro_postname,
            stats.mannwhitneyu(astro_pre, astro_post),'\n\n\n', )


    
def statistics_between_timepoints_prepost_only(astro_pre, astro_post, astro_prename, astro_postname):

    print(astro_prename + '  compared vs  ' + astro_postname,
            stats.mannwhitneyu(astro_pre, astro_post),'\n\n\n', )
    
    
    
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
            print(file.name, 'IT WORKS PEGGY!!! <3')
        
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

        
# might be important for nasa astro bar graphs

# # shape (51, 56)

# all_patients_df_copy = all_patients_df.drop([44, 45, 46], axis=0)

# #taking all_patients_df, removing the index & making a multi index of patient id and timepoint
# individ_cell_cols = all_patients_df_copy.reindex().set_index(['patient id', 'timepoint']) 

# #removing unnecessary columns
# individ_cell_cols.drop(['chr data', 'telo data', 'status', 'telo means', 'Q1', 'Q2-3', 'Q4'], axis=1, inplace=True)

# #exploding the series containing the individual telos
# explode_cells = individ_cell_cols['cell data'].apply(pd.Series)

# #transpose!
# explode_cells = explode_cells.reset_index(level=['patient id']).T

# print(explode_cells.shape)
# explode_cells.head(5)


# def evaluate_timepoint_presence(timepoint, df):
#     if df[timepoint] != '':
#         quartile = L_270
#         return quartile
    
#     elif L_180 != '':
#         quartile = L_180
#         return quartile
    
# for timepoint in [L_270, L_180, L_60, R_7, R_60, R_180, R_270]:
#     quartile = evaluate_timepoint_presence(L_270)
#     print(quartile)
#     break

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
        
        
def make_histograms_colored_by_quartile_for_astronauts(exploded_telos_df=None):

    astro_ids = ['5163', '2171', '1536', '7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']
    
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
            
            n_bins = 60

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

            graph_two_histograms(quartile_ref, n_bins, astro_L270, astro_R270,
                                                   name_L270, name_R270)
            
        plt.savefig(f'../individual telomere length histogram distributions/dso{astro_id_num} histogram of individual telomere length distributions.pdf')
            
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
    
            
def graph_four_histograms(quartile_ref, n_bins, df1, df2, df3, df4,
                                                name1, name2, name3, name4):
    
    n_bins = n_bins
    fig, axs = plt.subplots(2,2, sharey=True, tight_layout=False, figsize = (16, 12))
    sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
    
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df1, quartile_ref, name1, 0, 0)
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df2, quartile_ref, name2, 0, 1)
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df3, quartile_ref, name3, 1, 0)
    astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df4, quartile_ref, name4, 1, 1)
    
    
    
def graph_two_histograms(quartile_ref, n_bins, df1, df2,
                                               name1, name2):
    
    n_bins = n_bins
    fig, axs = plt.subplots(2, sharey=True, tight_layout=False, figsize = (16, 8))
    sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
    
    for ax in axs.flat:
        ax.label_outer()
    
    astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, df1, quartile_ref, name1, 0)
    astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, df2, quartile_ref, name2, 1)