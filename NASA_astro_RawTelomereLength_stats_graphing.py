def generate_histograms_and_dataframes_forTeloLengthData(patharg):

	"""
	opens raw telomere length count excel files from imageJ analyses and
	extracts the individual mean telomere lengths to make histograms
	pass the directory where the telomere length excel files (.xlsx) are located
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
			mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.drop(labels=[6, 193, 380, 567, 754, 941, 1128, 1315,
				1502, 1689, 1876, 2063, 2250, 2437, 2624, 2811, 2998, 3185, 3372, 3559, 3746, 3933, 4120, 4307, 4494, 4681, 4868, 5055, 5242, 5429])

			mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.iloc[7:5611]
			meantelos_str_toNaN = pd.to_numeric(mean_values_of_individual_telomere_lengths, errors='coerce')
			mean_individual_telos_cleaned = meantelos_str_toNaN.dropna(axis=0, how='any')
			mean_individ_df = mean_individual_telos_cleaned.to_frame(name=None)
			mean_individ_df = mean_individ_df[(np.abs(stats.zscore(mean_individ_df)) < 3).all(axis=1)]
			

			if ('5163' in file.name) or ('1536' in file.name):
				mean_individ_df_cy3Cal = mean_individ_df.div(59.86)

			elif '2171' in file.name:
				mean_individ_df_cy3Cal = mean_individ_df.div(80.5)

			elif '7673' in file.name:
				mean_individ_df_cy3Cal = mean_individ_df.div(2.11)

			elif '2479' in file.name:
				mean_individ_df_cy3Cal = mean_individ_df.div(2.18)

			elif '1261' in file.name:
				mean_individ_df_cy3Cal = mean_individ_df.div(2.16)

			else:
				mean_individ_df_cy3Cal = mean_individ_df

			file_name_trimmed = file.name.replace('.xlsx', '')
			dict_mean_individ_telos_dfs[file_name_trimmed] = mean_individ_df_cy3Cal
			

	astro_list_of_IDs = ['5163', '2171', '1536', '7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']
	timepointSeries = ['L-270', 'L-180', 'L-60', 'FD45', 'FD90', 'FD140', 'FD260', 'R+5', 'R+7' 'R+60', 'R+180', 'R+270']


	for idNO in astro_list_of_IDs:

	# 	#initialize blank list of timepoints
		astro_list_by_timepoint = []
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

	# 	#loop through dictionary keys
		# for i in dict_mean_individ_telos_dfs.keys():
		counter = 1
	# 		#loop through timepoints
		for j in timepointSeries:
			for i in dict_mean_individ_telos_dfs.keys():
				if (idNO in i) and j == 'L-270' and ('L-270' in i):
					astro_L270 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)
					astro_L270name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and j == 'L-180' and ('L-180' in i):
					astro_L180 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)
					astro_L180name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and j == 'L-60' and ('L-60' in i):
					astro_L60 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)
					astro_L60name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and (j == 'FD45' or j == 'FD90') and (j in i):
					astro_Mid1 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)					
					astro_Mid1name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and (j == 'FD140' or j == 'FD260') and (j in i):
					astro_Mid2 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)
					astro_Mid2name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and j == 'R+7' and (j in i):
					astro_R7 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)					
					astro_R7name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and j == 'R+60' and (j in i):
					astro_R60 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)						
					astro_R60name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and j == 'R+180' and (j in i):
					astro_R180 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)					
					astro_R180name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				elif (idNO in i) and j == 'R+270' and (j in i):
					astro_R270 = dict_mean_individ_telos_dfs[i]
					astro_list_by_timepoint.append(i)				
					astro_R270name = i.replace('mphase TeloFISH', '').replace('.xlsx', '')

				else:
					counter += 1


		
		if idNO == '5163' or idNO == '2171' or idNO == '1536':
			if (astro_L270.size > 25 or astro_L180.size > 25) and (astro_Mid1.size > 25 and astro_Mid2.size > 25 ) and (astro_R180.size > 25 or astro_R270.size > 25):

				n_bins = 60
				fig, axs = plt.subplots(2,2, sharey=True, tight_layout=True, figsize = (12.8, 9.6))

				astro_L270array = astro_L270.to_numpy()
				astro_L180array = astro_L180.to_numpy()
				astro_R180array = astro_R180.to_numpy()

				if astro_L270name != '': 
					if astro_R270name != '':
						# def astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone, axsNUMtwo):
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L270, astro_L270array, astro_L270name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L270array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L270array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R270, astro_L270array, astro_R270name, 1, 1)
						# statistics_between_timepoints(astro_pre, astro_mid1=None, astro_mid2=None, astro_post, 
						# 			  astro_prename, astro_mid1name=None, astro_mid2name=None, astro_postname)

						statistics_between_timepoints(astro_L270, astro_Mid1, astro_Mid2, astro_R270, 
									  astro_L270name, astro_Mid1name, astro_Mid2name, astro_R270name)

					elif astro_R270name == '':
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L270, astro_L270array, astro_L270name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L270array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L270array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R180, astro_L270array, astro_R180name, 1, 1)

						statistics_between_timepoints(astro_L270, astro_Mid1, astro_Mid2, astro_R180, 
									  astro_L270name, astro_Mid1name, astro_Mid2name, astro_R180name)

				elif astro_L270name == '':
					if astro_R270name == '':
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L180, astro_L180array, astro_L180name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L180array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L180array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R180, astro_L180array, astro_R180name, 1, 1)

						statistics_between_timepoints(astro_L180, astro_Mid1, astro_Mid2, astro_R180, 
									  astro_L180name, astro_Mid1name, astro_Mid2name, astro_R180name)

					elif astro_R270name != '':
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L180, astro_L180array, astro_L180name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L180array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L180array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R270, astro_L180array, astro_R270name, 1, 1)

						statistics_between_timepoints(astro_L180, astro_Mid1, astro_Mid2, astro_R270, 
									  astro_L180name, astro_Mid1name, astro_Mid2name, astro_R270name)

				else:
					continue

				# plt.savefig('Final telomere histogram non random sampling dso'+idNO+'.pdf')
				# plt.show()


		if idNO == '5163' or idNO == '2171' or idNO == '1536':

			if (astro_L270.size > 25 or astro_L180.size > 25) and (astro_Mid1.size > 25 and astro_Mid2.size > 25 ) and (astro_R180.size > 25 or astro_R270.size > 25):

				astro_L270 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_L270, 'rsamp')
				astro_L180 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_L180, 'rsamp')
				astro_Mid1 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_Mid1, 'rsamp')
				astro_Mid2 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_Mid2, 'rsamp')
				astro_R180 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_R180, 'rsamp')
				astro_R270 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_R270, 'rsamp')

				astro_L270array = astro_L270.to_numpy()
				astro_L180array = astro_L180.to_numpy()
				astro_Mid1array = astro_Mid1.to_numpy()
				astro_Mid2array = astro_Mid2.to_numpy()
				astro_R180array = astro_R180.to_numpy()
				astro_R270array = astro_R270.to_numpy()





				n_bins = 60
				fig, axs = plt.subplots(2,2, sharey=True, tight_layout=True, figsize = (12.8, 9.6))

				if astro_L270name != '': 
					if astro_R270name != '':
						# astronaut_histogram_sytlizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone):
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L270, astro_L270array, astro_L270name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L270array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L270array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R270, astro_L270array, astro_R270name, 1, 1)
						print('randomly sampled stats')
						statistics_between_timepoints(astro_L270, astro_Mid1, astro_Mid2, astro_R270, 
									  astro_L270name, astro_Mid1name, astro_Mid2name, astro_R270name)


					elif astro_R270name == '':
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L270, astro_L270array, astro_L270name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L270array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L270array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R180, astro_L270array, astro_R180name, 1, 1)
						print('randomly sampled stats')
						statistics_between_timepoints(astro_L270, astro_Mid1, astro_Mid2, astro_R180, 
									  astro_L270name, astro_Mid1name, astro_Mid2name, astro_R180name)


				elif astro_L270name == '':
					if astro_R270name == '':
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L180, astro_L180array, astro_L180name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L180array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L180array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R180, astro_L180array, astro_R180name, 1, 1)
						print('randomly sampled stats')
						statistics_between_timepoints(astro_L180, astro_Mid1, astro_Mid2, astro_R180, 
									  astro_L180name, astro_Mid1name, astro_Mid2name, astro_R180name)


					elif astro_R270name != '':
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_L180, astro_L180array, astro_L180name, 0, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid1, astro_L180array, astro_Mid1name, 0, 1)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_Mid2, astro_L180array, astro_Mid2name, 1, 0)
						astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astro_R270, astro_L180array, astro_R270name, 1, 1)
						print('randomly sampled stats')
						statistics_between_timepoints(astro_L180, astro_Mid1, astro_Mid2, astro_R270, 
									  astro_L180name, astro_Mid1name, astro_Mid2name, astro_R270name)


				else:
					continue

				
				# plt.savefig('Final telomere histogram random sampling dso'+idNO+'.pdf')
				# plt.show()



		if idNO in ['7673', '4819', '3228', '2494', '2479', '2381', '1261', '1062']:

			if (astro_L270.size > 25) and (astro_R270.size > 25):

				astro_L270array = astro_L270.to_numpy()
				astro_L180array = astro_L180.to_numpy()

				n_bins = 50
				fig, axs = plt.subplots(2, sharey=True, tight_layout=True, figsize = (12.8, 9.6))

				astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astro_L270, astro_L270array, astro_L270name, 0)
				astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astro_R270, astro_L270array, astro_R270name, 1)

				statistics_between_timepoints_prepost_only(astro_L270, astro_R270, astro_L270name, astro_R270name)

			else:
				continue

			# plt.savefig('nonsampled telomere histogram dso'+idNO+'.pdf')
			# plt.show()



			if (astro_L270.size > 25) and (astro_R270.size > 25):

				astro_L270 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_L270, 'rsamp')
				astro_R270 = gen_missing_values_andimpute_or_randomsampledown(30, 184, astro_R270, 'rsamp')

				astro_L270array = astro_L270.to_numpy()
			

				n_bins = 50
				fig, axs = plt.subplots(2, sharey=True, tight_layout=True, figsize = (12.8, 9.6))

				astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astro_L270, astro_L270array, astro_L270name, 0)
				astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astro_R270, astro_L270array, astro_R270name, 1)

				statistics_between_timepoints_prepost_only(astro_L270, astro_R270, astro_L270name, astro_R270name)

			else:
				continue

			# plt.savefig('Resampled telomere histogram dso'+idNO+'.pdf')
			# plt.show()
			



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

		


def astronaut_histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone, axsNUMtwo):

		astroarray = astroDF.to_numpy()

		N, bins, patches = axs[axsNUMone,axsNUMtwo].hist(astroarray, bins=n_bins, range=(0, 400), edgecolor='black')

		for a in range(len(patches)):
			if bins[a] <= np.quantile(astroquartile, 0.25):
				patches[a].set_facecolor('#fdff38')

			elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
				patches[a].set_facecolor('#d0fefe')

			elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
				patches[a].set_facecolor('#d0fefe')

			elif bins[a] > np.quantile(astroquartile, 0.75): 
				patches[a].set_facecolor('#ffbacd')


		axs[axsNUMone,axsNUMtwo].set_title('Histogram of ' + astroname + 's Telomeres')
		axs[axsNUMone,axsNUMtwo].set_xlabel('Bins of Individ. Telomeres')
		axs[axsNUMone,axsNUMtwo].set_ylabel('Freqs of Individ. Telomeres')
		axs[axsNUMone,axsNUMtwo].xaxis.set_major_locator(plt.MaxNLocator(12))


def astronaut_histogram_stylizer_divyBins_byQuartile_2Stacked(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone):

	astroarray = astroDF.to_numpy()

	N, bins, patches = axs[axsNUMone].hist(astroarray, bins=n_bins, range=(0, 400), edgecolor='black')

	for a in range(len(patches)):
		if bins[a] <= np.quantile(astroquartile, 0.25):
			patches[a].set_facecolor('#fdff38')

		elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
			patches[a].set_facecolor('#d0fefe')

		elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
			patches[a].set_facecolor('#d0fefe')

		elif bins[a] > np.quantile(astroquartile, 0.75): 
			patches[a].set_facecolor('#ffbacd')


	axs[axsNUMone].set_title('Histogram of Individual Telomeres for ' + astroname)
	axs[axsNUMone].set_xlabel('Bins of Individ. Telomeres')
	axs[axsNUMone].set_ylabel('Freqs of Individ. Telomeres')
	axs[axsNUMone].xaxis.set_major_locator(plt.MaxNLocator(19))


def find_max_count(astro_df1, astro_df2):
	astrolist = [astro_df1.size, astro_df2.size]
	return max(astrolist)



def make_missingValue_array_based_on_difference(n_cells, telosPercell, astro_df):

	missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
	range1 = np.arange(missing_data_difference)

	missingd = pd.DataFrame(np.nan, index=(range1), columns=['Mean Individ Telos'])
	print('missing data difference between 5520 telos and current astro is ', missing_data_difference)

	concat_ed = pd.concat([missingd, astro_df], sort=False)
	np.random.shuffle(concat_ed.to_numpy())
	print(missingd.shape, astro_df.shape, concat_ed.shape)
	
	return concat_ed


def gen_missing_values_andimpute_or_randomsampledown(n_cells, telosPercell, astro_df, option=None):

	if astro_df.size > 5520:
		astro_dfsampled = astro_df.sample(5520)
		return astro_dfsampled

	if astro_df.size > 25 and astro_df.size <= 2760:
		missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
		rsampled = astro_df.sample(missing_data_difference, replace=True, random_state=28)
		concat_ed = pd.concat([rsampled, astro_df], sort=False)
		np.random.shuffle(concat_ed.to_numpy())
		return concat_ed

	if astro_df.size > 25 and astro_df.size < 5520:
		missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )

		if option == 'rsamp':
			rsampled = astro_df.sample(missing_data_difference, random_state=28)
			concat_ed = pd.concat([rsampled, astro_df], sort=False)

			np.random.shuffle(concat_ed.to_numpy())
			return concat_ed

		elif option == 'kNN':
			range1 = np.arange(missing_data_difference)
			missingd = pd.DataFrame(np.nan, index=(range1), columns=['Mean Individ Telos'])
			concat_ed = pd.concat([missingd, astro_df], sort=False)

			np.random.shuffle(concat_ed.to_numpy())
			concat_ed = concat_ed.to_numpy().reshape(69, 80)
			imputer = KNNImputer(missing_values='NaN', n_neighbors=1, row_max_missing=1, col_max_missing=1, weights="uniform")

			x = imputer.fit_transform(concat_ed)
			x = x.reshape(5520,1)
			x = pd.DataFrame(x)

			return x

		else:
			return astro_df

	else:
		return astro_df





if __name__ == '__main__':
	import sys
	import os

	import numpy as np
	from numpy import array

	import pandas as pd
	from pandas import ExcelWriter
	from pandas import ExcelFile

	import matplotlib.pyplot as plt
	from matplotlib.ticker import StrMethodFormatter
	from matplotlib import colors
	from matplotlib.ticker import PercentFormatter

	from scipy import stats
	from missingpy import KNNImputer

	from statsmodels.graphics.gofplots import qqplot
	



	patharg = sys.argv[1]
	generate_histograms_and_dataframes_forTeloLengthData(patharg)
