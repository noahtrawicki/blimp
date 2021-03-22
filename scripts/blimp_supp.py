import pandas as pd
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import time

# ---- INITIALIZE VARIABLES ----
lil_del_dict_eth3 = []
rmv_msg = []
rmv_meta_list = []
output_sep = '--------------'

# ---- VARIABLES for read_Nu_data function (faster to define constants outside of functions to avoid looping)----

# Get indices reference and sample side measurements
# ref_b1_idx = np.linspace(0, 40, 21, dtype = int)
# ref_b2_idx = np.linspace(41, 81, 21, dtype = int)
# ref_b3_idx = np.linspace(82, 122, 21, dtype = int)
# ref_idx = np.concatenate([ref_b1_idx, ref_b2_idx, ref_b3_idx])

sam_b1_idx = np.linspace(1, 39, 20, dtype = int)
sam_b2_idx = np.linspace(42, 80, 20, dtype = int)
sam_b3_idx = np.linspace(83, 121, 20, dtype = int)
sam_idx = np.concatenate([sam_b1_idx, sam_b2_idx, sam_b3_idx])

# ---- ASSIGN PLOT DEFAULTS ----
sns.set_palette("colorblind")
pal = sns.color_palette()
plt.rcParams['font.sans-serif'] = "Gill Sans MT"
plt.rcParams["font.family"] = "sans-serif"
medium_font = 10
plt.rc('axes', labelsize=medium_font, labelweight = 'bold')
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# ---- READ IN PARAMETERS ('params.xlsx') ----

xls = pd.ExcelFile(Path.cwd().parents[0] / 'params.xlsx')
df_anc = pd.read_excel(xls, 'Anchors', index_col = 'Anchor')
df_const = pd.read_excel(xls, 'Constants', index_col = 'Name')
df_threshold = pd.read_excel(xls, 'Thresholds', index_col = 'Type')
df_rnm_rmv = pd.read_excel(xls, 'Rename_remove', index_col = 'UID')

long_term_d47_SD = df_threshold['Value'].loc['long_term_SD']
num_SD = df_threshold['Value'].loc['num_SD']
bad_count_thresh = df_threshold['Value'].loc['bad_count_thresh']
transducer_pressure_thresh = df_threshold['Value'].loc['transducer_pressure_thresh']
balance_high = df_threshold['Value'].loc['balance_high']
balance_low = df_threshold['Value'].loc['balance_low']

calc_a18O = df_const['Value'].loc['calc_a18O']
arag_a18O = df_const['Value'].loc['arag_a18O']
dolo_a18O = df_const['Value'].loc['dolo_a18O']

Nominal_D47 = df_anc.to_dict()['D47'] # Sets anchor values for D47crunch as dictionary {Anchor: value}

# ---- DEFINE FUNCTIONS ----

def calc_bern_temp(D47_value):
	''' Calculates D47 temp using calibration from Bernasconi et al. (2018) 25C '''		
	return (((0.0449 * 1000000) / (D47_value - 0.167))**0.5) - 273.15

def calc_MIT_temp(D47_value):
	''' Calculates D47 temp using preliminary calibration from Anderson et al. (2020) 90C'''
	return (((0.039 * 1000000) / (D47_value - 0.153))**0.5) - 273.15

def calc_Petersen_temp(D47_value):
	'''Calculates D47 temperature (C) using calibration from Petersen et al. (2019) 90C'''
	return (((0.0383 * 1000000) / (D47_value - 0.170))**0.5) - 273.15

def make_water(D47_T):
	'''Calculates fluid d18O based on D47 temperature from Kim and O'Neil (1999) '''
	thousandlna = 18.03 * (1e3 * (1/(D47_T + 273.15))) - 32.42
	a = np.exp((thousandlna/1000))
	eps = (a-1) * 1e3
	return eps

def thousandlna(mineral):
		'''Calculates 18O acid fractination factor to convert CO2 d18O to mineral d18O'''
		if mineral == 'calcite' or mineral == 'Calcite':
			#a = 1.00871 # Kim (2007)
			a = calc_a18O
		elif mineral == 'dolomite' or mineral == 'Dolomite':
			#a = 1.009926 #Rosenbaum and Sheppard (1986) from Easotope
			a = dolo_a18O
		elif mineral == 'aragonite' or mineral == 'Aragonite':
			#a = 1.0090901 # Kim (2007)
			a = arag_a18O
		else:
			a = calc_a18O
			
		return 1000*np.log(a)

def read_Nu_data(data_file, file_number, current_sample):
	'''INPUT ---- Nu data file (.txt)
   OUTPUT --- list of mean d45 to d49 (i.e. little delta) values'''

   # Deals with different .txt file formats starting at UID 9628
	if file_number > 9628: 
		n_skip = 31
	else:
   		n_skip = 30
	try:		
		df = pd.read_fwf(data_file, skiprows = n_skip, header = None) # read in file, skip 31 rows, no header
	except NameError:
		print('Data file not found for UID', file_number)

	df = df.drop(columns = [0]) # removes first column (full of zeros)
	df = df.dropna(how = 'any')
	df = df.astype('float64') # make sure data is read as floats; try moving this up
	df = df[(df.T != 0).any()] # remove all zeroes; https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame	
	df = df.reset_index(drop = 'True')

	df_zero = df.head(6).astype('float64') # first 6 rows are the "Blank" i.e. zero measurement for the whole replicate
	df_zero_mean = (df_zero.apply(np.mean, axis = 1)).round(21)

	df_mean = df.mean(axis = 1)	

	# Every 6th entry is a particular mass, starting with mass 49. Starts at 6 to avoid zero measurements.
	mass_49_index = np.arange(6, len(df), 6) 
	mass_48_index = np.arange(7, len(df), 6)
	mass_47_index = np.arange(8, len(df), 6)
	mass_46_index = np.arange(9, len(df), 6)
	mass_45_index = np.arange(10, len(df), 6)
	mass_44_index = np.arange(11, len(df), 6)

	m44 = df_mean[mass_44_index] - df_zero_mean[5] # subtract mass_44 zero measurement from each mass_44 meas
	m44 = m44.dropna()
	m44 = m44.reset_index(drop = True)

	m49 = df_mean[mass_49_index] - df_zero_mean[0]
	m49 = m49.dropna()
	m49 = m49.reset_index(drop = True)
	m49_44 = m49/m44 # calculate raw 49/44 ratio

	m48 = df_mean[mass_48_index] - df_zero_mean[1]
	m48 = m48.dropna()
	m48 = m48.reset_index(drop = True)
	m48_44 = m48/m44

	m47 = df_mean[mass_47_index] - df_zero_mean[2]
	m47 = m47.dropna()
	m47 = m47.reset_index(drop = True)
	m47_44 = m47/m44

	m46 = df_mean[mass_46_index] - df_zero_mean[3]
	m46 = m46.dropna()
	m46 = m46.reset_index(drop = True)
	m46_44 = m46/m44

	m45 = df_mean[mass_45_index] - df_zero_mean[4]
	m45 = m45.dropna()
	m45 = m45.reset_index(drop = True)
	m45_44 = m45/m44

	df_zero_corr = pd.DataFrame({'m44':m44, 'm45_44':m45_44,'m46_44':m46_44, 'm47_44':m47_44, 'm48_44':m48_44, 'm49_44':m49_44})
	
	# Calculate little deltas by correcting each sample side measurement to bracketing ref side measurements
	lil_del = []

	# samp = df_zero_corr.iloc[sam_idx]
	# samp = samp.reset_index(drop = True)
	# ref = df_zero_corr.iloc[ref_idx]
	# ref = ref.reset_index(drop = True)

	for i in df_zero_corr.columns:
	    for j in sam_idx:
	        s = df_zero_corr[i][j]
	        r_prev =  df_zero_corr[i][j-1]
	        r_next = df_zero_corr[i][j+1]
	        lil_del.append(((((s/r_prev) + (s/r_next))/2.)-1)*1000)

	d45 = lil_del[60:120]
	d46 = lil_del[120:180]
	d47 = lil_del[180:240]
	d48 = lil_del[240:300]
	d49 = lil_del[300:360]

	lil_del_dict = {'d45':d45, 'd46':d46,'d47':d47, 'd48':d48, 'd49':d49}		
	df_lil_del = pd.DataFrame(lil_del_dict)

	if 'ETH' in current_sample and '3' in current_sample: # this bit is to provide raw data for joyplots/etc.
	 	lil_del_dict_eth3.extend(d47)

	return df_lil_del

# def read_Nu_results(results_file, file_number):
	
# 	''' Not currently in use -- use read_Nu_data instead!!!
# 		INPUT: Nu Results .csv files
# 		OUTPUT: small delta values'''		
# 	if __name__ == "__main__":
# 		if file_number < 7980 or file_number > 8001: # Deals with different versions of Nu software that have slightly different formats
# 			skiprows = 7
# 		else:	
# 			skiprows = 6

# 		df_results = pandas.read_csv(results_file, skiprows = skiprows)
# 		df_results.rename(columns = {'Unnamed: 0':'rep',}, inplace = True) # First column is unnamed: this changes it to 'rep'
# 		results_file.replace('_', '0') # Gets rid of Y2K issue		

# 		try:
# 			df_results = df_results[df_results.Beam !='Major'] # In some versions of Nu software, headings are created for each block; this removes extraneous headings
# 			df_results = df_results[df_results.Beam !='Beam']

# 		except(AttributeError): # This block finds the Nu software version if the formatting from the normal way of doing things doesn't work
# 			skip_list = list(range(142))
# 			skip_list.remove(1)
# 			df_dummy = pandas.read_csv(results_file, skiprows = skip_list) 
# 			nu_vers_column = df_dummy.columns[5]
# 			nu_version = nu_vers_column[27:33]
# 			raise Exception('DataFrame has no attribute "Beam". Nu software version = ', nu_version)

# 		# Convert all data to floats (sometimes pulled in as strings)
# 		df_results = df_results.astype({'47':'float', 'Delta(45/44)':'float', 'Delta(46/44)':'float', 
# 			'Delta(47/44)':'float', 'Delta(48/44)':'float', 'Delta(49/44)':'float'})

# 		return df_results


def raw_to_D47crunch_fmt(results_file, file_number, current_sample, file_count, folder_name):
	'''INPUT: Nu Data .txt file and associated data
	   OUTPUT: List of small delta and metadata ready for D47crunch'''

	df_results = read_Nu_data(results_file, file_number, current_sample)

	SD_thresh = long_term_d47_SD*num_SD
	bad_count = 0 # Keeps track of bad cycles (cycles > 5 SD from sample mean)
	bad_rep_count = 0

	# Calculate median of all cycles.
	median_47 = df_results['d47'].median()
	 
	# # Removes any cycles that have D47 values > 5 SD away from sample mean. If more than 'bad_count_thresh' cycles violate this criteria, entire replicate is removed. 
	for i in range(len(df_results['d47'])):	
		if (df_results['d47'].iloc[i]) > ((median_47) + (SD_thresh)) or ((df_results['d47'].iloc[i]) < ((median_47) - (SD_thresh))):
			df_results['d47'].iloc[i] = np.nan	# 'Disables' cycle

			bad_count += 1

	#session = 'Session' + str(file_count) # Use this line to make session = run
	session = str(folder_name[:8])
	
	rmv_analyses = [] # analysis to be removed
	
	this_path = Path.cwd() / 'raw_data' / folder_name

	# This goes through batch summary data and checks values against thresholds from params.xlsx
	for i in os.listdir(this_path):
		if 'Batch Results.csv' in i and 'fail' not in os.listdir(this_path): # checks for results summary file, doesn't run if one of the samples failed (THIS DOESNT WORK)
			summ_file = Path.cwd() / 'raw_data' / folder_name / i
			#Read in Batch Results file, combine the two rows of column headers into one		
			df_results_summ = pd.read_csv(summ_file, encoding = 'latin1', skiprows = 3, header = [0,1])
			df_results_summ.columns = df_results_summ.columns.map('_'.join).str.strip()		

			#Get the index location of the row that corresponds to the given file number (i.e. replicate)
			curr_row = df_results_summ.loc[df_results_summ['Data_File'].str.contains(str(file_number))].index		
			
			transduc_press = float(df_results_summ['Transducer_Pressure'][curr_row])				
			samp_weight = float(df_results_summ['Sample_Weight'][curr_row])
			NuCarb_temp = float(df_results_summ['Ave_Temperature'][curr_row])
			pumpover = float(df_results_summ['MaxPumpOverPressure_'][curr_row])
			init_beam = float(df_results_summ['Initial_Sam Beam'][curr_row])
			balance = float(df_results_summ['Balance_%'][curr_row])
			vial_loc = float(df_results_summ['Vial_Location'][curr_row])

			meta_data_list = [file_number, transduc_press, samp_weight, NuCarb_temp, pumpover, init_beam, balance, vial_loc, bad_count]

			# Remove any replicates that fail thresholds, write a message to the terminal
			if transduc_press < transducer_pressure_thresh:
				rmv_analyses.append(file_number)
				rmv_msg.append((str(rmv_analyses[0]) + ' failed transducer pressure requirements (transducer_pressure = ' + str(round(transduc_press,1)) + ')' ))
			if balance > balance_high or balance < balance_low:
				rmv_analyses.append(file_number)
				rmv_msg.append((str(rmv_analyses[0]) + ' failed balance requirements (balance = ' +  str(round(balance,1)) + ')'))
			if bad_count > bad_count_thresh:
				rmv_analyses.append(file_number)
				rmv_msg.append((str(rmv_analyses[0]) + ' failed cycle-level reproducibility requirements (bad cycles = ' + str(bad_count) + ')'))
	# If replicate doesn't fail any thresholds, calculate the mean lil delta and return as a list
	if bad_count < bad_count_thresh and file_number not in rmv_analyses:
		d45_avg = format(df_results['d45'].mean(), 'f')
		d46_avg = format(df_results['d46'].mean(), 'f')
		d47_avg = format(df_results['d47'].mean(), 'f')
		d48_avg = format(df_results['d48'].mean(), 'f')
		d49_avg = format(df_results['d49'].mean(), 'f')
					
		data_list = [file_number, session, current_sample, d45_avg, d46_avg, d47_avg, d48_avg, d49_avg]

		return data_list, meta_data_list
	# If replicate fails any threshold, return list with nans for little deltas and add in metadata
	else: 
		data_list = [file_number, session, current_sample, np.nan, np.nan, np.nan, np.nan, np.nan]
		rmv_meta_list.append(meta_data_list)
		return None, None
	# Counts how many bad replicates in total... not currently used
	if bad_count >= bad_count_thresh:
		bad_rep_count +=1

def fix_names(df):
	'''Changes names of standards and samples to uniform entries based on conversion spreadsheet (names_to_change.csv)'''		

	df['Sample'] = df['Sample'].str.strip() # strip whitespace
	ntc_file = Path.cwd() / 'names_to_change.csv' 
	df_new = pd.read_csv(ntc_file, encoding = 'latin1')
	for i in range(len(df_new)):
		df['Sample']=df['Sample'].str.replace(df_new['old_name'][i], df_new['new_name'][i]) # replace old with new

	def change_anchor_name(old, new, d47_low, d47_high, d46_low, d46_high):
		'''Fixes mistake of labelling IAEA-C1 as IAEA-C2 or vice-versa'''
		for i in range(len(df)):
			if df['Sample'].iloc[i] == old: # If sample is labelled e.g., 'IAEA-C1'
				if float(df['d47'].iloc[i]) > d47_low and float(df['d47'].iloc[i]) < d47_high and float(df['d46'].iloc[i]) > d46_low and float(df['d46'].iloc[i]) < d46_high:
					df['Sample'].iloc[i] = new
					rmv_msg.append(old + ' changed to ' + new + ' for analysis ' + str(df['UID'].iloc[i]))

	change_anchor_name('IAEA-C2', 'IAEA-C1', 15, 18, 10, 13)
	change_anchor_name('IAEA-C1', 'IAEA-C2', 3, 6, -2, 1)

	dir_path_fixed = Path.cwd() / 'results' / 'raw_deltas.csv' # write new names to file
	df.to_csv(dir_path_fixed, index = False)


def run_D47crunch():
	# This pulls in Mathieu Daeron's D47crunch package (https://github.com/mdaeron/D47crunch),
	# passes crunched data to it, and makes output files
	import D47crunch

	results_path = Path.cwd() / 'results'

	data = D47crunch.D47data()
	data.Nominal_D47 = Nominal_D47
	data.ALPHA_18O_ACID_REACTION = 1.00871 # From Kim/Oneil 2007 GCA Eq. 3 with 70 deg C rxn temp
	data.read(Path.cwd() / 'results' / 'raw_deltas.csv') # INPUT	

	n_anal = len(data)
	n_samp = len({r["Sample"] for r in data})
	n_sess = len({r["Session"] for r in data})

	print(output_sep)
	print('Data contains:')
	print(n_anal, 'analyses')
	print(n_samp,  'samples')
	print(n_sess, 'sessions')
	print(output_sep) 

	data.wg()
	data.crunch()
	data.standardize()

	repeatability_all = data.repeatability['r_D47']
	rpt_d13C = data.repeatability['r_d13C_VPDB']
	rpt_d18O = data.repeatability['r_d18O_VSMOW']	

	data.table_of_sessions(verbose = True, print_out = True, dir = results_path, filename = 'sessions.csv', save_to_file = True)
	data.table_of_samples(verbose = True, print_out = True, dir = results_path, save_to_file = True, filename = 'samples.csv')
	data.table_of_analyses(print_out = False, dir = results_path, save_to_file = True, filename = 'analyses.csv')
	#data.plot_sessions(dir = Path.cwd() / 'plots' / 'session_plots')

	summ_dict = {'n_sessions':n_sess, 'n_samples':n_samp, 'n_analyses':n_anal,
	'Nominal_D47_ETH-1':data.Nominal_D47['ETH-1'], 'Nominal_D47_ETH-2':data.Nominal_D47['ETH-2'],
	'Nominal_D47_ETH-3':data.Nominal_D47['ETH-3'], 'Nominal_D47_ETH-4':data.Nominal_D47['ETH-4'],
	'Nominal_D47_IAEA-C2':data.Nominal_D47['IAEA-C2'], 'Nominal_D47_MERCK':data.Nominal_D47['MERCK'],
	'Reprod_d13C': rpt_d13C, 'Reprod_d18O':rpt_d18O, 'Reprod_D47':repeatability_all}

	df = pd.DataFrame(summ_dict, index=[0])
	df.to_csv(Path.cwd()/ 'results' / 'project_info.csv', index = False)

	print('Anchors are ', data.Nominal_D47)	

	print(output_sep)
	for i in rmv_msg: print(i) # print replicates that failed threshold
	print(output_sep)

	df = pd.DataFrame(rmv_meta_list, columns = ['UID', 'Transducer_Pressure', 'Sample_Weight', 'NuCarb_temp', 'Pumpover_Pressure', 'Initial_Sam', 'Balance', 'Vial_Location', 'Bad_count'])
	save_path =  Path.cwd() / 'results' / 'rmv_analyses.csv'
	df.to_csv(save_path, index = False)

	return repeatability_all

def add_metadata(dir_path, rptability, batch_data_list):
	'''Merges data from a user-specified spreadsheet to the output of D47crunch with "Sample" as key'''
	
	file = Path.cwd() / 'results' / 'samples.csv'
	file_meta = Path.cwd() / 'metadata.xlsx'
	
	df = pd.read_csv(file, encoding = 'latin1')

	if os.path.exists(file_meta):
		df_meta = pd.read_excel(file_meta)
		df = df.merge(df_meta, how = 'left')

	def calc_meas_95(N):
		return rptability/np.sqrt(N)

	df['95% CL'] = df['95% CL'].str[2:]
	df['95% CL analysis'] = list(map(calc_meas_95, df['N']))
	df['95% CL analysis'] = round(df['95% CL analysis'], 4)

	# Calc Bernasconi et al. (2018) temperature
	#df['Bern_2018_temp'] = round(df['D47'].map(calc_bern_temp), 2)

	# Calc Anderson et al. (2020?) temperature
	df['T_MIT'] = df['D47'].map(calc_MIT_temp)
	df['T_MIT'] = round(df['T_MIT'], 1)

	# Calc Petersen et al. (2019) temperature
	df['T_Petersen'] = df['D47'].map(calc_Petersen_temp).astype('float64')
	df['T_Petersen'] = round(df['T_Petersen'], 1)

	eps = df['T_MIT'].map(make_water)
	df['d18O_water_VSMOW'] = df['d18O_VSMOW'] - eps
	df['d18O_water_VSMOW'] = round(df['d18O_water_VSMOW'], 1)

	if 'Mineralogy' in df.columns:
		df['d18O_VPDB_mineral'] = round(((df['d18O_VSMOW'] - list(map(thousandlna, df['Mineralogy']))) - 30.92)/1.03092, 1) # convert from CO2 d18O (VSMOW) to mineral d18O (VPDB)
		df['d18O_water_VSMOW'] = df['d18O_VSMOW'] - eps - list(map(thousandlna, df['Mineralogy'])) # convert from CO2  d18O VSMOW to water d18O VSMOW
		df['d18O_water_VSMOW'] = round(df['d18O_water_VSMOW'], 1)
	else:
		df['d18O_VPDB_mineral'] = round(((df['d18O_VSMOW'] - 1000*np.log(1.00871) - 30.92)/1.03092), 1) # convert from CO2 d18O (VSMOW) to calcite d18O (VPDB) if mineralogy not specified
		df['d18O_water_VSMOW'] = df['d18O_VSMOW'] - eps - 1000*np.log(1.00871) # convert from CO2  d18O VSMOW to water d18O VSMOW
		df['d18O_water_VSMOW'] = round(df['d18O_water_VSMOW'], 1)


	df.to_csv(Path.cwd() / 'results' / 'summary.csv', index = False)

	df_anal = pd.read_csv(Path.cwd() / 'results' / 'analyses.csv')
	df_batch = pd.DataFrame(batch_data_list)

	df_anal = df_anal.merge(df_meta, how = 'left', on = 'Sample')

	df_anal['Transducer_Pressure'] = df_batch[1]
	df_anal['Sample_Weight'] = df_batch[2]
	df_anal['NuCarb_temp'] = df_batch[3]
	df_anal['Pumpover_Pressure'] = df_batch[4]
	df_anal['Init_Sam_beam'] = df_batch[5]
	df_anal['Balance'] = df_batch[6]
	df_anal['Vial_Location'] = df_batch[7]

	df_anal.to_csv(Path.cwd() / 'results' / 'analyses.csv', index = False)

	os.chdir(dir_path)

def plot_ETH_D47(repeatability_all):

	from matplotlib.lines import Line2D

	df = pd.read_csv(Path.cwd() / 'analyses.csv')

	df_anchor = df.loc[(df['Sample'] == 'ETH-1') | (df['Sample'] == 'ETH-2') | 
	(df['Sample'] == 'ETH-3') | (df['Sample'] == 'ETH-4') | (df['Sample'] == 'IAEA-C2') 
	| (df['Sample'] == 'MERCK')]# | (df['Sample'] == 'IAEA-C1')]

	df_anchor = df_anchor.reset_index(drop = True)

	#  ----- PLOT D47 ANCHORS -----
	fig, ax = plt.subplots()
	ax.scatter(df_anchor['D47'], df_anchor['Sample'], color = 'white', alpha = 0.8, edgecolor = 'black')

	ax.axhline(0.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(1.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(2.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(3.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(4.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(5.5, linestyle = '--', color = 'gray', alpha = 0.5)

	label = '> 3SD external reprod.'
	for i in Nominal_D47:
		ax.scatter(Nominal_D47[i], i, marker = 'd', color = pal[2], edgecolor = 'black')
		for j in range(len(df_anchor)):
			if df_anchor['Sample'][j] == i:
				if df_anchor['D47'][j] > (Nominal_D47[i] + 3*repeatability_all) or df_anchor['D47'][j] < (Nominal_D47[i] - 3*repeatability_all):
					ax.scatter(df_anchor['D47'][j], df_anchor['Sample'][j], color = 'red', alpha = 1, edgecolor = 'black', label = label)
					ax.text(df_anchor['D47'][j] + 0.005, df_anchor['Sample'][j], df_anchor['UID'][j], zorder = 6, family = 'sans-serif')
					label = None
	plt.xlabel('D47 I-CDES')
	plt.legend()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'anchor_D47.png')
	plt.close()

	# ---- PLOT OFFSET OF ANCHORS FROM NOMINAL ----- 
	fig, ax = plt.subplots()
	
	for j in range(len(df_anchor)):
		if df_anchor['Sample'][j] == 'ETH-1': col, nom_D47 = pal[0], Nominal_D47['ETH-1']
		if df_anchor['Sample'][j] == 'ETH-2': col, nom_D47 = pal[1], Nominal_D47['ETH-2']
		if df_anchor['Sample'][j] == 'ETH-3': col, nom_D47 = pal[2], Nominal_D47['ETH-3']
		if df_anchor['Sample'][j] == 'ETH-4': col, nom_D47 = pal[3], Nominal_D47['ETH-4']
		if df_anchor['Sample'][j] == 'IAEA-C2': col, nom_D47 = pal[4], Nominal_D47['IAEA-C2']
		if df_anchor['Sample'][j] == 'MERCK': col, nom_D47 = pal[5], Nominal_D47['MERCK']
					
		ax.scatter(df_anchor['UID'][j], df_anchor['D47'][j] - nom_D47, color = col, alpha = 1, edgecolor = 'black', zorder = 3)

	ax.axhline(repeatability_all, color = 'black')
	ax.axhline(-1*repeatability_all, color = 'black')
	leg_elem = [Line2D([0], [0], marker = 'o', markerfacecolor = pal[0], color='w', label = 'ETH-1'),
				Line2D([0], [0], marker = 'o', markerfacecolor = pal[1], color='w', label = 'ETH-2'),
				Line2D([0], [0], marker = 'o', markerfacecolor = pal[2], color='w', label = 'ETH-3'),
				Line2D([0], [0], marker = 'o', markerfacecolor = pal[3], color='w', label = 'ETH-4'),
				Line2D([0], [0], marker = 'o', markerfacecolor = pal[4], color='w', label = 'IAEA-C2'),
				Line2D([0], [0], marker = 'o', markerfacecolor = pal[5], color='w', label = 'MERCK'),
				Line2D([0], [0], color='black', lw=4, label='1SD reprod.')]
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)
	plt.xlabel('UID')
	plt.ylabel('D47 - Nominal')
	ax.legend(handles = leg_elem)
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'anchor_D47_offset.png')
	plt.close()

	# ---- PLOT D47RAW OFFSET FROM NOMINAL ---- 
	fig, ax = plt.subplots()
	
	for j in range(len(df_anchor)):
		if df_anchor['Sample'][j] == 'ETH-1': col, nom_D47 = pal[0], Nominal_D47['ETH-1']
		if df_anchor['Sample'][j] == 'ETH-2': col, nom_D47 = pal[1], Nominal_D47['ETH-2']
		if df_anchor['Sample'][j] == 'ETH-3': col, nom_D47 = pal[2], Nominal_D47['ETH-3']
		if df_anchor['Sample'][j] == 'ETH-4': col, nom_D47 = pal[3], Nominal_D47['ETH-4']
		if df_anchor['Sample'][j] == 'IAEA-C2': col, nom_D47 = pal[4], Nominal_D47['IAEA-C2']
		if df_anchor['Sample'][j] == 'MERCK': col, nom_D47 = pal[5], Nominal_D47['MERCK']
					
		ax.scatter(df_anchor['UID'][j], df_anchor['D47raw'][j] - nom_D47, color = col, alpha = 1, edgecolor = 'black')

	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)
	plt.xlabel('UID')
	plt.ylabel('D47raw - Nominal')
	ax.legend(handles = leg_elem[0:5])
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'anchor_D47raw_offset.png')
	plt.close()


	# ------ PLOT D47 ALL --------

	fig, ax = plt.subplots(figsize = (7, 12))
	ax.scatter(df['D47'], df['Sample'], alpha = 0.8, edgecolor = 'black', label = 'Unknown')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	for i in range(len(df)):
		if "ETH" in df.Sample.iloc[i] or "IAEA" in df.Sample.iloc[i] or "MERCK" in df.Sample.iloc[i]:
			ax.scatter(df.D47.iloc[i], df.Sample.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9, label = 'Anchor')
	plt.xlabel('D47 I-CDES')	
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'all_D47.png')
	plt.close()

	# ----- PLOT d13C/d18O -------

	d13C_median = df.groupby('Sample')['d13C_VPDB'].median()
		
	fig, ax = plt.subplots(figsize = (7, 12))

	for i in range(len(df)):
		samp = df['Sample'][i]
		ax.scatter(df['d13C_VPDB'][i] - d13C_median[samp], samp, color = pal[0], edgecolor = 'black')
	plt.xlabel('d13C VPDB offset from median')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'd13C.png')

	d18O_median = df.groupby('Sample')['d18O_VSMOW'].median()

	fig, ax = plt.subplots(figsize = (7, 12))
	for i in range(len(df)):
		samp = df['Sample'][i]
		ax.scatter(df['d18O_VSMOW'][i] - d18O_median[samp], samp, color = pal[0], edgecolor = 'black')
	plt.xlabel('d18O VSMOW offset from median')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'd18O.png')
	
def joy_plot():

	try:
		import joypy
		from matplotlib import cm

		rep = list(range(round(len(lil_del_dict_eth3)/60)))
		new_rep = [ele for ele in rep for i in range(60)]

		df = pd.DataFrame()
		df['d47'] = lil_del_dict_eth3
		df['rep'] = new_rep
		
		fig, axes = joypy.joyplot(df, by = 'rep', column = 'd47', ylim = 'own', ylabels = False, range_style = 'own',colormap = cm.Blues, linewidth = 0.1, x_range = [16,18])
		plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.9,
                			wspace=None, hspace=None)
		#plt.title('ETH-3 (outliers not removed)')
		plt.xlabel('d47')
		plt.savefig(Path.cwd().parents[0] / 'plots' / 'ETH-3_joyplot.png')

	except ModuleNotFoundError:
		print('Joypy module not found. Please install joypy.')

def cdv_plots():

	file = Path.cwd().parents[0] / 'results' / 'analyses.csv'
	df = pd.read_csv(file, encoding = 'latin1')
	file = Path.cwd().parents[0] / 'results' / 'rmv_analyses.csv'
	df_rmv = pd.read_csv(file, encoding = 'latin1')

	plt.figure(figsize=(10,7))
	use_baseline = 'n'

# Plots transducer pressure vs. sample weight for baseline and your run. 
	plt.subplot(2,2,1)
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	plt.xlim(350, 550)
	plt.ylim(10, 40)
	plt.xlabel('Sample weight')
	plt.ylabel("Transducer pressure (mbar)")
	if use_baseline == 'y':
		plt.scatter(df_baseline.sample_weight.iloc[3:], df_baseline.transducer_pressure.iloc[3:], color = baseline_col, alpha = baseline_opac, zorder = 3, label = 'Baseline')	
	plt.scatter(df.Sample_Weight.iloc[1:], df.Transducer_Pressure.iloc[1:], color = pal[1], alpha = 1, zorder = 6, label = 'Unknown')
	# Put dark circles around the ETH
	for i in range(len(df.Sample_Weight)): 
		if "ETH" in df.Sample.iloc[i] or "IAEA" in df.Sample.iloc[i] or "MERCK" in df.Sample.iloc[i]:
			plt.scatter(df.Sample_Weight.iloc[i], df.Transducer_Pressure.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9)
	plt.scatter(df.Sample_Weight.iloc[i], df.Transducer_Pressure.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9, label = 'Anchor') # dummy for legend
	plt.scatter(df_rmv.Sample_Weight, df_rmv.Transducer_Pressure, color = 'red', edgecolor = 'black', zorder = 12, label = 'Removed')
	for i in range(len(df_rmv)):
		plt.text(df_rmv.Sample_Weight.iloc[i] + 10, df_rmv.Transducer_Pressure.iloc[i], str(df_rmv.UID.iloc[i]), zorder = 6)
	plt.legend()

	# plots max pumpover pressure vs. sample weight for basline and your run
	plt.subplot(2,2,2)
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	if use_baseline == 'y':
		plt.scatter(df_baseline.sample_weight.iloc[3:], df_baseline.max_pump.iloc[3:], color = baseline_col, alpha = baseline_opac, zorder = 3, label = 'Baseline')
	plt.scatter(df.Sample_Weight.iloc[1:], df.Pumpover_Pressure.iloc[1:], color = pal[1], alpha = 1, zorder = 6, label = 'Unknown')
	
	plt.xlabel('Sample weight')
	plt.ylabel("Max pumpover pressure (mbar)")
	
	# Put dark circles around the ETH
	for i in range(len(df.Sample_Weight)):
		if "ETH" in df.Sample.iloc[i] or "IAEA" in df.Sample.iloc[i] or "MERCK" in df.Sample.iloc[i]:
			plt.scatter(df.Sample_Weight.iloc[i], df.Pumpover_Pressure.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9)
	plt.scatter(df.Sample_Weight.iloc[i], df.Pumpover_Pressure.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9, label = 'Anchor') # dummy for legend
	plt.scatter(df_rmv.Sample_Weight, df_rmv.Pumpover_Pressure, color = 'red', edgecolor = 'black', zorder = 12, label = 'Removed')
	plt.legend()

	# Plots Balance % vs. vial location for baseline and your run
	plt.subplot(2,2,3)
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	if use_baseline == 'y':
		plt.scatter(df_baseline.UID.iloc[3:], df_baseline.Balance.iloc[3:], color = baseline_col, alpha = baseline_opac, zorder = 3, label = 'Baseline')
	plt.scatter(df.UID.iloc[1:], df.Balance.iloc[1:], color = pal[1], alpha = 1, zorder = 6, label = 'Unknown')
	plt.xlabel("UID")
	plt.ylabel("Balance %")
	#plt.xlim(0, 50)
	
	# Put dark circles around the ETH
	for i in range(len(df.UID)):
		if "ETH" in df.Sample.iloc[i] or "IAEA" in df.Sample.iloc[i] or "MERCK" in df.Sample.iloc[i]:
			plt.scatter(df.UID.iloc[i], df.Balance.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9)
	
	plt.scatter(df.UID.iloc[i], df.Balance.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9, label = 'Anchor') # dummy for legend
	plt.scatter(df_rmv.UID, df_rmv.Balance, color = 'red', edgecolor = 'black', zorder = 12, label = 'Removed')
	for i in range(len(df_rmv)):
		plt.text(df_rmv.UID.iloc[i] + 1, df_rmv.Balance.iloc[i], str(df_rmv.UID.iloc[i]), zorder = 6)
	plt.legend()

	# Plots D49 vs. vial location for baseline and your run
	plt.subplot(2,2,4)
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	if use_baseline == 'y':
		plt.scatter(df_baseline.UID.iloc[3:], df_baseline.D49.iloc[3:], color = baseline_col, alpha = baseline_opac, zorder = 3, label = 'Baseline')
	plt.scatter(df.UID.iloc[1:], df.D49raw.iloc[1:], color = pal[1], alpha = 1, zorder = 6, label = 'Unknown')
	plt.xlabel("UID")
	#plt.xlim(0, 50)
	plt.ylabel('D49')
	
	# Put dark circles around the ETH
	for i in range(len(df.UID)):
		if "ETH" in df.Sample.iloc[i] or "IAEA" in df.Sample.iloc[i] or "MERCK" in df.Sample.iloc[i]:
			plt.scatter(df.UID.iloc[i], df.D49raw.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9)
	plt.scatter(df.UID.iloc[i], df.D49raw.iloc[i], color = pal[1], edgecolor = 'black', linewidth = .75, zorder = 9, label = 'Anchor') # dummy for legend

	plt.legend()

	#plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'cdv.png')

	plt.close()

	# ---- IAEA-C1 PLOT ----
	plt.figure(figsize = (6, 4))
	df = df.loc[df['Sample'] == 'IAEA-C1']
	plt.scatter(df['UID'], df['D47'], color = pal[0], edgecolor = 'black')
	plt.axhline(0.3, color = pal[2])
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	plt.xlabel('UID')
	plt.ylabel('D47')
	plt.title('IAEA-C1 D47')
	plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'IAEA-C1_reprod.png')



