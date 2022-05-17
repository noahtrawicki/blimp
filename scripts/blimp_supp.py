# --- VERSION 0.1.5 updated 20220517 by NTA ---

import pandas as pd
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import time



# ---- INITIALIZE VARIABLES ----
lil_del_dict_eth3 = []
lil_del_dict_eth3_UID = []
rmv_msg = []
rmv_meta_list = []
output_sep = '--------------'

# ---- VARIABLES for read_Nu_data function (faster to define constants outside of functions to avoid looping)----

# Get indices reference and sample side measurements
# ref_b1_idx = np.linspace(0, 40, 21, dtype = int)
# ref_b2_idx = np.linspace(41, 81, 21, dtype = int)
# ref_b3_idx = np.linspace(82, 122, 21, dtype = int)
# ref_idx = np.concatenate([ref_b1_idx, ref_b2_idx, ref_b3_idx])

# ---- ASSIGN PLOT DEFAULTS ----

sns.set_palette("colorblind")
pal = sns.color_palette()
medium_font = 10
plt.rc('axes', labelsize=medium_font, labelweight = 'bold')
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# ---- READ IN PARAMETERS ('params.xlsx') ----

xls = pd.ExcelFile(Path.cwd().parents[0] / 'params.xlsx')
df_anc = pd.read_excel(xls, 'Anchors', index_col = 'Anchor')
df_const = pd.read_excel(xls, 'Constants', index_col = 'Name')
df_threshold = pd.read_excel(xls, 'Thresholds', index_col = 'Type')
df_rnm = pd.read_excel(xls, 'Rename_by_UID')

long_term_d47_SD = df_threshold['Value'].loc['long_term_SD']
num_SD = df_threshold['Value'].loc['num_SD']
SD_thresh = long_term_d47_SD*num_SD # pulled from parameters file
bad_count_thresh = df_threshold['Value'].loc['bad_count_thresh']
transducer_pressure_thresh = df_threshold['Value'].loc['transducer_pressure_thresh']
balance_high = df_threshold['Value'].loc['balance_high']
balance_low = df_threshold['Value'].loc['balance_low']

calc_a18O = df_const['Value'].loc['calc_a18O']
arag_a18O = df_const['Value'].loc['arag_a18O']
dolo_a18O = df_const['Value'].loc['dolo_a18O']

Nominal_D47 = df_anc.to_dict()['D47'] # Sets anchor values for D47crunch as dictionary {Anchor: value}

# Every 6th entry is a particular mass, starting with mass 49. Starts at 6 to avoid zero measurements. 
# Used in 'read_Nu_data' function but defining here to avoid looping over definitions
mass_49_index = np.arange(6, 744, 6) 
mass_48_index = np.arange(7, 744, 6)
mass_47_index = np.arange(8, 744, 6)
mass_46_index = np.arange(9, 744, 6)
mass_45_index = np.arange(10, 744, 6)
mass_44_index = np.arange(11, 744, 6)
# ---- DEFINE FUNCTIONS USED TO CALCULATE D47/T/D18Ow ----

def calc_bern_temp(D47_value):
	''' Calculates D47 temp using calibration from Bernasconi et al. (2018) 25C '''		
	return (((0.0449 * 1000000) / (D47_value - 0.167))**0.5) - 273.15

def calc_MIT_temp(D47_value):
	''' Calculates D47 temp using Eq. 1 calibration from Anderson et al. (2021) 90C'''
	if D47_value > 0.153: #(prevents complex returns)
		return (((0.0391 * 1000000) / (D47_value - 0.154))**0.5) - 273.15
	else:
		return np.nan

def calc_Petersen_temp(D47_value):
	'''Calculates D47 temperature (C) using calibration from Petersen et al. (2019) 90C'''
	return (((0.0383 * 1000000) / (D47_value - 0.170))**0.5) - 273.15

def make_water_KON97(D47_T):
	'''Calculates fluid d18O based on D47 temperature from Kim and O'Neil (1997)'''
	thousandlna_KON97 = 18.03 * (1e3 * (1/(D47_T + 273.15))) - 32.42
	a_KON97 = np.exp((thousandlna_KON97/1000))
	eps_KON97 = (a_KON97-1) * 1e3

	return eps_KON97

def make_water_A21(D47_T):
	'''Calculates fluid d18O based on D47 temperature from Anderson et al. (2021)'''
	thousandlna_A21 = 17.5 * (1e3 * (1/(D47_T + 273.15))) - 29.1
	a_A21 = np.exp((thousandlna_A21/1000))
	eps_A21 = (a_A21-1) * 1e3
	
	return eps_A21

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

# ---- Define helper functions ----

# import fnmatch
# def find_file(pattern, path):  # DOESNT WORK YET -- MIGHT BE MORE EFFICIENT
# 	result = []
# 	for files in os.walk(path):
# 		print(files)
# 		if isinstance(files, str):
# 			print(files)
# 			if fnmatch.fnmatch(files, pattern):
# 				print('found')
# 				result.append(os.path.join(root, files))
# 	return result

# ---- READ AND CORRECT DATA ----

def read_Nu_data(data_file, file_number, current_sample, folder_name, run_type):
	'''
	PURPOSE: Read in raw voltages from Nu data file (e.g., Data_13553 ETH-1.txt), zero correct, calculate R values, and calculate little deltas; remove bad cycles and use batch data to remove bad replicates
	INPUTS: Path to Nu data file (.txt); analysis UID (e.g., 10460); sample name (e.g., ETH-1); and run type ('standard' or 'clumped')
	OUTPUT: List of mean d45 to d49 (i.e. little delta) values as Pandas dataframe
	'''
	session = str(folder_name[:8]) # creates name of session; first 8 characters of folder name are date of run start per our naming convention (e.g., 20211008 clumped apatite NTA = 20211008) 
	bad_count = 0 # Keeps track of bad cycles (cycles > 5 SD from sample mean)
	bad_rep_count = 0 # Keeps track of bad replicates

   # -- Read in file --
   # Deals with different .txt file formats starting at UID 1899, 9628 (Nu software updates)
	if file_number > 9628: n_skip = 31
	elif file_number < 1899: n_skip = 29
	else: n_skip = 30

	try:	
		df = pd.read_fwf(data_file, skiprows = n_skip, header = None) # Read in file, skip n_skip rows, no header
	except NameError:
		print('Data file not found for UID', file_number)

	# -- Clean up data -- 
	df = df.drop(columns = [0]) # removes first column (full of zeros)
	df = df.dropna(how = 'any')
	df = df.astype('float64') # make sure data is read as floats
	df = df[(df.T != 0).any()] # remove all zeroes; https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame	
	df = df.reset_index(drop = 'True')

	if len(df) != 744: # valid Nu results files should have 744 lines (at least as currently written); this will skip any files that violate that criteria
		print("Input file ", data_file, "is incomplete or incorrectly formatted; it has been skipped.")
		data_list = [file_number, session, current_sample, np.nan, np.nan, np.nan, np.nan, np.nan]
		return None, None

	# -- Read in blank i.e. zero measurement -- 
	df_zero = df.head(6).astype('float64') # first 6 rows are the "Blank" i.e. zero measurement; used to zero-correct entire replicate
	df_zero_mean = (df_zero.apply(np.mean, axis = 1)).round(21) # calculates mean of zero for each mass
	df_mean = df.mean(axis = 1)	# calculates the mean of each row (i.e., averages each individual measurement to calculate a cycle mean) 

	# -- Calculate R values --

	# subtract mass_44 zero measurement from each mass_44 meas
	m44 = (df_mean[mass_44_index] - df_zero_mean[5]).dropna().reset_index(drop = True)

	# For all masses, subtract zero measurement from actual measurement, and then calc raw 4X/49 ratio.
	def calc_4X44_ratio(mass_index, zero_index):
		m4x = (df_mean[mass_index] - df_zero_mean[zero_index]).dropna().reset_index(drop = True)
		return m4x/m44

	# Create a zero-corrected dataframe of R values	
	df_zero_corr = pd.DataFrame({'m44':m44, 'm45_44':calc_4X44_ratio(mass_45_index, 4),'m46_44':calc_4X44_ratio(mass_46_index, 3), 
								'm47_44':calc_4X44_ratio(mass_47_index, 2), 'm48_44':calc_4X44_ratio(mass_48_index, 1),
								'm49_44':calc_4X44_ratio(mass_49_index, 0)})

	# Calculate little deltas (d4X) by correcting each sample side measurement to bracketing ref side measurements
	lil_del = []

	# if clumped run, index locations of all sample side cycles are defined at top of script so they are not redefined with every analysis
	sam_b1_idx = np.linspace(1, 39, 20, dtype = int)
	sam_b2_idx = np.linspace(42, 80, 20, dtype = int)
	sam_b3_idx = np.linspace(83, 121, 20, dtype = int)
	sam_idx = np.concatenate([sam_b1_idx, sam_b2_idx, sam_b3_idx])

	# if standard run, index locations of sample side cycles are different
	if run_type == 'standard':
		sam_idx = np.linspace(1, 11, 6, dtype = int)

	# compare sample measurement to bracketing ref gas measurement
	for i in df_zero_corr.columns:
	    for j in sam_idx: # 'sam_idx' defined near top of script
	        # df_zero_corr[i][j] is the sample side
	        # df_zero_corr[i][j-1] is the previous ref side
	        # df_zero_corr[i][j+1] is the following ref side

	        lil_del.append(((((df_zero_corr[i][j]/df_zero_corr[i][j-1]) + (df_zero_corr[i][j]/df_zero_corr[i][j+1]))/2.)-1)*1000)

	# Define each little delta value by index position
	if run_type == 'clumped':	
		d45 = lil_del[60:120]
		d46 = lil_del[120:180]
		d47 = lil_del[180:240]
		d48 = lil_del[240:300]
		d49 = lil_del[300:360]

	elif run_type == 'standard':
		d45 = lil_del[6:12]
		d46 = lil_del[12:18]
		d47 = lil_del[18:24]
		d48 = lil_del[24:30]
		d49 = lil_del[30:36]

	lil_del_dict = {'d45':d45, 'd46':d46,'d47':d47, 'd48':d48, 'd49':d49}		
	df_lil_del = pd.DataFrame(lil_del_dict) # export to dataframe -- makes it easier for next function to handle

	if 'ETH' in current_sample and '3' in current_sample: # this bit is to provide raw data for joyplots/etc.
	 	lil_del_dict_eth3.extend(d47)
	 	# for i in range(len(lil_del_dict_eth3)):
	 	lil_del_dict_eth3_UID.append(file_number)	 


	batch_data_list = [file_number, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

	# Calculate median of all cycles.
	median_47 = df_lil_del['d47'].median()
	d47_pre_SE = df_lil_del['d47'].sem()
	 
	# -- FIND BAD CYCLES --  
	# Removes any cycles that have d47 values > 5 SD away from sample median. If more than 'bad_count_thresh' cycles violate this criteria, entire replicate is removed.
	if run_type == 'clumped':
		for i in range(len(df_lil_del['d47'])):

			# If d47 is outside threshold, remove little deltas of ALL masses for that cycle (implemented 20210819)	
			if (df_lil_del['d47'].iloc[i]) > ((median_47) + (SD_thresh)) or ((df_lil_del['d47'].iloc[i]) < ((median_47) - (SD_thresh))):
				df_lil_del['d45'].iloc[i] = np.nan # 'Disables' cycle; sets value to nan
				df_lil_del['d46'].iloc[i] = np.nan
				df_lil_del['d47'].iloc[i] = np.nan
				df_lil_del['d48'].iloc[i] = np.nan
				df_lil_del['d49'].iloc[i] = np.nan	

				bad_count += 1

	session = str(folder_name[:8]) # creates name of session; first 8 characters of folder name are date of run start per our naming convention (e.g., 20211008 clumped apatite NTA = 20211008) 
	
	d47_post_SE = df_lil_del['d47'].sem()

	rmv_analyses = [] # analysis to be removed
	
	this_path = Path.cwd() / 'raw_data' / folder_name

	# -- Find bad replicates -- 
	# This goes through batch summary data and checks values against thresholds from params.xlsx
	for i in os.listdir(this_path):
		if 'Batch Results.csv' in i and 'fail' not in os.listdir(this_path): # checks for and reads results summary file 
			summ_file = Path.cwd() / 'raw_data' / folder_name / i # i = e.g., 20210505 clumped dolomite apatite calibration 5 NTA Batch Results.csv
			df_results_summ = pd.read_csv(summ_file, encoding = 'latin1', skiprows = 3, header = [0,1])
			df_results_summ.columns = df_results_summ.columns.map('_'.join).str.strip()	# fixes weird headers of Nu Summary files

			#Get the index location of the row that corresponds to the given file number (i.e. replicate)
			curr_row = df_results_summ.loc[df_results_summ['Data_File'].str.contains(str(file_number))].index
			batch_data_list = [file_number, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, d47_pre_SE, d47_post_SE, bad_count]
			if len(curr_row) == 1 and run_type == 'clumped': # curr_row is Int64Index, which acts like a list. If prev line finds either 0 or 2 matching lines, it will skip this section.
				transduc_press = float(df_results_summ['Transducer_Pressure'][curr_row])				
				samp_weight = float(df_results_summ['Sample_Weight'][curr_row])
				NuCarb_temp = float(df_results_summ['Ave_Temperature'][curr_row])
				pumpover = float(df_results_summ['MaxPumpOverPressure_'][curr_row])
				init_beam = float(df_results_summ['Initial_Sam Beam'][curr_row])
				balance = float(df_results_summ['Balance_%'][curr_row])
				vial_loc = float(df_results_summ['Vial_Location'][curr_row])
				d13C_SE = float(df_results_summ['Std_Err.5'][curr_row])
				d18O_SE = float(df_results_summ['Std_Err.6'][curr_row])
				D47_SE = float(df_results_summ['Std_Err.7'][curr_row])

				batch_data_list = [file_number, transduc_press, samp_weight, NuCarb_temp, pumpover, init_beam, balance, vial_loc, d13C_SE, d18O_SE, D47_SE, d47_pre_SE, d47_post_SE, bad_count]

				# Remove any replicates that fail thresholds, compile a message that will be written to the terminal
				if transduc_press < transducer_pressure_thresh:
					rmv_analyses.append(file_number)
					rmv_msg.append((str(rmv_analyses[0]) + str(current_sample) + ' failed transducer pressure requirements (transducer_pressure = ' + str(round(transduc_press,1)) + ')' ))
				if balance > balance_high or balance < balance_low:
					rmv_analyses.append(file_number)
					rmv_msg.append((str(rmv_analyses[0]) + str(current_sample) + ' failed balance requirements (balance = ' +  str(round(balance,1)) + ')'))
				if bad_count > bad_count_thresh:
					rmv_analyses.append(file_number)
					rmv_msg.append((str(rmv_analyses[0]) + str(current_sample) + ' failed cycle-level reproducibility requirements (bad cycles = ' + str(bad_count) + ')'))

			break # Found a matching file? There only should be one, so stop here.

		else: # Couldn't find matching UID, or got confused. No batch summary data included.
			batch_data_list = [file_number, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, d47_pre_SE, d47_post_SE, bad_count]
	
	# If replicate doesn't fail any thresholds, calculate the mean lil delta and return as a list
	if bad_count < bad_count_thresh and file_number not in rmv_analyses:
		d45_avg = format(df_lil_del['d45'].mean(), 'f')
		d46_avg = format(df_lil_del['d46'].mean(), 'f')
		d47_avg = format(df_lil_del['d47'].mean(), 'f')
		d48_avg = format(df_lil_del['d48'].mean(), 'f')
		d49_avg = format(df_lil_del['d49'].mean(), 'f')
					
		data_list = [file_number, session, current_sample, d45_avg, d46_avg, d47_avg, d48_avg, d49_avg]

		return data_list, batch_data_list

	# If replicate fails any threshold, return list with nans for little deltas and add in metadata
	else: 
		data_list = [file_number, session, current_sample, np.nan, np.nan, np.nan, np.nan, np.nan]
		batch_data_list.append(current_sample)	
		rmv_meta_list.append(batch_data_list)
		return None, None

def fix_names(df):
	'''
	PURPOSE: Changes names of standards and samples to uniform entries based on conversion spreadsheet
	INPUT: Pandas Dataframe of little deltas, names_to_change tab in params.csv
	OUTPUT: Fully corrected, accurately named little deltas (raw_deltas.csv)'''
	
	df['Sample'] = df['Sample'].str.strip() # strip whitespace
	df_new = pd.read_excel(xls, 'Names_to_change')
	
	# rename based on name (names_to_change; i.e. EHT-1 --> ETH-1)
	for i in range(len(df_new)):
		df['Sample']=df['Sample'].str.replace(df_new['old_name'][i], df_new['new_name'][i])

	# rename samples based on UID (Rename_by_UID; i.e. whatever the name of 10155 is, change to 'ETH-1')
	if len(df_rnm) > 0: # check if there's anything to rename
		for i in range(len(df_rnm)):
			rnm_loc = np.where(df['UID'] == df_rnm['UID'][i])[0] # get index location of particular UID			
			df.loc[rnm_loc, 'Sample'] = df_rnm.loc[i, 'New_name'] # replace sample name in main df with sample name from rnm_rmv

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


def run_D47crunch(run_type):
	''' 
	PURPOSE: Calculate D47, d13C, d18O, using Mathieu Daeron's 'D47_crunch' package (https://github.com/mdaeron/D47crunch)	
	INPUT: Fully corrected little deltas ('raw_deltas.csv')
	OUTPUT: repeatability (1SD) of all calculated measurements; also writes 'sessions.csv', 'analyses.csv', and 'samples.csv',
	  '''
	import D47crunch
	
	results_path = Path.cwd() / 'results'

	data = D47crunch.D47data()
	data.Nominal_D47 = Nominal_D47
	data.ALPHA_18O_ACID_REACTION = calc_a18O # This is selected from the params file -- you can use whatever value you want in there.

	#values from Bernasconi et al 2018 Table 4
	if run_type == 'standard':
		data.SAMPLE_CONSTRAINING_WG_COMPOSITION = ('ETH-1', 2.02, -2.19) # oftentimes for standard racks, we don't use ETH-3, so this uses ETH-1 as the anchor
	elif run_type == 'clumped':
		data.SAMPLE_CONSTRAINING_WG_COMPOSITION = ('ETH-3', 1.71, -1.78)

	print('Sample constraining WG composition = ', data.SAMPLE_CONSTRAINING_WG_COMPOSITION)
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

	if run_type == 'clumped':
		data.standardize()
		repeatability_all = data.repeatability['r_D47']
		rpt_d13C = data.repeatability['r_d13C_VPDB']
		rpt_d18O = data.repeatability['r_d18O_VSMOW']

		# display and save session info as csv
		data.table_of_sessions(verbose = True, print_out = True, dir = results_path, filename = 'sessions.csv', save_to_file = True)
		
		sam = data.table_of_samples(verbose = True, print_out = True, save_to_file = False)
		analy = data.table_of_analyses(print_out = False, save_to_file = False)
		#data.plot_sessions(dir = Path.cwd() / 'plots' / 'session_plots') # Issue on everyones computer but Noah's...

		# send sample and analyses to pandas dataframe, clean up
		df_sam = pd.DataFrame(sam[1:], columns = sam[0]).replace(r'^\s*$', np.nan, regex=True)  #import as dataframe, replace empty strings with NaN
		df_sam['95% CL'] = df_sam['95% CL'].str[2:] # clean up plus/minus signs
		df_sam = df_sam.astype({'Sample':'str', 'N':'int32', 'd13C_VPDB':'float64', 'd18O_VSMOW':'float64', 'D47':'float64','SE':'float64', 
						'95% CL':'float64', 'SD':'float64', 'p_Levene':'float64'}) #recast types appropriately (all str by default)
		df_sam = df_sam.rename(columns = {'95% CL': 'CL_95_pct'})
		
		df_analy = pd.DataFrame(analy[1:], columns = analy[0]).replace(r'^\s*$', np.nan, regex=True)  #import as dataframe, replace empty strings with NaN
		df_analy = df_analy.astype({'UID':'int32', 'Session':'int32', 'Sample':'str', 'd13Cwg_VPDB':'float64', 
			'd18Owg_VSMOW':'float64', 'd45':'float64', 'd46':'float64', 'd47':'float64', 'd48':'float64', 
			'd49':'float64', 'd13C_VPDB':'float64', 'd18O_VSMOW':'float64', 'D47raw':'float64', 'D48raw':'float64',
       		'D49raw':'float64', 'D47':'float64'}) #recast types appropriately (all str by default)

		df_rmv = pd.read_excel('params.xlsx', 'Remove')
		manual_rmv = list(df_rmv.UID)
		manual_rmv_reason = list(df_rmv.Notes)

		# Create a little dictionary to store info about the project
		summ_dict = {'n_sessions':n_sess, 'n_samples':n_samp, 'n_analyses':n_anal,
		'Nominal_D47_ETH-1':data.Nominal_D47['ETH-1'], 'Nominal_D47_ETH-2':data.Nominal_D47['ETH-2'],
		'Nominal_D47_ETH-3':data.Nominal_D47['ETH-3'], 'Nominal_D47_ETH-4':data.Nominal_D47['ETH-4'],
		'Nominal_D47_IAEA-C2':data.Nominal_D47['IAEA-C2'], 'Nominal_D47_MERCK':data.Nominal_D47['MERCK'],
		'Reprod_d13C': rpt_d13C, 'Reprod_d18O':rpt_d18O, 'Reprod_D47':repeatability_all, 'Long_term_SD_threshold':long_term_d47_SD, 
		'Num_SD_threshold':num_SD, 'Bad_cycles_threshold':bad_count_thresh, 'Transducer_pressure_threshold': transducer_pressure_thresh,
		'Balance_high_threshold':balance_high,'Balance_low_threshold':balance_low, 'Manually_removed':manual_rmv, 'Manually_removed_reason': manual_rmv_reason}

		df_prj_summ = pd.DataFrame([summ_dict], index = [0])
		# df_prj_summ['Manually_removed'] = df_rmv['UID']
		# df_prj_summ['Manually_removed_reason'] = df_rmv['Notes']
		df_prj_summ.to_csv(Path.cwd()/ 'results' / 'project_info.csv', index = False)	

		print('Anchors are ', data.Nominal_D47)	# list anchors used and their nominal D47

		print(output_sep)
		for i in rmv_msg: print(i) # print replicates that failed threshold
		print(output_sep)

		print('Total # analyses removed = ', len(rmv_meta_list), '(', round((len(rmv_meta_list)/n_anal)*100,1), '% of total)')

		# For reps that failed, make csv with all parameters they could have failed on
		df = pd.DataFrame(rmv_meta_list, columns = ['UID', 'Transducer_Pressure', 'Sample_Weight', 'NuCarb_temp', 'Pumpover_Pressure', 'Initial_Sam', 'Balance', 'Vial_Location', 'd13C_SE (Nu)', 'd18O_SE (Nu)', 'D47_SE (Nu)', 'd47_pre_SE', 'd47_post_SE', 'Bad_count', 'Sample'])
		save_path =  Path.cwd() / 'results' / 'rmv_analyses.csv'
		df.to_csv(save_path, index = False)

		return df_sam, df_analy, repeatability_all

	# If it's a standard run, use Noah's reworking of Mathieu's code
	elif run_type == 'standard':
		table_of_analyses_std(data, print_out = False, dir = results_path, save_to_file = True, filename = 'analyses_bulk.csv')

		return np.nan

def table_of_analyses_std(data, dir = 'results', filename = 'analyses.csv', save_to_file = True, print_out = True):
        '''
        Print out an/or save to disk a table of analyses. Modified by NTA to just print out 'standard' (d13C and d18O) data

        __Parameters__

        + `dir`: the directory in which to save the table
        + `filename`: the name to the csv file to write to
        + `save_to_file`: whether to save the table to disk
        + `print_out`: whether to print out the table
        '''
        from D47crunch import make_csv
        out = [['UID','Session','Sample']]
       
        out[-1] += ['d13Cwg_VPDB','d18Owg_VSMOW','d45','d46','d47','d48','d49','d13C_VPDB','d18O_VSMOW']
        for r in data:
                out += [[f"{r['UID']}",f"{r['Session']}",f"{r['Sample']}"]]
                out[-1] += [
                        f"{r['d13Cwg_VPDB']:.3f}",
                        f"{r['d18Owg_VSMOW']:.3f}",
                        f"{r['d45']:.6f}",
                        f"{r['d46']:.6f}",
                        f"{r['d47']:.6f}",
                        f"{r['d48']:.6f}",
                        f"{r['d49']:.6f}",
                        f"{r['d13C_VPDB']:.6f}",
                        f"{r['d18O_VSMOW']:.6f}"]
        if save_to_file:
                if not os.path.exists(dir):
                        os.makedirs(dir)
                with open(f'{dir}/{filename}', 'w') as fid:
                        fid.write(make_csv(out))
        if print_out:
               	pass
        return out

def add_metadata(dir_path, rptability, batch_data_list, df, df_anal):
	'''
	PURPOSE: Merges sample metadata from 'Metadata' sheet in params.csv to the output of D47crunch with "Sample" as key;
	Calculate T, error on T, d18Ow (based on mineralogy)
	INPUT: 
	OUTPUT:
	'''	

	file_meta = Path.cwd() / 'params.xlsx'
	if os.path.exists(file_meta):
		df_meta = pd.read_excel(file_meta, 'Metadata')
		df = df.merge(df_meta, how = 'left')

	def calc_meas_95(N):
		return rptability/np.sqrt(N)

	df['CL_95_pct_analysis'] = list(map(calc_meas_95, df['N']))
	df['CL_95_pct_analysis'] = round(df['CL_95_pct_analysis'], 4)

	# Calc Bernasconi et al. (2018) temperature; I wouldn't recommend using this unless you also change the nominal D47 values of the anchors to Bern 2018 values
	#df['Bern_2018_temp'] = round(df['D47'].map(calc_bern_temp), 2)

	# Calc Anderson et al. (2021) temperature
	df['T_MIT'] = df['D47'].map(calc_MIT_temp)
	
	df['CL_95_pct'] = df['CL_95_pct'].astype('float64')
	df['SE'] = df['SE'].astype('float64')

	# Calculate upper and lower 95 CL temperatures (as the 'magnitude' of the error bar -- so 20 C sample with 10 C upper 95CL = 30 C upper 95 CL value)
	T_MIT_95CL_lower = df['D47'] + df['CL_95_pct']
	df['T_MIT_95CL_lower'] = round(abs(df['T_MIT'] - T_MIT_95CL_lower.map(calc_MIT_temp)), 1)

	T_MIT_95CL_upper = df['D47'] - df['CL_95_pct']
	df['T_MIT_95CL_upper'] = round(abs(df['T_MIT'] - T_MIT_95CL_upper.map(calc_MIT_temp)), 1)

	T_MIT_SE_lower = df['D47'] + df['SE']
	df['T_MIT_SE_lower'] = round(abs(df['T_MIT'] - T_MIT_SE_lower.map(calc_MIT_temp)), 1)

	T_MIT_SE_upper = df['D47'] - df['SE']
	df['T_MIT_SE_upper'] = round(abs(df['T_MIT'] - T_MIT_SE_upper.map(calc_MIT_temp)), 1)

	# Calc Petersen et al. (2019) temperature
	df['T_Petersen'] = df['D47'].map(calc_Petersen_temp).astype('float64')
	df['T_Petersen'] = round(df['T_Petersen'], 1)

	T_MIT_95CL_upper_val = df['T_MIT'] + df['T_MIT_95CL_upper']
	T_MIT_95CL_lower_val = df['T_MIT'] - df['T_MIT_95CL_lower']

	eps_KON97 = df['T_MIT'].map(make_water_KON97)
	eps_KON97_upper = T_MIT_95CL_upper_val.map(make_water_KON97)
	eps_KON97_lower = T_MIT_95CL_lower_val.map(make_water_KON97)

	eps_A21 = df['T_MIT'].map(make_water_A21)
	eps_A21_upper = T_MIT_95CL_upper_val.map(make_water_A21)
	eps_A21_lower = T_MIT_95CL_lower_val.map(make_water_A21)
		
	df['T_MIT'] = round(df['T_MIT'], 1)

	def calc_d18Ow(eps):

		if 'Mineralogy' in df.columns:
			df['d18O_VPDB_mineral'] = round(((df['d18O_VSMOW'] - list(map(thousandlna, df['Mineralogy']))) - 30.92)/1.03092, 1) # convert from CO2 d18O (VSMOW) to mineral d18O (VPDB)
			d18Ow_VSMOW = round(df['d18O_VSMOW'] - eps - list(map(thousandlna, df['Mineralogy'])),1) # convert from CO2  d18O VSMOW to water d18O VSMOW

		else:
			df['d18O_VPDB_mineral'] = round(((df['d18O_VSMOW'] - 1000*np.log(1.00871) - 30.92)/1.03092), 1) # convert from CO2 d18O (VSMOW) to calcite d18O (VPDB) if mineralogy not specified
			d18Ow_VSMOW = round(df['d18O_VSMOW'] - eps - 1000*np.log(1.00871),1) # convert from CO2  d18O VSMOW to water d18O VSMOW using Kim and O'Neil (1997)

		return d18Ow_VSMOW

	df['d18Ow_VSMOW_KON97'] = calc_d18Ow(eps_KON97)
	df['d18Ow_VSMOW_KON97_upper'] = round(abs(df['d18Ow_VSMOW_KON97'] - calc_d18Ow(eps_KON97_upper)), 1)
	df['d18Ow_VSMOW_KON97_lower'] = round(abs(df['d18Ow_VSMOW_KON97'] - calc_d18Ow(eps_KON97_lower)), 1)

	df['d18Ow_VSMOW_A21'] = calc_d18Ow(eps_A21)
	df['d18Ow_VSMOW_A21_upper'] = round(abs(df['d18Ow_VSMOW_A21'] - calc_d18Ow(eps_A21_upper)), 1)
	df['d18Ow_VSMOW_A21_lower'] = round(abs(df['d18Ow_VSMOW_A21'] - calc_d18Ow(eps_A21_lower)), 1)

	df.to_csv(Path.cwd() / 'results' / 'summary.csv', index = False)

	df_batch = pd.DataFrame(batch_data_list, columns = ['UID', 'Transducer_Pressure', 'Sample_Weight', 'NuCarb_temp','Pumpover_Pressure',
		'Init_Sam_beam', 'Balance', 'Vial_Location', 'd13C_SE (Nu)', 'd18O_SE (Nu)', 'D47_SE (Nu)', 'd47_pre_SE', 'd47_post_SE', 'Bad_Cycles'])

	df_anal['T_MIT'] = df_anal['D47'].map(calc_MIT_temp)
	eps_KON97 = df_anal['T_MIT'].map(make_water_KON97)
	eps_A21 = df_anal['T_MIT'].map(make_water_A21)

	df_anal['T_MIT'] = round(df_anal['T_MIT'], 1)

	if 'Mineralogy' in df_anal.columns:
		df_anal['d18O_VPDB_mineral'] = round(((df_anal['d18O_VSMOW'] - list(map(thousandlna, df_anal['Mineralogy']))) - 30.92)/1.03092, 1) # convert from CO2 d18O (VSMOW) to mineral d18O (VPDB)
		df_anal['d18O_water_VSMOW_KON97'] = round(df_anal['d18O_VSMOW'] - eps_KON97 - list(map(thousandlna, df_anal['Mineralogy'])),1) # convert from CO2  d18O VSMOW to water d18O VSMOW
		df_anal['d18O_water_VSMOW_A21'] = round(df_anal['d18O_VSMOW'] - eps_A21 - list(map(thousandlna, df_anal['Mineralogy'])),1) # convert from CO2  d18O VSMOW to water d18O VSMOW

	else:
		df_anal['d18O_VPDB_mineral'] = round(((df_anal['d18O_VSMOW'] - 1000*np.log(1.00871) - 30.92)/1.03092), 1) # convert from CO2 d18O (VSMOW) to calcite d18O (VPDB) if mineralogy not specified
		df_anal['d18O_water_VSMOW_KON97'] = round(df_anal['d18O_VSMOW'] - eps_KON97 - 1000*np.log(1.00871),1) # convert from CO2  d18O VSMOW to water d18O VSMOW
		df_anal['d18O_water_VSMOW_A21'] = round(df_anal['d18O_VSMOW'] - eps_A21 - 1000*np.log(1.00871),1) # convert from CO2  d18O VSMOW to water d18O VSMOW
		
	df_anal = df_anal.merge(df_meta, how = 'left', on = 'Sample')
	df_anal = df_anal.merge(df_batch, how = 'left', on = 'UID')

	n_bad_cycles = df_anal['Bad_Cycles'].sum()
	print('Total # bad cycles removed = ', n_bad_cycles, '(', round((n_bad_cycles/(len(df_anal)*60))*100, 2), '%)') # does not include bad cycles from disabled reps

	eth_loc = np.where(df_anal['Sample'] == 'ETH-4')	
	mbar_mg_eth = (df_anal['Transducer_Pressure'].iloc[eth_loc] / df_anal['Sample_Weight'].iloc[eth_loc]).median()
	df_anal['pct_evolved_carbonate'] = round(((df_anal['Transducer_Pressure'] / df_anal['Sample_Weight']) / mbar_mg_eth)*100, 1)

	df_anal.to_csv(Path.cwd() / 'results' / 'analyses.csv', index = False)
	to_earthchem(df_anal)

	os.chdir(dir_path)

def add_metadata_std():

	file_meta = Path.cwd() / 'params.xlsx'
	df_anal = pd.read_csv(Path.cwd() / 'results' / 'analyses_bulk.csv')

	if os.path.exists(file_meta):
		df_meta = pd.read_excel(file_meta, 'Metadata')
		df_anal= df_anal.merge(df_meta, how = 'left')	

	if 'Mineralogy' in df_anal.columns:
		df_anal['d18O_VPDB_mineral'] = round(((df_anal['d18O_VSMOW'] - list(map(thousandlna, df_anal['Mineralogy']))) - 30.92)/1.03092, 1) # convert from CO2 d18O (VSMOW) to mineral d18O (VPDB)

	else:
		df_anal['d18O_VPDB_mineral'] = round(((df_anal['d18O_VSMOW'] - 1000*np.log(1.00871) - 30.92)/1.03092), 1) # convert from CO2 d18O (VSMOW) to calcite d18O (VPDB) if mineralogy not specified

	df_anal.to_csv(Path.cwd() / 'results' / 'analyses_bulk.csv', index = False)
	

def to_earthchem(df_a):
	''' Formats analyses in format used for the EarthChem database'''

	df_ec = pd.DataFrame()

	df_ec['SampName'] = df_a['Sample']

	df_ec['SampCategory'] = 'null'
	df_ec['SampSubCategory'] = 'null'
	df_ec['SampNum'] = 'NA'
	df_ec['Mineralogy'] = df_a['Mineralogy']
	df_ec['Date'] = df_a['Session']
	df_ec['AnalysisID'] = df_a['UID']
	df_ec['RefYN'] = 'NA'
	df_ec['D47TE_SG_WD'] = 'NA'
	df_ec['MassSpec'] = 'Nu Perspective'
	# df_ec['FormT'] = df_a['Form_T']
	# df_ec['erFormT'] = df_a['err_Form_T']
	df_ec['rxnTemp'] = 70
	df_ec['Bad'] = 0

	for i in range(len(df_a['Sample'])):
		this = df_a['Sample'].loc[i] 
		
		if this == 'ETH-1' or this == 'ETH-2' or this == 'ETH-3' or this == 'ETH-4' or this == 'MERCK' or this == 'IAEA-C1' or this == 'IAEA-C2':

			df_ec['SampCategory'] = df_ec['SampCategory'].replace('null', 'carbSTD')
			df_ec['SampSubCategory'] = df_ec['SampSubCategory'].replace('null', this)

			if this != 'IAEA-C1':
				df_ec['RefYN'] = df_ec['RefYN'].replace('NA', 'Y')
			else:
				df_ec['RefYN'] = df_ec['RefYN'].replace('NA', 'N')
			if this == 'ETH-1':
				df_ec['D47TE_SG_WD'] = df_ec['D47TE_SG_WD'].replace('NA', 0.2052)
			if this == 'ETH-2':
				df_ec['D47TE_SG_WD'] = df_ec['D47TE_SG_WD'].replace('NA', 0.2085)
			if this == 'ETH-3':
				df_ec['D47TE_SG_WD'] = df_ec['D47TE_SG_WD'].replace('NA', 0.6132)
			if this == 'ETH-4':
				df_ec['D47TE_SG_WD'] = df_ec['D47TE_SG_WD'].replace('NA', 0.4505)
			if this == 'IAEA-C2':
				df_ec['D47TE_SG_WD'] = df_ec['D47TE_SG_WD'].replace('NA', 0.6409)
			if this == 'MERCK':
				df_ec['D47TE_SG_WD'] = df_ec['D47TE_SG_WD'].replace('NA', 0.5135)

		else:
			
			df_ec['SampCategory'] = df_ec['SampCategory'].replace('null','sample')
			# NB need subcategory for non-standards
			#df_ec['SampSubCategory'][i] = df_a['Method'][i]
			df_ec['RefYN'] = df_ec['RefYN'].replace('NA', 'N')
			#df_ec['D47TE_SG_WD'].loc[i] = 'NA'

	df_ec['AFF_WD'] = 'NA'
	df_ec['ARF_ID1'] = df_a['Session']
	df_ec['d45'] = df_a['d45']
	df_ec['d46'] = df_a['d46']
	df_ec['d47'] = df_a['d47']
	df_ec['d48'] = df_a['d48']
	df_ec['d49'] = df_a['d49']
	df_ec['d13C_wg_VPDB'] = df_a['d13Cwg_VPDB']
	df_ec['d18O_wg_VSMOW'] = df_a['d18Owg_VSMOW']
	df_ec['BRd13C'] = df_a['d13C_VPDB']
	df_ec['BRd18O'] = df_a['d18O_VSMOW']
	df_ec['BRD47'] = df_a['D47raw']
	df_ec['BRD48'] = df_a['D48raw']
	df_ec['BRD49'] = df_a['D49raw']
	df_ec['BRSlopeEGL'] = 'NA'
	df_ec['BRSlopeETF_WD'] = 'NA'
	df_ec['BRIntETF_WD'] = 'NA'
	df_ec['BRD47rfac_P_newAFF'] = df_a['D47']
	df_ec['d18Oac'] = 'NA' #df_a['d18O_mineral_VPDB']

	df_ec.to_csv(Path.cwd() / 'results' / 'analyses_earthchem_fmt.csv', index = False)

# ---- MAKE PLOTS ----

def plot_ETH_D47(repeatability_all, df):

	from matplotlib.lines import Line2D

	df = pd.read_csv(Path.cwd() / 'analyses.csv')

	df_anchor = df.loc[(df['Sample'] == 'ETH-1') | (df['Sample'] == 'ETH-2') | 
	(df['Sample'] == 'ETH-3') | (df['Sample'] == 'ETH-4') | (df['Sample'] == 'IAEA-C2') 
	| (df['Sample'] == 'MERCK')]# | (df['Sample'] == 'IAEA-C1')]

	df_anchor = df_anchor.reset_index(drop = True)

	#  ----- PLOT D47 ANCHORS -----
	fig, ax = plt.subplots()
	ax.scatter(df_anchor['D47'], df_anchor['Sample'], color = 'gray', alpha = 0.8, edgecolor = 'black')
	ax.axhline(0.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(1.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(2.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(3.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(4.5, linestyle = '--', color = 'gray', alpha = 0.5)
	ax.axhline(5.5, linestyle = '--', color = 'gray', alpha = 0.5)

	label = '> 3SD external reprod.; *not automatically disabled*'
	for i in Nominal_D47:
		ax.scatter(Nominal_D47[i], i, marker = 'd', color = 'white', edgecolor = 'black')
		for j in range(len(df_anchor)):
			if df_anchor['Sample'][j] == i:
				if df_anchor['D47'][j] > (Nominal_D47[i] + 3*repeatability_all) or df_anchor['D47'][j] < (Nominal_D47[i] - 3*repeatability_all):
					ax.scatter(df_anchor['D47'][j], df_anchor['Sample'][j], color = 'orange', alpha = 1, edgecolor = 'black', label = label)
					ax.text(df_anchor['D47'][j] + 0.005, df_anchor['Sample'][j], df_anchor['UID'][j], zorder = 6)
	plt.xlabel('D47 I-CDES')
	#plt.legend()
	plt.tight_layout()
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

	fig_ht = len(df)*0.02 + 1

	fig, ax = plt.subplots(figsize = (7, fig_ht))
	ax.scatter(df['D47'], df['Sample'], alpha = 0.8, color = 'white', edgecolor = 'black', label = 'Unknown')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	for i in range(len(df)):
		if "ETH" in df.Sample.iloc[i] or "IAEA" in df.Sample.iloc[i] or "MERCK" in df.Sample.iloc[i]:
			ax.scatter(df.D47.iloc[i], df.Sample.iloc[i], color = 'gray', edgecolor = 'black', linewidth = .75, zorder = 9, label = 'Anchor')
	plt.xlabel('D47 I-CDES')	
	plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'all_D47.png')

	plt.close()

	# ----- PLOT d13C/d18O -------

	d13C_median = df.groupby('Sample')['d13C_VPDB'].median()
		
	fig, ax = plt.subplots(figsize = (7, fig_ht))

	for i in range(len(df)):
		samp = df['Sample'][i]
		ax.scatter(df['d13C_VPDB'][i] - d13C_median[samp], samp, color = 'white', edgecolor = 'black')
	plt.xlabel('d13C VPDB offset from median')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'd13C.png')

	d18O_median = df.groupby('Sample')['d18O_VSMOW'].median()

	fig, ax = plt.subplots(figsize = (7, fig_ht))
	for i in range(len(df)):
		samp = df['Sample'][i]
		ax.scatter(df['d18O_VSMOW'][i] - d18O_median[samp], samp, color = 'white', edgecolor = 'black')
	plt.xlabel('d18O VSMOW offset from median')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)
	plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'd18O.png')
	plt.close()

def cdv_plots(df):

	file = Path.cwd().parents[0] / 'results' / 'analyses.csv'
	df = pd.read_csv(file, encoding = 'latin1')
	file = Path.cwd().parents[0] / 'results' / 'rmv_analyses.csv'
	df_rmv = pd.read_csv(file, encoding = 'latin1')

	plt.figure(figsize=(12,9))
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
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'cdv.png')
	plt.close()

	# ---- IAEA-C1 PLOT ----
	plt.figure(figsize = (6, 4))
	df = df.loc[df['Sample'] == 'IAEA-C1']
	plt.scatter(df['UID'], df['D47'], color = 'white', edgecolor = 'black')
	plt.axhline(0.3018, color = 'black')
	plt.grid(b=True, which='major', color='gray', linestyle='--', zorder = 0, alpha = 0.4)	
	plt.xlabel('UID')
	plt.ylabel('D47')
	plt.title('IAEA-C1 D47')
	plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'IAEA-C1_reprod.png')
	plt.close()

def d47_D47_plot(df):

	df_anchor = df.loc[(df['Sample'] == 'ETH-1') | (df['Sample'] == 'ETH-2') | 
	(df['Sample'] == 'ETH-3') | (df['Sample'] == 'ETH-4') | (df['Sample'] == 'IAEA-C2') 
	| (df['Sample'] == 'MERCK')]

	plt.figure(figsize=(10,7))
	plt.scatter(df['d47'], df['D47'], color = 'black', edgecolor = 'white')
	sns.scatterplot(data = df_anchor, x = 'd47', y = 'D47', hue = 'Sample', alpha = 1, edgecolor = 'black')
	
	plt.xlabel('d47')
	plt.ylabel('D47')
	plt.tight_layout()
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'd47_D47.png')
	plt.close()

def interactive_plots(df):

	try:
		from bokeh.io import output_file
		from bokeh.io import save
		from bokeh.plotting import figure, show, ColumnDataSource
		from bokeh.models.tools import HoverTool
		from bokeh.layouts import row
		import bokeh.models as bmo
		from bokeh.palettes import d3
		from bokeh.transform import jitter
		
	except ModuleNotFoundError:
		print("Bokeh package not found. Install Bokeh for interactive plots.")

# put this in try/except clause
	#df = pd.read_csv(Path.cwd().parents[0] / 'results' / 'analyses.csv')
	output_file(filename="D47_d47_interactive.html", title="D47_d47_interactive")

	df_anchor = df.loc[(df['Sample'] == 'ETH-1') | (df['Sample'] == 'ETH-2') | 
	(df['Sample'] == 'ETH-3') | (df['Sample'] == 'ETH-4') | (df['Sample'] == 'IAEA-C2') 
	| (df['Sample'] == 'MERCK')]

	data_anchors = ColumnDataSource.from_df(df_anchor)
	data_analyses = ColumnDataSource.from_df(df)

	TOOLTIPS = [("Sample name", "@Sample"),
				("UID", "@UID")]
	std_tools = ['pan,wheel_zoom,box_zoom,reset,hover']

	palette = d3['Category10'][len(df_anchor['Sample'].unique())]
	color_map = bmo.CategoricalColorMapper(factors=df_anchor['Sample'].unique(),
                                   palette=palette)

	f1 = figure(x_axis_label = 'd47',
				y_axis_label = 'D47',
				tools = std_tools,
				tooltips = TOOLTIPS)

	f1.scatter('d47', 'D47', source = data_analyses, color = 'black')
	f1.scatter('d47', 'D47', source = data_anchors, color = {'field':'Sample', 'transform':color_map})
	
	save(f1)

	output_file(filename="D47_all_interactive.html", title="D47_all_interactive")

	sample_names = pd.unique(df['Sample'])
	TOOLTIPS = [("Sample name", "@Sample"), ("UID", "@UID"), ("d13C", "@d13C_VPDB"), ("d18O_mineral", "@d18O_VPDB_mineral")]
	std_tools = ['pan,wheel_zoom,box_zoom,reset,hover']
	f2 = figure(x_axis_label = 'D47', y_axis_label = 'Sample', y_range = sample_names, tools = std_tools, tooltips = TOOLTIPS)
	f2.circle(x = 'D47', y=jitter('Sample', width=0.1, range=f2.y_range), source = data_analyses, color = 'white', line_color = 'black', size = 7)
	
	save(f2)


def joy_plot():

	sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
	sns.set_context()

	rep = list(range(round(len(lil_del_dict_eth3)/60)))
	new_rep = [ele for ele in rep for i in range(60)]

	df = pd.DataFrame()
	df['d47'] = lil_del_dict_eth3
	df['rep'] = new_rep
	
	cpal = sns.cubehelix_palette(10, rot=-.25, light=.7)
	g = sns.FacetGrid(df, row="rep", hue="rep", aspect=15, height=.5, palette=cpal)

	# Draw the densities in a few steps
	g.map(sns.kdeplot, "d47",
	      bw_adjust=.5, clip_on=False,
	      fill=True, alpha=1, linewidth=1.5)
	g.map(sns.kdeplot, "d47", clip_on=False, color="w", lw=2, bw_adjust=.5)

	# passing color=None to refline() uses the hue mapping
	g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

	# Define and use a simple function to label the plot in axes coordinates
	def label(x, color, label):
	    ax = plt.gca()
	    ax.text(0, .2, label, fontweight="bold", color=color,
	            ha="left", va="center", transform=ax.transAxes)

	g.map(label, "d47")

	# Set the subplots to overlap
	g.figure.subplots_adjust(hspace=-.75)

	# Remove axes details that don't play well with overlap
	g.set_titles("")
	g.set(yticks=[], ylabel="")
	g.despine(bottom=True, left=True)
	plt.xlabel('Uncorrected ETH-3 d47')
	plt.xlim(15, 17)
	plt.savefig(Path.cwd().parents[0] / 'plots' / 'ridgeplot_eth3.png')
	plt.style.use('default')  

	plt.close()


	