# clumped data vetting V2 (cdv2)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from bokeh.io import output_file
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.layouts import row

# -- GET AND CLEAN BATCH FILE --

dir_path = input("Enter the folder path of your results files (drag & drop): ").replace(r'\ ',' ').lstrip().rstrip().replace(r'\\ ', r'\ ').replace('"', '')
os.chdir(dir_path) # switch current directory to user-specified path
filenames = os.listdir() # list all the files in the working directory
Nu_Results_file = str([filename for i,filename in enumerate(filenames) if 'Results.csv' in filename]).replace("['",'').replace("']",'') # finds the name of the Nu Results file; converts from list to string; includes some additional character clean-up

df = pd.read_csv(Nu_Results_file, encoding = 'latin1', skiprows = 3, header = [0,1]) # read Batch file into pandas, skip first 3 rows, combine first two rows into single row of headers
df.columns = df.columns.map('_'.join).str.strip()	# fixes weird headers of Nu Summary files
df['Data_File'] = df['Data_File'].str.slice(start = 5, stop = 10) # changes from e.g., Data_11111 to 11111, makes easier for comparing to blimp

# -- INTERACTIVE PLOTS OF BATCH PARAMETERS WITH BOKEH --

output_file('bokeh_TEST.html') # where to save file (will be in current directory)
data = ColumnDataSource.from_df(df) # data source

TOOLTIPS = [("Sample name", "@Sample_Name"),
			("UID", "@Data_File")] # what to display with hover
std_tools = ['pan,wheel_zoom,box_zoom,reset,hover'] # which tools to include

f1 = figure(x_axis_label = 'Sample_Weight', 
	y_axis_label = 'Transducer pressure (mbar)',
	tools= std_tools,
	tooltips=TOOLTIPS)
f1.scatter('Sample_Weight', 'Transducer_Pressure', source = data)

f2 = figure(x_axis_label = 'Vial Number', 
	y_axis_label = 'Balance %',
	tools=std_tools,
	tooltips=TOOLTIPS)
f2.scatter('Vial_Location', 'Balance_%', source = data)

show(row(f1, f2))

exit()

# -- GET AND CLEAN LITTLE DELTA VALUES FROM DATA.TXT FILES --

d47_crunch_fmt = []
batch_data_list = []

# NOTES FOR NOAH -- FIX PATH; CALL FUNCTION; 

def read_Nu_data(data_file, file_number, current_sample, folder_name, run_type):
	'''
	PURPOSE: Read in raw voltages from Nu data file (e.g., Data_13553 ETH-1.txt), zero correct, calculate R values, and calculate little deltas
	INPUTS: Path to Nu data file (.txt); analysis UID (e.g., 10460); sample name (e.g., ETH-1); and run type ('standard' or 'clumped')
	OUTPUT: List of mean d45 to d49 (i.e. little delta) values as Pandas dataframe
	'''
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

	# -- Read in blank i.e. zero measurement -- 
	df_zero = df.head(6).astype('float64') # first 6 rows are the "Blank" i.e. zero measurement; used to zero-correct entire replicate
	df_zero_mean = (df_zero.apply(np.mean, axis = 1)).round(21) # calculates mean of zero for each mass

	df_mean = df.mean(axis = 1)	# calculates the mean of each row (i.e., averages each individual measurement to calculate a cycle mean) 

	# Every 6th entry is a particular mass, starting with mass 49. Starts at 6 to avoid zero measurements. 
	mass_49_index = np.arange(6, len(df), 6) 
	mass_48_index = np.arange(7, len(df), 6)
	mass_47_index = np.arange(8, len(df), 6)
	mass_46_index = np.arange(9, len(df), 6)
	mass_45_index = np.arange(10, len(df), 6)
	mass_44_index = np.arange(11, len(df), 6)

	# -- Calculate R values --

	# subtract mass_44 zero measurement from each mass_44 meas
	m44 = df_mean[mass_44_index] - df_zero_mean[5]
	m44 = m44.dropna()
	m44 = m44.reset_index(drop = True)

	# For all masses, subtract zero measurement from actual measurement, and then calc 4X/49 ratio.
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

	# Create a zero-corrected dataframe of R values
	df_zero_corr = pd.DataFrame({'m44':m44, 'm45_44':m45_44,'m46_44':m46_44, 'm47_44':m47_44, 'm48_44':m48_44, 'm49_44':m49_44})
	
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
					rmv_msg.append((str(rmv_analyses[0]) + ' failed transducer pressure requirements (transducer_pressure = ' + str(round(transduc_press,1)) + ')' ))
				if balance > balance_high or balance < balance_low:
					rmv_analyses.append(file_number)
					rmv_msg.append((str(rmv_analyses[0]) + ' failed balance requirements (balance = ' +  str(round(balance,1)) + ')'))
				if bad_count > bad_count_thresh:
					rmv_analyses.append(file_number)
					rmv_msg.append((str(rmv_analyses[0]) + ' failed cycle-level reproducibility requirements (bad cycles = ' + str(bad_count) + ')'))

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