import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path

rd_path = Path.cwd() / 'raw_data'
results_path = Path.cwd() / ('results')
plot_path = Path.cwd() / ('plots')


sam_b1_idx = np.linspace(1, 39, 20, dtype = int)
sam_b2_idx = np.linspace(42, 80, 20, dtype = int)
sam_b3_idx = np.linspace(83, 121, 20, dtype = int)
sam_idx = np.concatenate([sam_b1_idx, sam_b2_idx, sam_b3_idx])

if os.path.isdir(Path.cwd().parents[0]/ 'raw_data' / 'single_analyses'):
	print('Single analysis found.')
	for item in os.listdir(Path.cwd().parents[0]/ 'raw_data' / 'single_analyses'):
		if 'Data' in item and '.txt' in item and '.fail' not in item: # if there's a Data file not within a run folder
			file_number = int(item[5:10])
			samp_name = str(item[10:-4])

			print(item)

			if file_number > 9628:
				n_skip = 31
			elif file_number < 1899:
				n_skip = 29
			else:
		   		n_skip = 30

			try:	
				df = pd.read_fwf(Path.cwd().parents[0]/ 'raw_data' / 'single_analyses' / str(item), skiprows = n_skip, header = None) # read in file, skip 31 rows, no header
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

			 # subtract mass_44 zero measurement from each mass_44 meas
			m44 = df_mean[mass_44_index] - df_zero_mean[5]
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

			# Create a zero-corrected dataframe of R values 
			df_zc = pd.DataFrame({'m44':m44, 'm45_44':m45_44,'m46_44':m46_44, 'm47_44':m47_44, 'm48_44':m48_44, 'm49_44':m49_44})
			print('Data is zero corrected')

			lil_del = []

			for i in df_zc.columns:
			    for j in sam_idx: # 'sam_idx' defined near top of script
			        s = df_zc[i][j]
			        r_prev =  df_zc[i][j-1]
			        r_next = df_zc[i][j+1]
			        lil_del.append(((((s/r_prev) + (s/r_next))/2.)-1)*1000)

	# Define each little delta value by index position
	
			d45 = lil_del[60:120]
			d46 = lil_del[120:180]
			d47 = lil_del[180:240]
			d48 = lil_del[240:300]
			d49 = lil_del[300:360]

			lil_del_dict = {'d45':d45, 'd46':d46,'d47':d47, 'd48':d48, 'd49':d49}		
			df_lil_del = pd.DataFrame(lil_del_dict)
			df_lil_del.to_csv(Path.cwd().parents[0] / 'raw_data' / 'single_analyses' / 'lil_d.csv')

			print('little deltas calculated')
			print('ready to plot')

			y = 'm44'
			df_zc_sam = df_zc.iloc[sam_idx]
			plt.scatter(df_zc.index, df_zc[y], color = 'blue', s = 15, label = 'ref')
			plt.scatter(df_zc_sam.index, df_zc_sam[y], color = 'red', s = 15, label = 'sam')
			title = y + ' ' + str(file_number) + samp_name
			plt.title(title)
			plt.ylabel(y)
			plt.xlabel('cycle')
			plt.legend()
			plt.savefig(Path.cwd().parents[0] / 'raw_data' / 'single_analyses' / (title + '.png'))

			plt.close()

			x = list(range(len(d47)))
			plt.scatter(x, d47, color = 'blue', s = 15)			
			title =  'd47 ' + str(file_number) + ' ' + samp_name
			plt.title(title)
			plt.ylabel('d47')
			plt.xlabel('cycle')
			plt.savefig(Path.cwd().parents[0] / 'raw_data' / 'single_analyses' / (title + '.png'))


