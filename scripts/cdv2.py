# clumped data vetting V2 (cdv2)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import blimp_supp as b_s
from bokeh.io import output_file
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.layouts import row


dir_path = Path.cwd().parents[0]
os.chdir(dir_path)

rd_path = Path.cwd() / 'raw_data'
results_path = Path.cwd() / ('results')
plot_path = Path.cwd() / ('plots')

if os.path.isdir(results_path):
	pass 
else:
	os.mkdir(results_path)

if os.path.isdir(plot_path):
	pass 
else:
	os.mkdir(plot_path)

d47_crunch_fmt = []
batch_data_list = []
fold_count = 0

def plot_lil_d(file_n, samp_name, lil_d):
	pass

# -- THIS GRABS NU SUMMARY FILE and PLOTS -- 

def read_and_plot_summary(folder_name):

	this_path = Path.cwd() / 'raw_data' / folder_name
	for i in os.listdir(this_path):
		print(i)
		if 'Batch Results.csv' in i:# and 'fail' not in os.listdir(this_path): # checks for results summary file
			summ_file = Path.cwd() / 'raw_data' / folder_name / i
			df_summ = pd.read_csv(summ_file, encoding = 'latin1', skiprows = 3, header = [0,1])
			df_summ.columns = df_summ.columns.map('_'.join).str.strip()	# fixes weird headers of Nu Summary files
			print(df_summ.columns)
			df_summ['Data_File'] = df_summ['Data_File'].str.slice(start = 5, stop = 10)

			break


	# -- PLOTS --

	output_file('bokeh_TEST.html')

	df = pd.read_csv(Path.cwd() / 'results' / 'analyses.csv')

	data = ColumnDataSource.from_df(df_summ)
	data_analyses = ColumnDataSource.from_df(df)


	TOOLTIPS = [("Sample name", "@Sample_Name"),
				("UID", "@Data_File")]
	std_tools = ['pan,wheel_zoom,box_zoom,reset,hover']

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

	TOOLTIPS = [("Sample name", "@Sample"),
			("UID", "@UID")]



	f3 = figure(x_axis_label = 'd47',
				y_axis_label = 'D47',
				tools = std_tools,
				tooltips = TOOLTIPS)
	f3.scatter('d47', 'D47', source = data_analyses)

	show(row(f1, f3))

	exit()


		# -- THIS MAKES PLOTS --

# -- THIS CALCULATES LITTLE DELTAS FROM SCRATCH
if os.path.isdir(rd_path): # If there is a raw data folder...
	print('Crunching folders: ')
	for folder in os.listdir(rd_path):
		print(folder)	
		if 'clumped' in folder or 'Clumped' in folder:			
			# for file in os.listdir(Path.cwd() / 'raw_data' / folder):
			# 	if 'Data' in file and '.txt' in file and '.fail' not in file:
			# 		file_path = rd_path / folder / file					
			# 		if os.path.getsize(file_path) > 225000:	# Make sure .txt file is complete					
			# 			file_n = int(file[5:10])
			# 			samp_name = str(file[10:-4])													
			# 			lil_d = b_s.read_Nu_data(Path.cwd() / 'raw_data' / folder / file, file_n, samp_name)
			# 			plot_lil_d(file_n, samp_name, lil_d)					
			# 			d47_crunch_fmt.append(lil_d)

			read_and_plot_summary(folder)

			fold_count += 1