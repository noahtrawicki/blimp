
import os
import blimp_supp as b_s
import pandas as pd
from pathlib import Path
#import make_pdf as mk_pdf

output_sep = '--------------'

print(output_sep)

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

df_rmv = pd.read_excel('rmv.xlsx')
manual_rmv = list(df_rmv.UID)

if os.path.isdir(rd_path): # If there is a raw data folder...
	print('Crunching folders: ')
	for folder in os.listdir(rd_path):
		print(folder)	
		if 'clumped' in folder or 'Clumped' in folder:			
			for file in os.listdir(Path.cwd() / 'raw_data' / folder):
				if 'Data' in file and '.txt' in file and '.fail' not in file:
					file_path = rd_path / folder / file					
					if os.path.getsize(file_path) > 225000:	# Make sure .txt file is complete					
						file_n = int(file[5:10])
						samp_name = str(file[10:-4])
						if file_n not in manual_rmv:								
							lil_d, batch_data = b_s.raw_to_D47crunch_fmt(file_path, file_n, samp_name, fold_count, folder)
							
							if lil_d != None:		
								d47_crunch_fmt.append(lil_d)
							if batch_data != None:
								batch_data_list.append(batch_data)

			fold_count += 1


df_d47 = pd.DataFrame(d47_crunch_fmt, columns = ['UID', 'Session', 'Sample', 'd45', 'd46', 'd47', 'd48', 'd49'])
b_s.fix_names(df_d47)

rptability = b_s.run_D47crunch()


b_s.add_metadata(results_path, rptability, batch_data_list)

print(output_sep)
print('Repeatability for all samples is ', round(rptability, 3)*1000, 'ppm' )

b_s.plot_ETH_D47(rptability)

b_s.joy_plot()
b_s.cdv_plots()
#mk_pdf.run_mk_pdf()
print(output_sep)
print('SCRIPT COMPLETE.')
