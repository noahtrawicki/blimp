# Script for processing D47 data from raw delta values (without any metadata)

#from blimp_supp import run_D47crunch
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

raw_deltas_file = "C:/Users/noaha/Documents/MIT/mass_spec/dol_calib/Clumpycrunch dolomites_muller_from_Stefano_NTA_processed_irrel_removed.csv"

xls = pd.ExcelFile(Path.cwd().parents[0] / 'params.xlsx')
df_anc = pd.read_excel(xls, 'Anchors', index_col = 'Anchor')
Nominal_D47 = df_anc.to_dict()['D47'] # Sets anchor values for D47crunch as dictionary {Anchor: value}
output_sep = '================='



def run_D47crunch(run_type, raw_deltas_file):
	''' 
	PURPOSE: Calculate D47, d13C, d18O, using Mathieu Daeron's 'D47_crunch' package (https://github.com/mdaeron/D47crunch)	
	INPUT: Fully corrected little deltas ('raw_deltas.csv')
	OUTPUT: repeatability (1SD) of all calculated measurements; also writes 'sessions.csv', 'analyses.csv', and 'samples.csv',
	  '''
	import D47crunch
	
	results_path = Path.cwd() / 'results'

	xls = pd.ExcelFile(Path.cwd().parents[0] / 'params.xlsx')
	df_anc = pd.read_excel(xls, 'Anchors', index_col = 'Anchor')
	Nominal_D47 = df_anc.to_dict()['D47'] # Sets anchor values for D47crunch as dictionary {Anchor: value}


	data = D47crunch.D47data()
	data.Nominal_D47 = Nominal_D47

	print('Anchors are ', data.Nominal_D47)	# list anchors used and their nominal D47
	data.ALPHA_18O_ACID_REACTION = 1.00909 # This is selected from the params file -- you can use whatever value you want in there.

	#values from Bernasconi et al 2018 Table 4
	if run_type == 'standard':
		data.SAMPLE_CONSTRAINING_WG_COMPOSITION = ('ETH-1', 2.02, -2.19) # oftentimes for standard racks, we don't use ETH-3, so this uses ETH-1 as the anchor
	elif run_type == 'clumped':
		data.SAMPLE_CONSTRAINING_WG_COMPOSITION = ('ETH-3', 1.71, -1.78)

	print('Sample constraining WG composition = ', data.SAMPLE_CONSTRAINING_WG_COMPOSITION)
	#data.read(Path.cwd() / 'results' / raw_deltas_file) # INPUT
	data.read(raw_deltas_file) # INPUT	

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

		return df_sam, df_analy, repeatability_all



df_sam, df_analy, rptability = run_D47crunch('clumped', raw_deltas_file)
print('Repeatability for all samples is ', round(rptability, 3)*1000, 'ppm' )

df_meta = pd.read_excel('C:/Users/noaha/Documents/MIT/mass_spec/dol_calib/dolo_Muller_from_stefano_metadata.xlsx')

df_sam = df_sam.merge(df_meta, how = 'left', on = 'Sample')
df_analy = df_analy.merge(df_meta, how = 'left', on = 'Sample')

df_sam.to_csv(Path.cwd() / 'results' / 'samples_Muller_w_bad_analy.csv', index = False)
df_analy.to_csv(Path.cwd() / 'results' / 'analyses_Muller_w_bad_analy.csv', index = False)


from matplotlib.lines import Line2D

df = pd.read_csv(Path.cwd() / 'results' / 'analyses_Muller.csv')

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
			if df_anchor['D47'][j] > (Nominal_D47[i] + 3*rptability) or df_anchor['D47'][j] < (Nominal_D47[i] - 3*rptability):
				ax.scatter(df_anchor['D47'][j], df_anchor['Sample'][j], color = 'orange', alpha = 1, edgecolor = 'black', label = label)
				ax.text(df_anchor['D47'][j] + 0.005, df_anchor['Sample'][j], df_anchor['UID'][j], zorder = 6)
plt.xlabel('D47 I-CDES')
#plt.legend()
plt.tight_layout()
plt.savefig(Path.cwd().parents[0] / 'plots' / 'anchor_D47_Muller.png')
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

	# --- D47_all_interactive ----

	output_file(filename="D47_all_interactive.html", title="D47_all_interactive")

	sample_names = pd.unique(df['Sample'])
	TOOLTIPS = [("Sample name", "@Sample"), ("UID", "@UID")]
	std_tools = ['pan,wheel_zoom,box_zoom,reset,hover']
	f2 = figure(x_axis_label = 'D47', y_axis_label = 'Sample', y_range = sample_names, tools = std_tools, tooltips = TOOLTIPS)
	f2.circle(x = 'D47', y=jitter('Sample', width=0.1, range=f2.y_range), source = data_analyses, color = 'white', line_color = 'black', size = 7)

	save(f2)

	# --- D47_raw_nominal_interactive ---

	output_file(filename="D47_raw_nominal_interactive.html", title="D47_raw_nominal_interactive")
	TOOLTIPS = [("Sample name", "@Sample"), ("UID", "@UID")]
	std_tools = ['pan,wheel_zoom,box_zoom,reset,hover']

	df_anchor = df_anchor.reset_index(drop = True)
	f3 = figure(x_axis_label = 'UID', y_axis_label = 'D47raw-nominal', tools = std_tools, tooltips = TOOLTIPS)

	# for j in range(len(df_anchor)):
	# 	if df_anchor['Sample'][j] == 'ETH-1': col, nom_D47 = pal[0], Nominal_D47['ETH-1']
	# 	if df_anchor['Sample'][j] == 'ETH-2': col, nom_D47 = pal[1], Nominal_D47['ETH-2']
	# 	if df_anchor['Sample'][j] == 'ETH-3': col, nom_D47 = pal[2], Nominal_D47['ETH-3']
	# 	if df_anchor['Sample'][j] == 'ETH-4': col, nom_D47 = pal[3], Nominal_D47['ETH-4']
	# 	if df_anchor['Sample'][j] == 'IAEA-C2': col, nom_D47 = pal[4], Nominal_D47['IAEA-C2']
	# 	if df_anchor['Sample'][j] == 'MERCK': col, nom_D47 = pal[5], Nominal_D47['MERCK']

	f3.scatter('UID', 'D47', source = data_anchors, color = 'black')
	save(f3)

interactive_plots(df)

print('done')
