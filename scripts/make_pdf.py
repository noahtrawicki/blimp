def run_mk_pdf():
	import pdflatex
	from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
	    Plot, Figure, Matrix, Alignat, Command, Center, SmallText
	from pylatex.utils import italic, NoEscape, bold
	import os
	import pandas as pd
	from pathlib import Path


	df = pd.read_csv(Path.cwd().parents[0] / 'results' / 'project_info.csv', encoding = 'latin1')
	df_summ = pd.read_csv(Path.cwd().parents[0] / 'results' / 'summary.csv', encoding = 'latin1')
	doc = Document()

	#df_analy = pd.read_csv('C:/users/noaha/Desktop/blimp_test/results/analyses.csv', encoding = 'latin1')
	geometry_options = {'tmargin': '1in', 'lmargin': '1in', 'rmargin': '1in'}
	doc = Document(geometry_options=geometry_options)

	n_sess = df['n_sessions'][0]
	n_samp = df['n_samples'][0]
	n_anal = df['n_analyses'][0]


	doc.preamble.append(Command('title', 'Bergmann Lab clumped isotope analysis report'))
	doc.preamble.append(Command('author','MIT Carbonate Research Laboratory'))
	doc.append(NoEscape(r'\maketitle'))


	with doc.create(Center()) as centered:
		with centered.create(Tabular('r l', pos = 'c', booktabs = True)) as table:
		    #table.add_hline()
		    table.add_row((bold('Number of sessions: '), n_sess))
		    #table.add_hline(1, 2)
		    table.add_row((bold('Number of samples: '), n_samp))
		    table.add_row((bold('Number of analyses: '), n_anal))
		    #table.add_hline()
		    table.add_row((bold(NoEscape(r'Nominal $\Delta_{47}$ ETH-1: ')), df['Nominal_D47_ETH-1'][0]))
		    table.add_row((bold(NoEscape(r'Nominal $\Delta_{47}$ ETH-2: ')), df['Nominal_D47_ETH-2'][0]))
		    table.add_row((bold(NoEscape(r'Nominal $\Delta_{47}$ ETH-3: ')), df['Nominal_D47_ETH-3'][0]))
		    table.add_row((bold(NoEscape(r'Nominal $\Delta_{47}$ ETH-4: ')), df['Nominal_D47_ETH-4'][0]))
		    table.add_row((bold(NoEscape(r'Nominal $\Delta_{47}$ IAEA-C2: ')), df['Nominal_D47_IAEA-C2'][0]))
		    table.add_row((bold(NoEscape(r'Nominal $\Delta_{47}$ MERCK: ')), df['Nominal_D47_MERCK'][0]))
		    #table.add_row((bold('Measured D47 IAEA-C1: '), n_samp))
		    #table.add_hline()
		    table.add_row((bold(NoEscape(r'External reproducibility ($\delta^{13}C_{VPDB}$): ')), round(df['Reprod_d13C'][0],4)))
		    table.add_row((bold('External reproducibility (d18O_VSMOW): '), round(df['Reprod_d18O'][0],4)))
		    table.add_row((bold(NoEscape(r'External reproducibility ($\Delta_{47}$): ')), round(df['Reprod_D47'][0],4)))
		    #table.add_hline()
	    
	anchor_D47_path =  Path.cwd().parents[0] / 'plots' / 'anchor_D47.png'
	with doc.create(Figure(position='hbt!')) as plot:
		plot.add_image(str(anchor_D47_path), width = '320px')
		plot.add_caption(NoEscape(r'Final $\Delta_{47}$ of anchors; green diamonds = nominal D47, red = 2SD outside nominal.'))
	
	doc.append(NoEscape(r'\begin{footnotesize}'))
	with doc.create(Center()) as centered:
	    with centered.create(Tabular('r | l | l | l | l | l | l | l', pos = 'c', booktabs = True, row_height = 0.5)) as data_table:
	        header_row1 = ["Sample", "N", NoEscape(r"$\delta^{13}C_{VPDB}$"), NoEscape(r"Min. $\delta^{18}O_{VPDB}$"), 
	        NoEscape(r"$\Delta_{47}$"), "95% CL", NoEscape(r"T($^{\circ}$C; MIT)"), NoEscape(r"Fluid $\delta^{18}O$")]

	        data_table.add_row(header_row1, mapper=[bold])
	        data_table.add_hline()
	        #data_table.end_table_header()
	        
	        for i in range(len(df_summ)):
	            data_table.add_row(df_summ['Sample'][i], df_summ['N'][i], df_summ['d13C_VPDB'][i],
	                df_summ['d18O_VPDB_mineral'][i], round(df_summ['D47'][i],4), df_summ['95% CL'][i],
	                df_summ['T_MIT'][i], df_summ['d18O_water_VSMOW'][i])
	            data_table.add_hline(color = 'gray')
	doc.append(NoEscape(r'\end{footnotesize}'))

	        #data_table.add_caption(NoEscape(r'Summary of isotopic results. $\Delta_{47}$ reported in the I-CDES 90$^{\circ}$C reference frame.'))

	#with doc.create(Center()) as centered:	

	cdv_path = Path.cwd().parents[0] / 'plots' / 'cdv.png'
	with doc.create(Figure(position='ht!')) as plot:
	    plot.add_image(str(cdv_path), width = '500px')
	    plot.add_caption('Gas prep line values.')

	joyplot_path = Path.cwd().parents[0] / 'plots' / 'ETH-3_joyplot.png'
	with doc.create(Figure(position='ht!')) as plot:
	    plot.add_image(str(joyplot_path), width = '300px')
	    plot.add_caption('Distribution of cycle d47 for ETH-3 over project.')

	c1_path = Path.cwd().parents[0] / 'plots' / 'IAEA-C1_reprod.png'
	with doc.create(Figure(position='ht!')) as plot:
	    plot.add_image(str(c1_path), width = '300px')
	    plot.add_caption('Final D47 I-CDES value for IAEA-C1 (expected D47 = 0.30)')

	allD47_path = Path.cwd().parents[0] / 'plots' / 'all_D47.png'
	with doc.create(Figure(position='ht!')) as plot:
	    plot.add_image(str(allD47_path), width = '300px')
	    plot.add_caption('Final D47 I-CDES for all analyses.')
	    

	with doc.create(Section('Methods')):
		with doc.create(Subsection('Mass spectrometry')):
			doc.append('These are the mass spectrometry methods.')
		with doc.create(Subsection('Data processing')):
			doc.append('These are the data processing methods.')
	    

	doc.generate_pdf(Path.cwd().parents[0] / 'results' / 'summary', clean_tex = False,compiler='pdfLaTeX')
