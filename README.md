# blimp
Bergmann Lab Isotopic Measurement Pipeline for analysis of D47 data from a Nu Perspective IRMS.

This repository contains scripts to reduce raw voltage data from a Nu Perspective IRMS, calculate isotopic compositions using Mathieu Daeron's ['D47crunch' library](https://github.com/mdaeron/D47crunch), and create a suite of tables and figures that summarize data.

## Requirements
The Anaconda distribution of Python 3 will have all the packages you will need (pandas, numpy, pathlib, matplotlib, os, seaborn **except** [D47crunch](https://github.com/mdaeron/D47crunch) (can be installed with 'pip install D47crunch'), and, optionally, [joypy](https://github.com/sbebo/joypy) ('pip install joypy').

Blimp has been tested on Nu Stable versions 1.69.5, 1.71.2, and 1.71.5. For now, blimp is only reliably usable for data from the Bergmann Lab Nu Perspective IRMS.

## Using blimp
1. Drag and drop data folders output by the Nu Perspective software (containing raw data files and batch summary) into 'raw_data' folder.
2. If you would like to change the default parameters (nominal anchor values, analysis removal thresholds, a18O values, etc.), add metadata to your samples (e.g., sample mineralogy), or change the original names of analyses (e.g., 'MySmaple' --> 'MySample'), you should update the corresponding spreadsheet ('params.xlsx'; 'metadata.xlsx'; 'names_to_change.csv').
3. Navigate to '.../blimp_test/scripts' and run 'blimp_main.py' in Terminal (Mac), Windows Command Prompt (Windows), or Anaconda Command Prompt.
4. Find processed data in 'results' folder and plots in 'plots' folder.


