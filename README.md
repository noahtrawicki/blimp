# blimp
Bergmann Lab Isotopic Measurement Pipeline for analysis of D47 data from a Nu Perspective IRMS.

This repository contains scripts to reduce raw voltage data from a Nu Perspective IRMS, calculate isotopic compositions using Mathieu Daeron's ['D47crunch' library](https://github.com/mdaeron/D47crunch), and create a suite of tables and figures that summarize data.

## Requirements
It is strongly recommended to run blimp within a conda environment using the provided 'environment.yml' file (instructions below). It's likely that you used Anaconda or Miniconda when you originally installed Python, but if not, install Miniconda here https://docs.conda.io/en/latest/miniconda.html.

If you don't use the conda environment, you will need pandas, numpy, pathlib, matplotlib, os, seaborn, and [D47crunch](https://github.com/mdaeron/D47crunch) (can be installed on the command line with 'pip install D47crunch').

Blimp has been tested on Nu Stable versions 1.69.5, 1.71.2, and 1.71.5. For now, blimp is only reliably usable for data from the Bergmann Lab Nu Perspective IRMS.

## Using blimp

1. Download and unzip blimp.
2. Copy data folder(s) created by the Nu Perspective software (containing raw data files and batch summary) and paste into the 'raw_data' folder. Ideally, all analyses therein will be associated with the project of interest.
3. If you would like to change the default parameters (nominal anchor values, analysis removal thresholds, a18O values, etc.), add metadata to your samples (e.g., sample mineralogy), or change the original names of analyses (e.g., 'MySmaple' --> 'MySample'), you should update the 'params.xlsx' worksheet.
4. Navigate to '.../blimp' using the Anaconda Command Prompt (use 'cd' command to change directories). 
5. Type 'conda env create -n blimp-env --file environment.yml'; this will create an 'environment' in which to run blimp.
6. Type 'conda activate blimp-env'.
7. Navigate to the '.../blimp/scripts' folder ('cd scripts').
8. Run blimp by typing 'python blimp_main.py' in your Anaconda Command Prompt. Blimp will take ~20 seconds -- 2 min to run, depending on how many analyses are processed and your machine.
9. Find processed data in 'results' folder and plots in 'plots' folder.

If you want a LATEX file/pdf with a summary of key results and plots, talk to Noah -- you will need to install (at least) pdflatex, Perl, and a LATEX interpreter.
