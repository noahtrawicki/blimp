import os
import pandas as pd
from pathlib import Path

# Directory where 'params' file is stored
def define_params_location():
	params_loc = input("Enter path to params file from %s (omit first /):\n" %Path.cwd())
	params_dir = Path.cwd() / params_loc
	return params_dir

def check_params_location():
	params_dir = define_params_location() # user inputs directory containing params file
	# check that params file exists:
	# first print the directory that the user has directed us to:
	params_yn = input('Is the params file in %s? (Y/N):\n' %params_dir)

	return params_yn, params_dir

params_yn = 'N'

while params_yn == 'N':
	params_yn, params_dir = check_params_location()

if Path.is_file(params_dir / 'params.xlsx'):
	params = params_dir / 'params.xlsx'
else:
	print('params.xlsx not found in %s' %params_dir)

print(params)
