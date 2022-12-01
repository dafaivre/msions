"""
This module contains functions that are useful for interacting with
Hardklor output files in Python.
"""
import pandas as pd
from typing import Union


def hk2df(hk_file: str, by_int: bool = False) -> pd.DataFrame:
	"""
	Read a Hardklor tab-delimited file to a pandas DataFrame.
	
	After import, all columns that can be converted to a numeric data
	type will be.
	
	Parameters
	----------
	hk_file : str
		The Hardklor tab-delimited file to read.
	by_int: bool
		Sort data by intensity.
		
	Returns
	-------
	pd.DataFrame
		A pandas DataFrame of the input file.

	Examples
	-------
	>>> import msions.hardklor as hk
	>>> hk.hk2df("test.hk")	
	"""
	# open file
	with open(hk_file, "r") as open_file:
		scan_num = 0
		rt = 0.0
		pep_arrays = []

		# read file, keep scan number, retention time, and all peptide info
		for line in open_file:
			if line[0] == 'S':               
				scan_info = line.strip().split()
				scan_num = int(scan_info[1])
				rt = float(scan_info[2])        
			else:
				pep_info = line.strip().split()[1:]
				pep_info.extend([scan_num, rt])
				pep_arrays.append(pep_info)

		# create data frame from info
		pep_df = pd.DataFrame(pep_arrays,
							  columns=['mass', 'charge',
									   'intensity', 'base_peak',
									   'window', 'unk',
									  'mod', 'corr', 'scan_num', 'rt'])
		# change data types
		pep_df = pep_df.astype({'mass': 'float', 'charge': 'int64',
								'intensity': 'int64', 'base_peak': 'float',
							   'scan_num': 'int64', 'rt': 'float'})

		# sort by intensity if true
		if by_int:
			pep_df.sort_values(by="intensity", ascending=False, inplace=True)
			pep_df.reset_index(drop=True, inplace=True)
		
		# calculate m/z
		pep_df['mz'] = (pep_df['mass']+pep_df['charge']*1.00728)/pep_df['charge']

		# return data frame of info
		return pep_df


def summarize_df(hk_input: Union[pd.DataFrame, str], full_ms1_df: pd.DataFrame = None) -> pd.DataFrame:
	"""
	Summarize the TIC in each scan from a Hardklor pandas DataFrame or Hardklor tab-delimited file.
	
	If an additional pandas DataFrame is provided with the MS1 scan information,
	the ion injection time will be mapped to each scan.
	
	Parameters
	----------
	hk_input : pd.Dataframe or str
		The Hardklor pandas DataFrame or Hardklor tab-delimited file.
	full_ms1_df: pd.DataFrame
		The pandas DataFrame containing the MS1 scan information.
		
	Returns
	-------
	pd.DataFrame
		A summarized pandas DataFrame.

	Examples
	-------
	>>> import msions.hardklor as hk
	>>> hk_df = hk.hk2df("test.hk")
	>>> hk.summarize_df(hk_df)	
	"""
	# if it's a hardklor file
	if isinstance(hk_input, str):
		# create hardklor data frame
		hk_df = hk2df(hk_input)
		
	# if it's a data frame already
	else:
		hk_df = hk_input

	# group data frame by scan number & rt
	byscan_rt = hk_df.groupby(["scan_num", "rt"])
	
	# sum by grouping
	sum_group = pd.DataFrame(byscan_rt["intensity"].aggregate(sum))
	
	# rename column
	sum_group.rename(columns={'intensity': 'TIC'}, inplace=True)

	# reset indices
	sum_group.reset_index(level=0, inplace=True)
	sum_group.reset_index(level=0, inplace=True)

	# if complete data frame is given
	if full_ms1_df is not None:
		# merge ion injection time into data frame
		sum_group = pd.merge(sum_group, full_ms1_df[['scan_num', 'IT']], on='scan_num')

		# calculate ions per scan
		# ions per scan = ion current (for scan) * inject time /1000
		sum_group["ions"] = sum_group['TIC']*sum_group['IT']/1000

	# return data frame
	return sum_group
