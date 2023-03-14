"""
This module contains functions that are useful for interacting with
Kronik output files in Python.
"""
import pandas as pd
import numpy as np
from typing import Union


def simple_df(kro_input: Union[pd.DataFrame, str], cv: Union[int, str] = None, topN: int = None, bestInt_thresh: float = None,
			  sumInt_thresh: float = None, remove1: bool = False, by_int: bool = False) -> pd.DataFrame:
	"""
	Create a simplified Kronik pandas DataFrame.
	
	The DataFrame can be filtered by topN intensity values and/or by removing +1 charges.
	
	Parameters
	----------
	kro_input : pd.Dataframe or str
		The Kronik pandas DataFrame or Kronik tab-delimited file.
	cv : int or str
		CV value associated with the dataset or "given" for already present
	topN : int
		Only include features with topN summed intensity.
	bestInt_thresh: float
		Only include features with apex intensity above intensity threshold.
	sumInt_thresh: float
		Only include features with summed intensity above intensity threshold.
	remove1 : bool
		Remove +1 charges from DataFrame.
	by_int: bool
		Sort data by summed intensity.
		
	Returns
	-------
	pd.DataFrame
		A pandas DataFrame of the input file.

	Examples
	-------
	>>> import msions.kronik as kro
	>>> kro.simple_df("test.kro")	
	"""
	# if it's a kronik file
	if isinstance(kro_input, str):
		# create kronik data frame
		kro_df = pd.read_csv(kro_input, header=0, sep='\t')

		# if dataset has a FAIMS CV
		if cv is not None:
			# if cv is an integer
			if isinstance(cv, int):
				# define CV associated with scans
				kro_df['CV'] = cv

			# if cv is not equal to given
			else: 
				assert cv == "given", "CV is not an integer or 'given.' Please check input."

			# select columns of interest
			df_short = kro_df.loc[:, ["First Scan","Last Scan", "Num of Scans", "Monoisotopic Mass", 
									  "Charge", "Best Intensity", "Summed Intensity", "Best RTime", "CV"]]

		else:
			# select columns of interest
			df_short = kro_df.loc[:, ["First Scan","Last Scan", "Num of Scans", "Monoisotopic Mass", 
									  "Charge", "Best Intensity", "Summed Intensity", "Best RTime"]]
		
	# if it's a data frame already
	else:
		kro_df = kro_input

		# if dataset has a FAIMS CV
		if cv is not None:
			# select columns of interest
			df_short = kro_df.loc[:, ["First Scan","Last Scan", "Num of Scans", "Monoisotopic Mass", 
									  "Charge", "Best Intensity", "Summed Intensity", "Best RTime", "CV"]]

		else:
			# select columns of interest
			df_short = kro_df.loc[:, ["First Scan","Last Scan", "Num of Scans", "Monoisotopic Mass", 
									  "Charge", "Best Intensity", "Summed Intensity", "Best RTime"]]

	# rename to remove spaces
	df_short.rename(columns={'First Scan':'first_scan','Last Scan':'last_scan',
                             'Num of Scans':'num_scans',
                             'Monoisotopic Mass':'mass', 'Charge':'charge',
                             'Best Intensity':'best_int',
                             'Summed Intensity':'sum_int', 'Best RTime':'best_rt'}, inplace=True)

	# round retention time
	df_short['best_rt'] = df_short['best_rt'].round(4)

	# remove features with +1 charge
	if remove1:
		df_short = df_short[df_short.charge != 1]

		df_short.reset_index(drop=True, inplace=True)

	# sort DataFrame by summed intensity of features
	if by_int:
		df_short.sort_values(by="sum_int", ascending=False, inplace=True)
		
		df_short.reset_index(drop=True, inplace=True)

	# filter DataFrame to only include topN features
	if topN is not None:
		if by_int:
			df_short = df_short.iloc[0:topN, ]
		else:
			df_short.sort_values(by="sum_int", ascending=False, inplace=True)
		
			df_short.reset_index(drop=True, inplace=True)
		
			df_short = df_short.iloc[0:topN, ]

			df_short.sort_values(by="best_rt", inplace=True)

			df_short.reset_index(drop=True, inplace=True)

	# filter DataFrame to only include features with apex intensity above threshold
	if bestInt_thresh is not None:
		df_short = df_short[df_short["best_int"] >= bestInt_thresh]

		df_short.reset_index(drop=True, inplace=True)
			
	# filter DataFrame to only include features with summed intensity above threshold
	if sumInt_thresh is not None:
		df_short = df_short[df_short["sum_int"] >= sumInt_thresh]

		df_short.reset_index(drop=True, inplace=True)

	# calculate m/z for each feature
	df_short['mz'] = (df_short['mass']+df_short['charge']*1.00728)/df_short['charge']

	# calculate retention time in seconds
	df_short['best_rt_s'] = df_short['best_rt']*60

	return df_short


def filter_df(df, start=0, stop=None) -> pd.DataFrame:
	"""
	Filter a pandas DataFrame containing Kronik data with a start and stop time.
	
	Parameters
	----------
	df : pd.DataFrame
		pandas DataFrame containing Kronik data.
	start : float
		Starting time to use to filter the DataFrame.
	stop : float
		Ending time to use to filter the DataFrame.

	Returns
	-------
	pd.DataFrame
		A filtered pandas DataFrame.

	Examples
	-------
	>>> import msions.kronik as kro
	>>> kro_df = kro.simple_df("test.kro")
	>>> kro.filter_df(kro_df, start=15.0)	
	"""
	# if there is not a stop time
	if stop == None:
		# make the stop time the last retention time
		stop = np.sort(df.best_rt)[-1]
    
	return df.loc[df.best_rt.between(start, stop)] #may need to add .copy() to prevent SettingwithCopyWarning


def match_rt_mass(ref_row: pd.Series, other_df: pd.DataFrame, rt_diff: float = None) -> int:
	""" 
	Match Kronik output with itself.

	Parameters 
	---------- 
	ref_row : str 
		The row of data to match.
	other_df : pd.DataFrame
		The other DataFrame to match.
	rt_diff : float
		Retention time difference window to use to search for a match.

	Returns 
	------- 
	int
		Number of matches in DataFrame.

	Examples 
	------- 
	>>> from msions.kronik import simple_df
	>>> from msions.kronik import match_rt_mass 
	>>> kro_df = simple_df("test.kro") 
	>>> redund_df = kro_df.copy()
	>>> redund_df["redund"] = redund_df.apply(match_rt_mass, axis=1, other_df=kro_df, rt_diff=1) 
	""" 
	# define info to match
	mass2match = ref_row.mass
	charge2match = ref_row.charge
	rt2match = ref_row.best_rt

	# if a retention time difference is given
	if rt_diff is not None:
		# filter DataFrame
		small_df = filter_df(other_df, rt2match-rt_diff, rt2match+rt_diff)
	else:
		small_df = other_df

    # only search scans that match	
	small_df = small_df.loc[small_df.charge == charge2match]

    # look for mass that matches
	small_df = small_df[np.isclose(small_df.mass, mass2match, rtol=5e-6)]

	# return number of rows (subtracting 1 for self-match)
	return small_df.shape[0] - 1
