"""
This module contains functions that are useful for interacting with
Kronik output files in Python.
"""
import pandas as pd
from typing import Union


def simple_df(kro_input: Union[pd.DataFrame, str], topN: int = None, bestInt_thresh: float = None,
			  sumInt_thresh: float = None, remove1: bool = False, by_int: bool = False) -> pd.DataFrame:
	"""
	Create a simplified Kronik pandas DataFrame.
	
	The DataFrame can be filtered by topN intensity values and/or by removing +1 charges.
	
	Parameters
	----------
	kro_input : pd.Dataframe or str
		The Kronik pandas DataFrame or Kronik tab-delimited file.
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
		
	# if it's a data frame already
	else:
		kro_df = kro_input
	
	# select columns of interest
	df_short = kro_df.loc[:, ["Monoisotopic Mass", "Charge", "Best Intensity", "Summed Intensity",
								 "First RTime", "Last RTime", "Best RTime"]]

	# rename to remove spaces
	df_short.rename(columns={'Monoisotopic Mass': 'mass', 'Charge': 'charge',
							 'Best Intensity': 'best_int', 'Summed Intensity': 'sum_int',
							 'First RTime': 'first_rt', 'Last RTime': 'last_rt',
							 'Best RTime': 'best_rt'}, inplace=True)

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
