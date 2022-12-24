"""
This module contains utility functions that are useful for interacting with
MS data.
"""
import numpy as np
import pandas as pd
import math
from typing import List

def bin_list(start: float, end: float, bin_size: float, bin_mult: float = 1) -> List[float]:
	"""
	Create a list of bin edges.
	
	Parameters
	----------
	start : float
		The starting bin value.
	end : float
		The ending bin value.
	bin_size : float
		The size of the bins.		
	bin_mult : float
		A multiplier to adjust bin sizing.

	Returns
	-------
	List[float]
		A list of bin edges.
	
	Examples
	-------
	>>> from msions.utils import bin_list
	>>> mz_bin_size = 4
	>>> mz_bin_mult = 1.0005
	>>> mz_start = 399
	>>> mz_end = 1005
	>>> bin_mz_list = bin_list(mz_start, mz_end, mz_bin_size, mz_bin_mult)
	"""    
	bin_list = []

	for num in np.arange(start, math.ceil(end+(bin_mult-1)*end+1), bin_size*bin_mult):
		bin_list.append(num)

	return bin_list

def bin_data(df: pd.DataFrame, type: str, bin_rt_list: List[float] = None, bin_mz_list: List[float] = None) -> pd.DataFrame:
	"""
	Bin a pandas DataFrame using list(s) of bin edges.
	
	Parameters
	----------
	df : pd.DataFrame
		The pandas DataFrame of data.
	type : str
		Type of binning ("rt", "mz", "both").
	bin_rt_list : List[float]
		List of retention time bin edges.		
	bin_mz_list : List[float]
		List of m/z bin edges.		

	Returns
	-------
	pd.DataFrame
		The binned pandas DataFrame.
	
	Examples
	-------
	>>> from msions.utils import bin_list
	>>> from msions.utils import bin_data
	>>> mz_bin_size = 4
	>>> mz_bin_mult = 1.0005
	>>> mz_start = 399
	>>> mz_end = 1005
	>>> bin_mz_list = bin_list(mz_start, mz_end, mz_bin_size, mz_bin_mult)
	>>> bin_data(peak_df, type="mz", bin_mz_list=bin_mz_list)
	"""    
	if type == "rt":
		# create bin column
		df['bin_rt'] = pd.cut(df.rt, bin_rt_list, right=False)
		# sum intensities into bins
		df_binned = df.groupby(['mz','bin_rt'], as_index=False)[['ips']].sum()       

	elif type == "mz":
		# create bin column
		df['bin_mz'] = pd.cut(df.mz, bin_mz_list, right=False)
		# sum intensities into bins
		df_binned = df.groupby(['rt','bin_mz'], as_index=False)[['ips']].sum()

	elif type == "both":
		# create bin columns
		df['bin_rt'] = pd.cut(df.rt, bin_rt_list, right=False)
		df['bin_mz'] = pd.cut(df.mz, bin_mz_list, right=False)
		# sum intensities into bins
		df_binned = df.groupby(['bin_rt','bin_mz'], as_index=False)[['ips']].sum()

	return df_binned