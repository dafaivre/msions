"""
This module contains functions that are useful for interacting with
EncyclopeDIA files in Python.
"""
import sqlite3
import pandas as pd
import numpy as np


def dia_df(input_elib: str) -> pd.DataFrame:
	"""
	Create a pandas DataFrame from an EncyclopeDIA elib output
	
	Parameters
	----------
	input_elib : str
		The input elib file.
		
	Returns
	-------
	pd.DataFrame
		A pandas DataFrame containing PrecursorMz, PrecursorCharge, PeptideModSeq, PeptideSeq, 
		RtInSeconds, RTInSecondsStart, and RTInSecondsStop columns.
	
	Examples
	-------
	>>> from msions.encyclopedia import dia_df
	>>> dia_df("test.elib")
	"""
	# create connection object
	elib_connection = sqlite3.connect(input_elib)

	# create DataFrame with SQL query
	encyclo_df = pd.read_sql_query("""SELECT PrecursorMz, PrecursorCharge, PeptideModSeq, PeptideSeq, RtInSeconds, RTInSecondsStart, RTInSecondsStop
									FROM entries""", elib_connection)

	# close connection
	elib_connection.close()

	# return data frame
	return encyclo_df


def match_hk(ref_row, other_df): 
	""" 
	Match EncyclopeDIA elib output to Hardklor output

	Parameters 
	---------- 
	ref_row : str 
		The row of data to match.
	other_df : pd.DataFrame
		The other DataFrame to match.

	Returns 
	------- 
	pd.DataFrame 
		A pandas DataFrame containing PrecursorMz, PrecursorCharge, PeptideModSeq, PeptideSeq,  
		RtInSeconds, RTInSecondsStart, and RTInSecondsStop columns. 

	Examples 
	------- 
	>>> from msions.encyclopedia import match_hk 
	>>> import msions.hardklor as hk 
	>>> import msions.encyclopedia as encyclo
	>>> hk_df = hk.hk2df("test.hk") 
	>>> encyclo_df = encyclo.dia_df("test.elib") 
	>>> hk_df["in_encyclo"] = hk_df.apply(match_hk, axis=1, other_df=encyclo_df) 
	""" 
	# define info to match
	mz2match = ref_row.mz 
	charge2match = ref_row.charge 
	rt2match = ref_row.rt_s 

    # re-assign data frame to use previously written code
	small_df = other_df 

    # only search scans that match
	small_df = small_df.loc[(small_df.RTInSecondsStart <= rt2match) & (small_df.RTInSecondsStop >= rt2match) & (small_df.PrecursorCharge == charge2match)] 

    # look for mass that matches
	small_df = small_df[np.isclose(small_df.PrecursorMz, mz2match, rtol=5e-6)] 

	# return number of matches
	return(small_df.shape[0])
