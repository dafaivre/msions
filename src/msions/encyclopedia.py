"""
This module contains functions that are useful for interacting with
EncyclopeDIA files in Python.
"""
import sqlite3
import pandas as pd

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
