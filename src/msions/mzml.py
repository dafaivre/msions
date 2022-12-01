"""
This module contains functions that are useful for interacting with
mzML files in Python.
"""
import pymzml
import pandas as pd
from typing import List, Union


def tic_df(input_mzml: str) -> pd.DataFrame:
	"""
	Find the TIC and injection time for each scan in an mzML file.
	
	Parameters
	----------
	input_mzml : str
		The input mzML file.
		
	Returns
	-------
	pd.DataFrame
		A pandas DataFrame containing the retention time, TIC, and injection time for each scan.
	
	Examples
	-------
	>>> from msions.mzml import tic_df
	>>> test_tic_df = tic_df("test.mzML")
	"""
	# create Reader object
	run = pymzml.run.Reader(input_mzml)

	# create array
	ms1_tic = []

	# record scan time, TIC, & injection time
	for spectrum in run:
		if spectrum.ms_level == 1:
			ms1_tic.append([spectrum.ID, spectrum.scan_time[0], spectrum.TIC,
							spectrum.get_element_by_path(['scanList', 'scan', 'cvParam'])[2].get('value')])

	# create dataframe
	ms1_df = pd.DataFrame(ms1_tic,
						  columns=['scan_num', 'rt', 'TIC', 'IT'])

	ms1_df['rt'] = ms1_df['rt'].round(4)

	# update column data types
	ms1_df['IT'] = ms1_df['IT'].astype("float")

	# calculate ions per scan
	# ions per scan = ion current (for scan) * inject time /1000
	ms1_df["ions"] = ms1_df['TIC']*ms1_df['IT']/1000

	# return data frame
	return ms1_df
