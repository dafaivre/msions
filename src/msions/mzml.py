"""
This module contains functions that are useful for interacting with
mzML files in Python.
"""
import pymzml
import pandas as pd
from typing import List, Union


def tic_df(input_mzml: str, level="1") -> pd.DataFrame:
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
	tic_lst = []

    #record scan, scan time, TIC, & injection time
	if level == "1":
		for spectrum in run:
			if spectrum.ms_level == 1:
				tic_lst.append([spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')])
	elif level == "2":
		for spectrum in run:
			if spectrum.ms_level == 2:
				tic_lst.append([spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')])
	elif level == "all":
		for spectrum in run:
			tic_lst.append([spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
							spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')])    

	# create dataframe
	tic_df = pd.DataFrame(tic_lst,
						  columns=['scan_num', 'rt', 'TIC', 'IT'])

	tic_df['rt'] = tic_df['rt'].round(4)

	# update column data types
	tic_df['IT'] = tic_df['IT'].astype("float")

	# calculate ions per scan
	# ions per scan = ion current (for scan) * inject time /1000
	tic_df["ions"] = tic_df['TIC']*tic_df['IT']/1000

	# return data frame
	return tic_df
