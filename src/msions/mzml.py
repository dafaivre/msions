"""
This module contains functions that are useful for interacting with
mzML files in Python.
"""
import pymzml
import pandas as pd
from typing import List, Union


def tic_df(input_mzml: str, level: str = "1", include_ms1_info: bool = False) -> pd.DataFrame:
	"""
	Find the TIC and injection time for each scan in an mzML file.
	
	Parameters
	----------
	input_mzml : str
		The input mzML file.
	level : str
		Level of MS scan (1 or 2)
	include_ms1_info : bool
		Returns MS1 scan number, m/z, and intensity associated with precursor analyzed in MS2
		(requires level="2")
		
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
				if include_ms1_info:
					tic_lst.append([int(spectrum.get_element_by_path(['precursorList', 'precursor'])[0].get('spectrumRef').split('=')[3]),
									float(spectrum.get_element_by_path(['precursorList', 'precursor', 
																		'selectedIonList', 'selectedIon', 'cvParam'])[0].get('value')),
									float(spectrum.get_element_by_path(['precursorList', 'precursor', 
																		'selectedIonList', 'selectedIon', 'cvParam'])[2].get('value')),
									spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
									spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')])
				else:
					tic_lst.append([spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
									spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')])
	elif level == "all":
		for spectrum in run:
			tic_lst.append([spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
							spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')])    

	# create dataframe
	if level == "2" and include_ms1_info:
		tic_df = pd.DataFrame(tic_lst,
							columns=['ms1_scan', 'ms1_mz', 'ms1_int','ms2_scan', 'rt', 'TIC', 'IT'])
	else:
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


def peak_df(input_mzml: str) -> pd.DataFrame:
	""" 
	Create a pandas DataFrame containing the m/z, 
	ion current, and retention time for all MS1 peaks.
	
	Parameters
	----------
	input_mzml : str
		The input mzML file.
		
	Returns
	-------
	pd.DataFrame
		A pandas DataFrame containing the m/z, ion current, and retention time for all MS1 peaks.

	Examples 
	------- 
	>>> from msions.mzml import peak_df
	>>> peak_df("test.mzML")
	""" 
	# create run object
	run = pymzml.run.Reader(input_mzml)

	# initiate peak DataFrame
	peak_df = pd.DataFrame(columns=["mz", "ips", "rt"])

	# loop through spectra
	for spectra in run:
		if spectra.ms_level == 1:
			peak_array = pd.DataFrame(spectra.peaks("centroided")).rename(columns={0: "mz", 1: "ips"})
			peak_array["mz"] = peak_array["mz"].round(4)
			peak_array["rt"] = spectra.scan_time[0]
			peak_df = pd.concat([peak_df, peak_array])	

	return peak_df