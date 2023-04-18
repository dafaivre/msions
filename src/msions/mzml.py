"""
This module contains functions that are useful for interacting with
mzML files in Python.
"""
import pymzml
import pandas as pd
from typing import List, Union


def tic_df(input_mzml: str, level: str = "1", include_ms1_info: bool = False, faims: bool = False) -> pd.DataFrame:
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
	faims : bool
		Returns CV associated with each scan.
		
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

    # record scan, scan time, TIC, & injection time
	# if examining MS1
	if level == "1":
		# for each spectrum
		for spectrum in run:
			info_lst = []
			if spectrum.ms_level == 1:
				try:
					if spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('name') == "ion injection time":
						it_val = spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')
						info_lst = [spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
									it_val]
					elif spectrum.get_element_by_path(['scanList','scan','cvParam'])[3].get('name') == "ion injection time":
						it_val = spectrum.get_element_by_path(['scanList','scan','cvParam'])[3].get('value')
						info_lst = [spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
									it_val]
				except:	
					"File has different formatting than expected. Ion injection time needs to be re-defined."
				
				if faims:
					try:
						if spectrum.get_element_by_path(['scanList','scan','cvParam'])[7].get('name') == "FAIMS compensation voltage":
							cv_val = spectrum.get_element_by_path(['scanList','scan','cvParam'])[7].get('value')
							info_lst.append(int(float(cv_val)))
					except:
						"File has different formatting than expected. CV needs to be re-defined."
				tic_lst.append(info_lst)
	# if examining MS2
	elif level == "2":
		# for each spectrum
		for spectrum in run:
			info_lst = []
			if spectrum.ms_level == 2:
				if include_ms1_info:
					# define triggered MS1 scan number
					info_lst = [int(spectrum.get_element_by_path(['precursorList', 'precursor'])[0].get('spectrumRef').split('=')[3]),
									# define triggered MS1 m/z
									float(spectrum.get_element_by_path(['precursorList', 'precursor', 
																		'selectedIonList', 'selectedIon', 'cvParam'])[0].get('value')),
									# define triggered MS1 intensity
									float(spectrum.get_element_by_path(['precursorList', 'precursor', 
																		'selectedIonList', 'selectedIon', 'cvParam'])[2].get('value')),
									spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
									# define injection time (IT)
									spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')]
				else:
					info_lst = [spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')]
				if faims:
					info_lst.append(int(float(spectrum.get_element_by_path(['cvParam'])[7].get('value'))))
				tic_lst.append(info_lst)
	elif level == "all":
		for spectrum in run:
			info_lst = []
			if spectrum.ms_level == 2:
				if include_ms1_info:
					# define triggered MS1 scan number
					info_lst = [int(spectrum.get_element_by_path(['precursorList', 'precursor'])[0].get('spectrumRef').split('=')[3]),
								# define triggered MS1 m/z
								float(spectrum.get_element_by_path(['precursorList', 'precursor', 
																	'selectedIonList', 'selectedIon', 'cvParam'])[0].get('value')),
								# define triggered MS1 intensity
								float(spectrum.get_element_by_path(['precursorList', 'precursor', 
																	'selectedIonList', 'selectedIon', 'cvParam'])[2].get('value')),
								spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								# define injection time (IT)
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')]
				else:
					info_lst = [spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')]
				if faims:
					info_lst.append(int(float(spectrum.get_element_by_path(['cvParam'])[7].get('value'))))
			else:
				if include_ms1_info:
					info_lst = [-1,
								-1,
								-1,
								spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')]
				else:
					info_lst = [spectrum.ID, spectrum.scan_time[0], spectrum.TIC, 
								spectrum.get_element_by_path(['scanList','scan','cvParam'])[2].get('value')]				
				if faims:
					info_lst.append(int(float(spectrum.get_element_by_path(['cvParam'])[7].get('value'))))
			tic_lst.append(info_lst)

	# create dataframe
	if (level == "2" or level == "all") and include_ms1_info:
		if faims:
			tic_df = pd.DataFrame(tic_lst,
								  columns=['ms1_scan', 'ms1_mz', 'ms1_int','scan_num', 'rt', 'TIC', 'IT', 'CV'])
			# round ms1_mz
			tic_df['ms1_mz'] = tic_df['ms1_mz'].round(4)
		else:
			tic_df = pd.DataFrame(tic_lst, 
								  columns=['ms1_scan', 'ms1_mz', 'ms1_int','scan_num', 'rt', 'TIC', 'IT'])
			# round ms1_mz
			tic_df['ms1_mz'] = tic_df['ms1_mz'].round(4)
	else:
		if faims:
			tic_df = pd.DataFrame(tic_lst,
								  columns=['scan_num', 'rt', 'TIC', 'IT', 'CV'])	
		else:	
			tic_df = pd.DataFrame(tic_lst,
								  columns=['scan_num', 'rt', 'TIC', 'IT'])

	# round retention time
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