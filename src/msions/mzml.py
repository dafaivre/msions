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


def find_precursorscan(ms2_scannum: int, pymzml_input) -> List[Union[int, float, float]]: 
	""" 
	Find the MS1 precursor scan info for a PSM scan number.  
	This function can be applied to the PSM pandas DataFrame 
	'scan_num' column.
	
	Parameters 
	---------- 
	ms2_scannum : int 
		The MS2 scan number. 
	pymzml_input : pymzml.run.Reader(mzML_file) or str 
		The mzML Reader object or file. 

	Returns 
	------- 
	List[int, float, float] 
		List of precursor scan number, precursor m/z, and precursor intensity. 

	Examples 
	------- 
	>>> from msions.mzml import find_precursorscan
	>>> from msions.mzml import tic_df 
	>>> psm_xml_df = perc.psms2df("test.xml")
	>>> run = "test.mzML"
	>>> psm_xml_df['ms1_scan'], psm_xml_df['ms1_mz'], psm_xml_df['ms1_intensity'] =  
			zip(*psm_xml_df['scan_num'].apply(find_precursorscan, pymzml_input=run)) 
	""" 
	# check if run is provided as a file 
	if isinstance(pymzml_input, str): 
		# create reader object 
		pymzml_run = pymzml.run.Reader(pymzml_input) 

	# if already a reader object 
	else: 
		pymzml_run = pymzml_input 

	# define spectrum reader object 
	spectrum = pymzml_run[ms2_scannum] 

	# find MS1 precursor scan for MS2 
	ms1_precursorscan = int(spectrum.get_element_by_path(['precursorList', 'precursor'])[0].get('spectrumRef').split('=')[3]) 
	ms1_mz = float(spectrum.get_element_by_path(['precursorList', 'precursor', 
												 'selectedIonList', 'selectedIon', 'cvParam'])[0].get('value')) 
	ms1_intensity = float(spectrum.get_element_by_path(['precursorList', 'precursor', 
														'selectedIonList', 'selectedIon', 'cvParam'])[2].get('value')) 

	# return MS1 scan number, precursor m/z, and precursor intensity
	return [ms1_precursorscan, ms1_mz, ms1_intensity]
