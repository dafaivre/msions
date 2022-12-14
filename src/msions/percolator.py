"""
This module contains functions that are useful for interacting with
XMLs, such as Percolator output.
"""
import xml.etree.ElementTree as ET
import pandas as pd
from typing import Union, List
import numpy as np


def parse_psms(xmlfile: str) -> List[dict]:
	"""
	Parse the PSMs in an XML file.
	
	Parameters
	----------
	xmlfile : str
		The XML file.
		
	Returns
	-------
	List[dict]
		A list of dictionaries containing PSM information.

	Examples
	-------
	>>> from msions.percolator import parse_psms
	>>> parse_psms("test.xml")
	"""
	# create element tree object
	tree = ET.parse(xmlfile)

	# get root element
	root = tree.getroot()

	# create empty list for PSMs
	psms = []

	# define root string
	root_string = './{http://per-colator.com/percolator_out/15}psms/{http://per-colator.com/percolator_out/15}psm'

	# iterate PSMs
	for psm in root.findall(root_string):

		# empty PSM dictionary
		psm_dict = {}

		psm_dict['{http://per-colator.com/percolator_out/15}psm_id'] = psm.attrib['{http://per-colator.com/percolator_out/15}psm_id']

		# iterate child elements of PSM
		for child in psm:

			# record PSM information in dictionary
			if child.tag == '{http://per-colator.com/percolator_out/15}protein_id':
				if '{http://per-colator.com/percolator_out/15}protein_id' in psm_dict.keys():
					psm_dict['{http://per-colator.com/percolator_out/15}protein_id'].append(child.text)
					# might need .encode(utf8)?
				else:
					psm_dict['{http://per-colator.com/percolator_out/15}protein_id'] = [child.text]
					# might need encode?
			elif child.tag == '{http://per-colator.com/percolator_out/15}peptide_seq':
				psm_dict['{http://per-colator.com/percolator_out/15}peptide_seq'] = (child.attrib['seq'])
			else:
				psm_dict[child.tag] = child.text
				# might need encode?

		# append PSM dictionary to PSM list
		psms.append(psm_dict)

	# return PSM list
	return psms


def parse_peps(xmlfile: str) -> List[dict]:
	"""
	Parse the peptides in an XML file.

	Parameters
	----------
	xmlfile : str
		The XML file.
		
	Returns
	-------
	List[dict]
		A list of dictionaries containing peptide information.

	Examples
	-------
	>>> from msions.msxml import parse_peps
	>>> parse_peps("test.xml")
	"""
	# create element tree object
	tree = ET.parse(xmlfile)

	# get root element
	root = tree.getroot()

	# create empty list for peptides
	peptides = []

	root_string = './{http://per-colator.com/percolator_out/15}peptides/{http://per-colator.com/percolator_out/15}peptide'

	# iterate peptides
	for peptide in root.findall(root_string):

		# empty peptide dictionary
		pep = {}

		pep['{http://per-colator.com/percolator_out/15}peptide_id'] = peptide.attrib['{http://per-colator.com/percolator_out/15}peptide_id']

		# iterate child elements of peptide
		for child in peptide:

			# record peptide information in dictionary
			if child.tag == '{http://per-colator.com/percolator_out/15}psm_ids':
				for grand_child in child:
					if '{http://per-colator.com/percolator_out/15}psm_ids' in pep.keys():
						pep['{http://per-colator.com/percolator_out/15}psm_ids'].append(grand_child.text)
						# might need .encode(utf8)?
					else:
						pep['{http://per-colator.com/percolator_out/15}psm_ids'] = [grand_child.text]
						# might need encode?
			else:
				pep[child.tag] = child.text
				# might need encode?

		# append peptide dictionary to peptides list
		peptides.append(pep)

	# return peptides list
	return peptides


def psms2df(xml_input: Union[List[dict], str]) -> pd.DataFrame:
	"""
	Create a pandas DataFrame of PSM XML information.

	Parameters
	----------
	xml_input : list[dict] or str
		The PSM list of dictionaries or the XML file.
		
	Returns
	-------
	pd.DataFrame
		The pandas DataFrame of PSM information.

	Examples
	-------
	>>> from msions.percolator import psms2df
	>>> psms2df("test.xml")
	""" 
	# if it's an XML file
	if isinstance(xml_input, str):
		# create list of dictionaries
		psm_xml = parse_psms(xml_input)

	# if it's a list of dictionaries already
	else:
		psm_xml = xml_input

	# define prefix
	prefix = '{http://per-colator.com/percolator_out/15}'

	# initiate array
	xml_lst = []

	# iterate through peptide xml
	for psm in psm_xml:
		seq = psm[prefix+'peptide_seq']
		prots = ','.join(psm[prefix+'protein_id'])
		q_val = psm[prefix+'q_value']
		exp_mass = psm[prefix+'exp_mass']
		calc_mass = psm[prefix+'calc_mass']
		scan = psm[prefix+'psm_id'].strip().split('_')[2]
		xml_lst.append([seq, prots, q_val, exp_mass, calc_mass, scan])

	# create pandas DataFrame
	xml_df = pd.DataFrame(xml_lst, 
						  columns=['peptide', 'protein_s',
								   'q_value',
								   'exp_mass', 'calc_mass',
								   'scan_num'])

	# change data types
	xml_df = xml_df.astype({'q_value': 'float',
							'exp_mass': 'float',
							'calc_mass': 'float',
							'scan_num': 'int64'})

	# return data frame
	return xml_df	


def peps2df(xml_input: Union[List[dict], str]) -> pd.DataFrame:
	"""
	Create a pandas DataFrame of peptide XML information.

	Parameters
	----------
	xml_input : list[dict] or str
		The peptide list of dictionaries or the XML file.
		
	Returns
	-------
	pd.DataFrame
		The pandas DataFrame of peptide information.

	Examples
	-------
	>>> from msions.percolator import peps2df
	>>> peps2df("test.xml")
	""" 
	# if it's an XML file
	if isinstance(xml_input, str):
		# create list of dictionaries
		pep_xml = parse_peps(xml_input)

	# if it's a list of dictionaries already
	else:
		pep_xml = xml_input

	# define prefix
	prefix = '{http://per-colator.com/percolator_out/15}'

	# initiate array
	xml_lst = []

	# iterate through peptide xml
	for peptide in pep_xml:
		seq = peptide[prefix+'peptide_id']
		q_val = peptide[prefix+'q_value']
		exp_mass = peptide[prefix+'exp_mass']
		calc_mass = peptide[prefix+'calc_mass']
		prot = peptide[prefix+'protein_id']
		for psm in peptide[prefix+'psm_ids']:
			scan = psm.strip().split('_')[2]
			xml_lst.append([seq, q_val, exp_mass, calc_mass, prot, scan])

	# create pandas DataFrame
	xml_df = pd.DataFrame(xml_lst,
						  columns=['peptide', 'q_value',
								   'exp_mass', 'calc_mass',
								   'protein', 'scan_num'])

	# change data types
	xml_df = xml_df.astype({'q_value': 'float',
							'exp_mass': 'float',
							'calc_mass': 'float',
							'scan_num': 'int64'})

	# return data frame
	return xml_df


def id_scans(perc_target, ms2_tic_df):
	"""
	Create a column saying whether an MS2 was identified

	Parameters
	----------
	perc_target : str
		The TXT file of percolator output.
	ms2_tic_df: pd.DataFrame
		The pandas DataFrame of MS2 scan information.

	Examples
	-------
	>>> from msions.percolator import id_scans
	>>> ms2_tic_df = mzml.tic_df("test.mzML", level="2")
	>>> id_scans("test.percolator.target.peptides.txt", ms2_tic_df)
	""" 
	# create DataFrame of percolator results for run
	perc_df = pd.read_csv(perc_target, header=0, sep="\t")

	# remove q-values greater than 0.01 and sort by scan
	perc_sig = perc_df[perc_df['percolator q-value'] < 0.01].sort_values(by="scan").reset_index(drop=True)

    # create list for whether MS2 was ID'd
	ms2_tic_df["IDd"] = np.isin(ms2_tic_df["scan_num"], perc_sig['scan'])
