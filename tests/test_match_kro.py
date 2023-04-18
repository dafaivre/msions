from msions.percolator import match_kro
import msions.kronik as kro
import msions.percolator as perc
import msions.mzml as mzml  
import pandas as pd


def test_match_kro():
	"""Test match function for Kronik and Percolator XML output"""
	kro_df = kro.simple_df(pd.read_csv("tests/DDA_match.kro", 
										header=0, sep='\t')).sort_values(by="best_rt").reset_index(drop=True)
	perc_df = perc.psms2df("tests/psm_fixture.pout.xml")
	ms_df = mzml.tic_df("tests/large_fixtures/short_DDA_file.mzML", level="all", include_ms1_info=True)
	match_kro(kro_df, perc_df, ms_df)
	expected_matches = 7
	kro_matches = sum(kro_df.ID_d)
	perc_matches = sum(perc_df.in_kro)
	assert kro_matches == expected_matches, "Kronik and Percolator match function did not work properly. Please look at Kronik input"
	assert perc_matches == expected_matches, "Kronik and Percolator match function did not work properly. Please look at Percolator input"
	