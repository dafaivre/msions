from msions.mzml import find_precursorscan
import msions.percolator as perc

def test_find_precursorscan():
	"""Test precursor scan function"""
	expected_data =[53645, 962.962768554688, 7377221.03125]
	psm_xml_df = perc.psms2df("tests/psm_fixture.pout.xml")
	run = "tests/large_fixtures/short_DDA_file.mzML"
	actual_data = psm_xml_df['scan_num'].apply(find_precursorscan, pymzml_input=run)[0]
	assert actual_data == expected_data, "Function did not find the correct precursor scans."