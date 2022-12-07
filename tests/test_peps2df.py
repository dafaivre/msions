from msions.percolator import peps2df
from msions.percolator import parse_peps

def test_peps2df():
	"""Test Peptide DataFrame creation from list of dictionaries or XML file"""
	# test string input for XML file
	expected_type = "DataFrame"
	actual_df = peps2df("tests/pep_fixture.pout.xml")
	actual_type_str = type(actual_df).__name__
	assert actual_type_str == expected_type, "DataFrame was not created correctly. Check format of file."

	# test list of dictionaries input
	actual_type_dicts = type(peps2df(parse_peps("tests/pep_fixture.pout.xml"))).__name__
	assert actual_type_dicts == expected_type, "DataFrame was not created correctly. Check format of input."	
	
	# test expected rows and columns
	expected_rows = 141
	expected_columns = 6
	actual_rows = actual_df.shape[0]
	actual_columns = actual_df.shape[1]
	assert actual_rows == expected_rows, "Number of rows in DataFrame is unexpected."
	assert actual_columns == expected_columns, "Number of columns in DataFrame is unexpected."