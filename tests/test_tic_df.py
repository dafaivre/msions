from msions.mzml import tic_df

def test_tic_df():
	"""Test DataFrame creation from an mzML file"""
	expected_type = "DataFrame"
	actual_type = type(tic_df("tests/mzml_fixture.mzML")).__name__
	assert actual_type == expected_type, "DataFrame was not created correctly. Check format of file."

