from msions.mzml import peak_df

def test_peak_df():
	"""Test DataFrame creation from an mzML file"""
	expected_type = "DataFrame"
	expected_rows = 1329
	ms1_peaks = peak_df("tests/mzml_fixture.mzML")
	actual_type = type(ms1_peaks).__name__
	actual_rows = ms1_peaks.shape[0]
	assert (actual_type == expected_type) and (actual_rows == expected_rows), "DataFrame was not created correctly. Check format of file."

