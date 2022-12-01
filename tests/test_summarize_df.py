from msions.hardklor import summarize_df
from msions.hardklor import hk2df
from msions.mzml import tic_df

def test_summarize_df():
	"""Test summarized DataFrame"""
	expected_rows = 2  ## testing function without MS1 DataFrame
	expected_columns = 5  ## testing function with MS1 DataFrame
	hk_df = hk2df("tests/hk_fixture.hk")
	actual_rows = summarize_df(hk_df).shape[0]  ## testing function without MS1 DataFrame
	actual_columns = summarize_df(hk_df, tic_df("tests/mzml_fixture.mzML")).shape[1]  ## testing function with MS1 DataFrame
	assert actual_rows == expected_rows, "DataFrame was not summarized correctly."
	assert actual_columns == expected_columns, "Ion injection times were not added properly."
