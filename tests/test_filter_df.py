from msions.kronik import simple_df
from msions.kronik import filter_df

def test_filter_df():
	"""Filter simplified Kronik DataFrame"""
	kro_df = simple_df("tests/kro_fixture.kro")
	expected_rows = 5175
	actual_rows = filter_df(kro_df, start=30, stop=45).shape[0]
	assert actual_rows == expected_rows, "DataFrame was not filtered correctly."