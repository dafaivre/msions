from msions.kronik import simple_df

def test_simple_df():
	"""Test simplified DataFrame creation from a Kronik file or DataFrame"""
	expected_type = "DataFrame"
	actual_type = type(simple_df("tests/kro_fixture.kro")).__name__
	assert actual_type == expected_type, "DataFrame was not created correctly. Check format of file."