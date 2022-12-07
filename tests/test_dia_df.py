from msions.encyclopedia import dia_df

def test_dia_df():
	"""Test DataFrame creation from an elib file"""
	expected_type = "DataFrame"
	actual_type = type(dia_df("tests/large_fixtures/elib_fixture.elib")).__name__
	assert actual_type == expected_type, "DataFrame was not created correctly. Check format of file."

