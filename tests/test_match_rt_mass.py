from msions.kronik import simple_df
from msions.kronik import match_rt_mass

def test_filter_df():
	"""Test match function for Kronik DataFrames"""
	kro_df = simple_df("tests/kro_fixture.kro")
	redund_df = kro_df.copy().iloc[0:1000, ]
	expected_matches = 1042
	actual_matches = sum(redund_df.apply(match_rt_mass, axis=1, other_df=kro_df, rt_diff=1))
	assert actual_matches == expected_matches, "Kronik match function did not work properly."
	