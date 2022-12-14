from msions.encyclopedia import match_hk
from msions.encyclopedia import dia_df 
import msions.hardklor as hk 


def test_match_hk():
	"""Test match function for Hardklor and EncyclopeDIA output"""
	hk_df = hk.hk2df("tests/hk_fixture.hk") 
	encyclo_df = dia_df("tests/large_fixtures/elib_fixture.elib") 
	expected_matches = 140
	actual_matches = sum(hk_df.apply(match_hk, axis=1, other_df=encyclo_df))
	assert actual_matches == expected_matches, "Hardklor and EncyclopeDIA match function did not work properly."
	