from msions.percolator import parse_peps

def test_parse_peps():
	"""Test list of dictionaries containing parsed PSMs"""
	expected_len = 10
	actual_len = len(parse_peps("tests/pep_fixture.pout.xml"))
	assert actual_len == expected_len, "Peptides were not parsed correctly."
