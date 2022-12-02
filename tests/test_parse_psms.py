from msions.msxml import parse_psms

def test_parse_psms():
	"""Test list of dictionaries containing parsed PSMs"""
	expected_len = 10
	actual_len = len(parse_psms("tests/psm_fixture.pout.xml"))
	assert actual_len == expected_len, "PSMs were not parsed correctly."
