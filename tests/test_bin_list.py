from msions.utils import bin_list


def test_bin_list():
	"""Test function that creates bin edges"""
	mz_start = 399
	mz_end = 1005
	mz_bin_size = 4
	mz_bin_mult = 1.0005
	expected_num_edges = 152
	actual_num_edges = len(bin_list(mz_start, mz_end, mz_bin_size, mz_bin_mult))
	assert actual_num_edges == expected_num_edges, "List of bin edges was not created properly."