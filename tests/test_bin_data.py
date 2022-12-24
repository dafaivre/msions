from msions.utils import bin_list
# from msions.utils import bin_data
import msions.hardklor as hk

def test_bin_data():
	"""Test function that creates a binned DataFrame"""
	# hk_df = hk.hk2df("tests/hk_fixture.hk") 
	mz_start = 399
	mz_end = 1005
	mz_bin_size = 4
	mz_bin_mult = 1.0005
	bin_mz_list = bin_list(mz_start, mz_end, mz_bin_size, mz_bin_mult)
	# bin_data(peak_df, type="mz", bin_mz_list=bin_mz_list)
	expected_num_edges = 152
	actual_num_edges = len(bin_mz_list)
	assert actual_num_edges == expected_num_edges, "List of bin edges was not created properly."

