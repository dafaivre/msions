from msions.utils import bin_list
from msions.utils import bin_data
import pandas as pd
import pymzml

def test_bin_data():
	"""Test function that creates a binned DataFrame"""
	# create run object
	run = pymzml.run.Reader("tests/mzml_fixture.mzML")

	# initiate peak DataFrame
	peak_df = pd.DataFrame(columns=["mz", "ips", "rt"])

	# loop through spectra
	for spectra in run:
		if spectra.ms_level == 1:
			peak_array = pd.DataFrame(spectra.peaks("centroided")).rename(columns={0: "mz", 1: "ips"})
			peak_array["mz"] = peak_array["mz"].round(4)
			peak_array["rt"] = spectra.scan_time[0]
			peak_df = pd.concat([peak_df, peak_array])
	
	# initiate bin variables
	mz_start = 399
	mz_end = 1005
	mz_bin_size = 4
	mz_bin_mult = 1.0005
	bin_mz_list = bin_list(mz_start, mz_end, mz_bin_size, mz_bin_mult)

	# test binning function
	bin_df = bin_data(peak_df, type="mz", bin_mz_list=bin_mz_list)
	expected_num_rows = 302
	actual_num_rows = bin_df.shape[0]
	assert actual_num_rows == expected_num_rows, "Binned DataFrame was not created properly."
