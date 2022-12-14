from msions.percolator import id_scans
import msions.mzml as mzml

def test_tic_df():
	"""Test DataFrame creation from an mzML file"""
	expected_ids = 352
	ms2_tic_df = mzml.tic_df("tests/large_fixtures/short_DDA_file.mzML", level="2")
	id_scans("tests/large_fixtures/DDA_percolator.target.peptides.txt", ms2_tic_df)
	actual_ids = sum(ms2_tic_df["IDd"])
	assert actual_ids == expected_ids, "Scans were not identified properly."
	