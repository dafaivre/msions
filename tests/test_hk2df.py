from msions.hardklor import hk2df
import numpy as np
import pandas as pd

def test_hk2df():
	"""Test DataFrame creation from a Hardklor file"""
	expected_type = "DataFrame"
	actual_type = type(hk2df("tests/hk_fixture.hk")).__name__
	assert actual_type == expected_type, "DataFrame was not created correctly. Check format of file."

def test_byint():
	"""Test DataFrame creation and sorting from a Hardklor file"""
	expected_type = "DataFrame"
	expected_int = 589898
	actual = hk2df("tests/hk_fixture.hk", by_int=True)
	actual_type = type(actual).__name__
	actual_int = actual.loc[0,"intensity"]
	assert actual_type == expected_type, "DataFrame was not created correctly. Check format of file."
	assert actual_int == expected_int, "Test Failed. DataFrame may not be sorted correctly."
