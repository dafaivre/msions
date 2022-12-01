from msions.hardklor import summarize_df
from msions.hardklor import hk2df
import numpy as np
import pandas as pd

def test_summarize_df():
	"""Test summarized DataFrame"""
	expected_rows = 2
	actual_rows = summarize_df(hk2df("tests/hk_fixture.hk")).shape[0]
	assert actual_rows == expected_rows, "DataFrame was not summarized correctly."
    
