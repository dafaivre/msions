from msions.msplot import plot_ions
from msions.mzml import tic_df

def test_plot_ions():
    """Test plotting of ions"""
    ms1_df = tic_df("tests/mzml_fixture.mzml")
    fig = plot_ions(ms1_df)
    assert fig is None, "Plot did not generate properly"
