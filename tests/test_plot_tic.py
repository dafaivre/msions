from msions.msplot import plot_tic
from msions.mzml import tic_df

def test_plot_tic():
    """Test plotting of TIC"""
    ms1_df = tic_df("tests/mzml_fixture.mzml")
    fig = plot_tic(ms1_df)
    assert fig is None, "Plot did not generate properly"
