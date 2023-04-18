from msions.msplot import plot_data
from msions.mzml import tic_df

def test_plot_data():
    """Test plotting of TIC"""
    ms1_df = tic_df("tests/mzml_fixture.mzml")
    string_output = plot_data("tests/mzml_fixture.mzml", color="black", 
                              fig_params=[12,8,500])
    fig = plot_data(ms1_df, color="black", fig_params=[12,8,500])
    assert fig is None, "Plot did not generate properly"
    assert string_output == fig, "File input was not processed correctly."

def test_return_dfs():
    df, sumid_feat_df, id_df, feat_df = plot_data("tests/mzml_fixture.mzml", 
              "tests/hk_fixture.hk", 
              "tests/large_fixtures/elib_fixture.elib", 
              method="DIA", data_type="both", 
              return_dfs=True)
    expected_df_rows = 2
    actual_df_rows = df.shape[0]
    assert actual_df_rows == expected_df_rows, "MS1 DataFrame was not created correctly."
    expected_sumid_val = 40.3431
    actual_sumid_val = sumid_feat_df.IT[1].round(4)
    assert actual_sumid_val == expected_sumid_val, "Identified Hardklor features were not summarized correctly."
    expected_id_rows = 32686
    actual_id_rows = id_df.shape[0]
    assert actual_id_rows == expected_id_rows, "EncyclopeDIA DataFrame was not created correctly."
    expected_feat_rows = 383
    actual_feat_rows = feat_df.shape[0]
    assert actual_feat_rows == expected_feat_rows, "Hardklor DataFrame was not created correctly."

