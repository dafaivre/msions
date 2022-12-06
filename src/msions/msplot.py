"""
This module contains functions that are useful for plotting MS data in Python.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # for despine of plots
from typing import List


def plot_tic(df: pd.DataFrame, color: str = "#1f77b4", xlabel: str = "Time (min)", 
			 ylabel: str = "Total Ion Current", no_labels: bool = False, alpha: float = 1.0, 
			 ymax: float = None, fig_params: List[float] = None):
	"""
	Plots TIC against retention time from a pandas DataFrame.

	Parameters
	----------
	df : pd.DataFrame
		The pandas DataFrame containing retention time and TIC.
	color: str
		The color for the line plot.
	xlabel: str
		Title for x-axis.
	ylabel: str
		Title for y-axis.
	no_labels: bool
		Removes ticks and labels.
	alpha: float
		Changes the alpha value for the line plot.
	ymax: float
		Sets the maximum value of the y-axis.
	fig_params: List[float]
		Sets the figure size and optionally the dpi.
	
	Examples
	-------
	>>> from msions.msplot import plot_tic
	>>> from msions.mzml import tic_df
	>>> import matplotlib.pyplot as plot
	>>> ms1_df = tic_df("test.mzML")
	>>> plot_tic(ms1_df)
	>>> plt.show()
	""" 
	# change figure size
	if fig_params is not None:
		if len(fig_params) == 2:
			plt.figure(figsize=(fig_params[0], fig_params[1]))
		elif len(fig_params) == 3:
			plt.figure(figsize=(fig_params[0], fig_params[1]), dpi=fig_params[2])

	# plot TIC
	plt.plot(df['rt'], df['TIC'], color=color, alpha=alpha)

	# gives scientific notation
	plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
	
	# x- and y-axis start at 0
	plt.xlim(left=0)
	plt.ylim(bottom=0)

	# defines max of y-axis
	if ymax is not None:
		plt.ylim(0, ymax)

	# removes right side & top of plot
	sns.despine()

	# set x- and y-axis titles
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)

	# removes labels		
	if no_labels:
		# get rid of x- and y-axis titles
		plt.xlabel(None)
		plt.ylabel(None)
		
		# gives numbers instead of scientific notation, needed if trying to get rid of all labels
		plt.ticklabel_format(style="plain")

		# draw ticks and labels
		plt.tick_params(
			axis='both',
			which='both',  # both major and minor ticks are affected
			bottom=False,
			top=False,   
			left=False,
			right=False,
			labelbottom=False,
			labelleft=False,
			labeltop=False,
			labelright=False) 


def plot_ions(df: pd.DataFrame, color: str = "#1f77b4", xlabel: str = "Time (min)", 
			  ylabel: str = "Ions", no_labels: bool = False, alpha: float = 1.0, 
			  ymax: float = None, fig_params: List[float] = None):
	"""
	Plots ions against retention time from a pandas DataFrame.
	
	Parameters
	----------
	df : pd.DataFrame
		The pandas DataFrame containing retention time and ion counts.
	color: str
		The color for the line plot.
	xlabel: str
		Title for x-axis.
	ylabel: str
		Title for y-axis.
	no_labels: bool
		Removes ticks and labels.
	alpha: float
		Changes the alpha value for the line plot.
	ymax: float
		Sets the maximum value of the y-axis.
	fig_params: List[float]
		Sets the figure size and optionally the dpi.
	
	Examples
	-------
	>>> from msions.msplot import plot_ions
	>>> from msions.mzml import tic_df
	>>> import matplotlib.pyplot as plot
	>>> ms1_df = tic_df("test.mzML")
	>>> plot_ions(ms1_df)
	>>> plt.show()
	""" 
	# change figure size
	if fig_params is not None:
		if len(fig_params) == 2:
			plt.figure(figsize=(fig_params[0], fig_params[1]))
		elif len(fig_params) == 3:
			plt.figure(figsize=(fig_params[0], fig_params[1]), dpi=fig_params[2])

	# plot TIC
	plt.plot(df['rt'], df['ions'], color=color, alpha=alpha)

	# gives scientific notation
	plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

	# x- and y-axis start at 0
	plt.xlim(left=0)
	plt.ylim(bottom=0)

	# defines max of y-axis
	if ymax is not None:
		plt.ylim(0, ymax)

	# removes right side & top of plot
	sns.despine()

	# set x- and y-axis titles
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)

	# removes labels		
	if no_labels:
		# get rid of x- and y-axis titles
		plt.xlabel(None)
		plt.ylabel(None)

		# gives numbers instead of scientific notation, needed if trying to get rid of all labels
		plt.ticklabel_format(style="plain")

		# draw ticks and labels
		plt.tick_params(
			axis='both',
			which='both',  # both major and minor ticks are affected
			bottom=False,
			top=False,   
			left=False,
			right=False,
			labelbottom=False,
			labelleft=False,
			labeltop=False,
			labelright=False)
