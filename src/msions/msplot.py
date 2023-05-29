"""
This module contains functions that are useful for plotting MS data in Python.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # for despine of plots
from msions.mzml import tic_df
from msions.encyclopedia import dia_df
from msions.hardklor import hk2df
from msions.encyclopedia import match_hk
from msions.hardklor import summarize_df
import numpy as np
from typing import List, Union


def plot_data(mzml_input: Union[pd.DataFrame, str], 
			  feat_input: Union[pd.DataFrame, str] = None, 
			  id_input: Union[pd.DataFrame, str] = None,
			  method: str = None,
			  data_type = "TIC",
			  stats = None,
			  return_dfs = False,
			  color: Union[str, List[str]] = ["black", "#1f77b4"], 
			  no_labels: bool = False, alpha: float = 1.0, 
			  fig_params: List[float] = None):
	"""
	Plots TIC against retention time.

	Parameters
	----------
	mzml_input : pd.DataFrame or str
		The pandas DataFrame containing retention time and TIC or the mzML file.
	feat_input : pd.DataFrame or str
		The pandas DataFrame containing retention time and TIC or the Hardklor/Kronik file.
	id_input : pd.DataFrame or str
		The pandas DataFrame containing retention time and TIC or the EncyclopeDIA/Percolator file.
	method : str
		Type of acquisition method (e.g., "DIA")
	data_type : str
		Data chosen for plot ("TIC", "ions", "both")
	stats : bool
		Print stats for data
	return_dfs : bool
		Return calculated DataFrames if True
	color: List[str]
		The lists of colors for the line plots.
	no_labels: bool
		Removes ticks and labels.
	alpha: float
		Changes the alpha value for the line plot.
	fig_params: List[float]
		Sets the figure size and optionally the dpi.
	
	Examples
	-------
	>>> from msions.msplot import plot_tic
	>>> from msions.mzml import tic_df
	>>> import matplotlib.pyplot as plot
	>>> ms1_df = tic_df("test.mzML")
	>>> plot_data(ms1_df)
	>>> plt.show()
	""" 
	# create blank variables for returning
	df = ""
	feat_df = ""
	id_df = ""
	sumid_feat_df = ""
	title_txt = []

	# if it's an mzML file
	if isinstance(mzml_input, str):
		# create mzML data frame
		df = tic_df(mzml_input)
		
	# if it's a data frame already
	else:
		df = mzml_input	

	if isinstance(color, str):
		color = [color]

	# change figure size
	if fig_params is not None:
		if len(fig_params) == 2:
			plt.figure(figsize=(fig_params[0], fig_params[1]))
		elif len(fig_params) == 3:
			plt.figure(figsize=(fig_params[0], fig_params[1]), dpi=fig_params[2])

	if data_type.lower() == "ions":
		if stats.lower() == "print":
			# find total ions across all 
			print("Total # of ions: %.2e" % sum(df.ions))
			
		elif stats.lower() == "title":
			# find total ions across all 
			title_txt.append("Total # of ions: %.2e\n" % sum(df.ions)) 
		
		if fig_params is None:
			# define figure size
			plt.figure(figsize=(10,8))

		# plot ions
		plt.plot(df['rt'], df['ions'], color=color[0], alpha=alpha)	
	
	elif data_type.lower() == "both":
		if stats == "print":
			# find totals across all scans
			print("Total Ion Current (TIC): %.2e \t Total # of ions: %.2e" % (sum(df.TIC), sum(df.ions)))

		elif stats.lower() == "title":
			# find total ions across all 
			title_txt.append("Total Ion Current (TIC): %.2e\n" % sum(df.TIC)) 
			title_txt.append("Total # of ions: %.2e\n" % sum(df.ions))

		if fig_params is None:
			plt.figure(figsize=(16, 6))

		# plot TIC
		plt.subplot(1, 2, 1)
		plt.plot(df['rt'], df['TIC'], color=color[0], alpha=alpha)
	
		# plot ions
		plt.subplot(1, 2, 2)
		plt.plot(df['rt'], df['ions'], color=color[0], alpha=alpha)	

	else:
		if stats == "print":
			# find total ion current across all 
			print("Total Ion Current (TIC): %.2e" % sum(df.TIC))

		elif stats.lower() == "title":
			# find total ions across all 
			title_txt.append("Total Ion Current (TIC): %.2e\n" % sum(df.TIC))

		if fig_params is None:
			# define figure size
			plt.figure(figsize=(10,8))

		# plot TIC
		plt.plot(df['rt'], df['TIC'], color=color[0], alpha=alpha)

	# if feature file and ID file is given
	if feat_input is not None and id_input is not None:
		if isinstance(feat_input, str):
			# check if Hardklor file
			if feat_input[-2:] == "hk":
				feat_df = hk2df(feat_input)
		else:
			feat_df = feat_input

		if isinstance(id_input, str):
			# check if EncyclopeDIA file
			if id_input[-4:] == "elib":
				id_df = dia_df(id_input)
		else:
			id_df = id_input

		if method == "DIA":
			# if scans have been summed already
			if len(np.unique(feat_df.scan_num)) == len(feat_df.scan_num):
				sumid_feat_df = feat_df

			# if features have been matched already
			elif "in_encyclo" in feat_df.columns:
				# create DataFrame of only identified features
				id_feat_df = feat_df[feat_df["in_encyclo"] > 0].reset_index(drop=True)			

				# summarize identified features DataFrame
				sumid_feat_df = summarize_df(id_feat_df, full_ms1_df=df)
			else:
				# find Hardklor/encyclopeDIA match
				feat_df["in_encyclo"] = feat_df.apply(match_hk, axis=1, other_df=id_df)

				# create DataFrame of only identified features
				id_feat_df = feat_df[feat_df["in_encyclo"] > 0].reset_index(drop=True)

				# summarize identified features DataFrame
				sumid_feat_df = summarize_df(id_feat_df, full_ms1_df=df)

			if data_type.lower() == "ions":
				if stats.lower() == "print":
					# find identified ions
					print("Ions mapped to peptides: %.2e" % sum(sumid_feat_df.ions))

					# calculate ratio of identified ions to total ions
					print("%.1f%% of the signal" % float(sum(sumid_feat_df.ions)/sum(df.ions)*100))

					# print number of peptide IDs
					print ("Number of peptide IDs: %.0f" % len(id_df))

				elif stats.lower() == "title":
					title_txt[0] += "Ions mapped to peptides: %.2e\n" % sum(sumid_feat_df.ions)
					title_txt[0] += "%.1f%% of the signal\n" % float(sum(sumid_feat_df.ions)/sum(df.ions)*100)
					title_txt[0] += "Number of peptide IDs: %.0f\n" % len(id_df)

				# plot ID'd ions
				plt.plot(sumid_feat_df['rt'], sumid_feat_df['ions'], color=color[1], alpha=alpha)

			elif data_type.lower() == "both":
				if stats.lower() == "print":
					# find identified signals
					print("ID'd TIC: %.2e \t\t\t Ions mapped to peptides: %.2e" % (sum(sumid_feat_df.TIC), sum(sumid_feat_df.ions)))

					# calculate ratio of identified signals to total signal
					print("%.1f%% of the signal \t\t\t %.1f%% of the signal" % (float(sum(sumid_feat_df.TIC)/sum(df.TIC)*100),
						  float(sum(sumid_feat_df.ions)/sum(df.ions)*100)))

					# print number of peptide IDs
					print ("Number of peptide IDs: %.0f" % len(id_df))

				elif stats == "title":
					title_txt[0] += "ID'd TIC: %.2e\n" % sum(sumid_feat_df.TIC)
					title_txt[0] += "%.1f%% of the signal\n" % float(sum(sumid_feat_df.TIC)/sum(df.TIC)*100)
					title_txt[0] += "Number of peptide IDs: %.0f\n" % len(id_df)
					title_txt[1] += "Ions mapped to peptides: %.2e\n" % sum(sumid_feat_df.ions)
					title_txt[1] += "%.1f%% of the signal\n\n" % float(sum(sumid_feat_df.ions)/sum(df.ions)*100)
					# title_txt[1] += "Number of peptide IDs: %.0f" % len(id_df)				
			

				# plot TIC
				plt.subplot(1, 2, 1)
				plt.plot(sumid_feat_df['rt'], sumid_feat_df['TIC'], color=color[1], alpha=alpha)

				# plot ions
				plt.subplot(1, 2, 2)
				plt.plot(sumid_feat_df['rt'], sumid_feat_df['ions'], color=color[1], alpha=alpha)	

			else:
				if stats == "print":
					# find identified total ion current
					print("ID'd TIC: %.2e" % sum(sumid_feat_df.TIC))

					# calculate ratio of identified ion current to total ion current
					print("%.1f%% of the signal" % float(sum(sumid_feat_df.TIC)/sum(df.TIC)*100))

					# print number of peptide IDs
					print ("Number of peptide IDs: %.0f" % len(id_df))

				elif stats == "title":
					title_txt[0] += "ID'd TIC: %.2e\n" % sum(sumid_feat_df.TIC)
					title_txt[0] += "%.1f%% of the signal\n" % float(sum(sumid_feat_df.TIC)/sum(df.TIC)*100)
					title_txt[0] += "Number of peptide IDs: %.0f\n" % len(id_df)


				# plot ID'd TIC
				plt.plot(sumid_feat_df['rt'], sumid_feat_df['TIC'], color=color[1], alpha=alpha)

	if data_type.lower() == "ions":
		plt.xticks(fontsize=14)
		plt.xlabel("Time (min)", fontsize=18)
		plt.yticks(fontsize=14)
		plt.ylabel("Ions", fontsize=18)
		if stats == "title":
			plt.title(title_txt[0], loc="left", fontsize=18)

	elif data_type.lower() == "both":

		plt.subplot(1, 2, 1)
		plt.xticks(fontsize=14)
		plt.xlabel("Time (min)", fontsize=18)
		plt.yticks(fontsize=14)
		plt.ylabel("Total Ion Current", fontsize=18)
		if stats == "title":
			plt.title(title_txt[0], loc="left", fontsize=18)

		# gives scientific notation
		plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

		# removes right side & top of plot
		sns.despine()
		
		# x- and y-axis start at 0
		plt.xlim(left=0)
		plt.ylim(bottom=0)

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

		plt.subplot(1, 2, 2)
		plt.xticks(fontsize=14)
		plt.xlabel("Time (min)", fontsize=18)
		plt.yticks(fontsize=14)
		plt.ylabel("Ions", fontsize=18)

		if stats == "title":
			plt.title(title_txt[1], loc="left", fontsize=18)

	else:
		plt.xticks(fontsize=14)
		plt.xlabel("Time (min)", fontsize=18)
		plt.yticks(fontsize=14)
		plt.ylabel("Total Ion Current", fontsize=18)
		if stats == "title":
			plt.title(title_txt[0], loc="left", fontsize=18)

	# gives scientific notation
	plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

	# removes right side & top of plot
	sns.despine()
	
	# x- and y-axis start at 0
	plt.xlim(left=0)
	plt.ylim(bottom=0)

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

	if return_dfs:
		return df, sumid_feat_df, id_df, feat_df
