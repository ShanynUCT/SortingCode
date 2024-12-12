import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import sys


plt.rcParams['axes.linewidth'] = 1
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['DejaVu Serif']
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
# axis labels bold
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams["legend.fancybox"] = True
plt.rcParams["legend.edgecolor"] = 'white'
plt.rcParams["legend.handlelength"] = 2
plt.rcParams["legend.handleheight"] = 0
plt.rcParams["legend.handletextpad"] = 0
plt.rcParams["legend.markerscale"] = 1
#incude tick marks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 2.5
plt.rcParams['ytick.minor.size'] = 2.5
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.minor.width'] = 1
#show minor ticks
plt.rcParams['xtick.minor.visible'] = True

# read in three text files from "/Users/shanyn/Documents/PhD/Codes/SortingCode/steppedSortBeam/" with names insync_R27.txt, outsync_R27.txt, and prompt_R27.txt
#first column is energy and then counts, space separated
#plot the energy vs counts for each file on the same plot

data_dir = "/Users/shanyn/Documents/PhD/Codes/SortingCode/steppedSortBeam/"
insync_data = 'insync_R27.txt'
outsync_data = 'outsync_R27.txt'
prompt_data = 'prompt_R27.txt'

insync_data_df = pd.read_csv(data_dir + insync_data, sep=' ', header=None)
insync_data_df.columns = ['energy', 'counts']
insync_data_df = insync_data_df.reset_index(drop=True)

outsync_data_df = pd.read_csv(data_dir + outsync_data, sep=' ', header=None)
outsync_data_df.columns = ['energy', 'counts']
outsync_data_df = outsync_data_df.reset_index(drop=True)

prompt_data_df = pd.read_csv(data_dir + prompt_data, sep=' ', header=None)
prompt_data_df.columns = ['energy', 'counts']
prompt_data_df = prompt_data_df.reset_index(drop=True)

# plot the energy vs counts in 4000 bins between 0 and 8000 to increase binning
fig, ax = plt.subplots()
ax.hist(insync_data_df['energy'], bins=2000, weights=insync_data_df['counts'], histtype='step', label='Insync', color='b')
ax.hist(outsync_data_df['energy'], bins=2000, weights=outsync_data_df['counts'], histtype='step', label='Outsync', color='r')
ax.hist(prompt_data_df['energy'], bins=2000, weights=prompt_data_df['counts'], histtype='step', label='Prompt', color='k')
ax.set_xlabel('Energy (keV)', fontsize=16, fontweight='bold')
ax.set_ylabel('Counts', fontsize=16, fontweight='bold')
plt.legend(loc='upper right', fontsize=16)
plt.ylim(0, 500)
plt.show()
