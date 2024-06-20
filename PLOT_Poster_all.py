from cmath import nan
import collections
from ctypes import sizeof
from operator import index

from tokenize import Double
from unittest import skip
from venv import create
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot as up
import csv as csv
import io as io
import sys, os
from math import sin, cos, pi, log, floor
import cProfile
import re
import pylab
import pyroot as pr
import decimal
import time as time
%matplotlib inline
#%jsroot on
from sklearn.preprocessing import MinMaxScaler
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

#--------------------------------------------------------------------------------
#file_dir = sys.argv[1]
file_dir = /home/shan/Documents/PhD/ExperimentalResults/220615/run43-66MeV_proton_beam-combo_run-D0_315deg_22cm-Carbon_Target_1cm_depth-20min/TimeSyncResults/POSTER/
save_dir = file_dir + '/' + 'POSTER/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))
#--------------------------------------------------------------------------------


# ---------------------------- 1. Read in TimeDiff data from the LaBrPOLARIS.py script before concatenation -----------------------------------
timediff_dir = os.path.join( file_dir, '../' )
filename_timediff_LaBr0 = timediff_dir + '/' + 'LaBrData_Det0_df_BeforeCONCAT_forTD.csv'
filename_timediff_LaBr1 = timediff_dir + '/' + 'LaBrData_Det1_df_BeforeCONCAT_forTD.csv'
filename_timediff_LaBr2 = timediff_dir + '/' + 'LaBrData_Det2_df_BeforeCONCAT_forTD.csv'
filename_timediff_LaBr3 = timediff_dir + '/' + 'LaBrData_Det3_df_BeforeCONCAT_forTD.csv'

filename_timediff_POLARIS = timediff_dir + '/' + 'mod_event_data_df_BeforeCONCAT_forTD.csv'

TD_LaBr0 = pd.read_csv(filename_timediff_LaBr0, sep = ',', header = 0)
TD_LaBr0.columns = ['Timestamp_LaBr',  'EnergySlow_LaBr_1stOrderLowRange',  'Detector', 'x',   'y',   'z']

TD_LaBr1 = pd.read_csv(filename_timediff_LaBr1, sep = ',', header = 0)
TD_LaBr1.columns = ['Timestamp_LaBr',  'EnergySlow_LaBr_1stOrderLowRange',  'Detector','x',   'y',   'z']

TD_LaBr2 = pd.read_csv(filename_timediff_LaBr2, sep = ',', header = 0)
TD_LaBr2.columns = ['Timestamp_LaBr',  'EnergySlow_LaBr_1stOrderLowRange', 'Detector',   'x',   'y',   'z']

TD_LaBr3 = pd.read_csv(filename_timediff_LaBr3, sep = ',', header = 0)
TD_LaBr3.columns = ['Timestamp_LaBr',  'EnergySlow_LaBr_1stOrderLowRange', 'Detector',   'x',   'y',   'z']

TD_POLARIS = pd.read_csv(filename_timediff_POLARIS, sep = ',', header = 0)
TD_POLARIS.columns = ['Timestamp_LaBr',  'EnergySlow_LaBr_1stOrderLowRange', 'Detector',   'x',   'y',   'z']

# find the frequency of the energy counts per unit timestamp
CpuT_LaBr0 = TD_LaBr0[(TD_LaBr0['EnergySlow_LaBr_1stOrderLowRange'] > 0) & (TD_LaBr0['EnergySlow_LaBr_1stOrderLowRange'] < 10000)]
CpuT_LaBr0 = CpuT_LaBr0.groupby(['EnergySlow_LaBr_1stOrderLowRange', 'Timestamp_LaBr']).size().reset_index(name = 'counts')
TD_LaBr0['Frequency'] = TD_LaBr0.groupby('Timestamp_LaBr')['Timestamp_LaBr'].transform('count')

CpuT_LaBr1 = TD_LaBr1[(TD_LaBr1['EnergySlow_LaBr_1stOrderLowRange'] > 0) & (TD_LaBr1['EnergySlow_LaBr_1stOrderLowRange'] < 10000)]
CpuT_LaBr1 = CpuT_LaBr1.groupby(['EnergySlow_LaBr_1stOrderLowRange', 'Timestamp_LaBr']).size().reset_index(name = 'counts')
TD_LaBr1['Frequency'] = TD_LaBr1.groupby('Timestamp_LaBr')['Timestamp_LaBr'].transform('count')

CpuT_LaBr2 = TD_LaBr2[(TD_LaBr2['EnergySlow_LaBr_1stOrderLowRange'] > 0) & (TD_LaBr2['EnergySlow_LaBr_1stOrderLowRange'] < 10000)]
CpuT_LaBr2 = CpuT_LaBr2.groupby(['EnergySlow_LaBr_1stOrderLowRange', 'Timestamp_LaBr']).size().reset_index(name = 'counts')
TD_LaBr2['Frequency'] = TD_LaBr2.groupby('Timestamp_LaBr')['Timestamp_LaBr'].transform('count')

CpuT_LaBr3 = TD_LaBr3[(TD_LaBr3['EnergySlow_LaBr_1stOrderLowRange'] > 0) & (TD_LaBr3['EnergySlow_LaBr_1stOrderLowRange'] < 10000)]
CpuT_LaBr3 = CpuT_LaBr3.groupby(['EnergySlow_LaBr_1stOrderLowRange', 'Timestamp_LaBr']).size().reset_index(name = 'counts')
TD_LaBr3['Frequency'] = TD_LaBr3.groupby('Timestamp_LaBr')['Timestamp_LaBr'].transform('count')

TD_POLARIS['Frequency'] = TD_POLARIS.groupby('Timestamp_LaBr')['Timestamp_LaBr'].transform('count')

# plot the frequency of the energy counts per unit timestamp 
# from the TD_LaBr df find the counts per unit time (1s) and plot it vs EnergySlow_LaBr_1stOrderLowRange


timediffLaBr0entries_plot = ((TD_LaBr0['Timestamp_LaBr'])*10e-8*1000).diff()
timediffLaBr1entries_plot = ((TD_LaBr1['Timestamp_LaBr'])*10e-8*1000).diff()
timediffLaBr2entries_plot = ((TD_LaBr2['Timestamp_LaBr'])*10e-8*1000).diff()
timediffLaBr3entries_plot = ((TD_LaBr3['Timestamp_LaBr'])*10e-8*1000).diff()
timediffPOLARISentries_plot = ((TD_POLARIS['Timestamp_LaBr'])*10e-8*1000).diff()
mintime = timediffPOLARISentries_plot.min()
meantimePOLARIS = timediffPOLARISentries_plot.mean()
meantimePOLARIS = round(meantimePOLARIS, 2)
meantimeLaBr0 = timediffLaBr0entries_plot.mean()
meantimeLaBr0 = round(meantimeLaBr0, 2)
meantimeLaBr1 = timediffLaBr1entries_plot.mean()
meantimeLaBr1 = round(meantimeLaBr1, 2)
meantimeLaBr2 = timediffLaBr2entries_plot.mean()
meantimeLaBr2 = round(meantimeLaBr2, 2)
meantimeLaBr3 = timediffLaBr3entries_plot.mean()
meantimeLaBr3 = round(meantimeLaBr3, 2)

# plot the previous two figures on the same axes to compare the time differences between consecutive entries, displaying the mean value for each
plt.figure(figsize = (30,15))
plt.hist(timediffPOLARISentries_plot, bins = 1000, histtype = 'step', label = 'POLARIS', range=(mintime, 2), color = 'royalblue', linewidth = 2)
plt.hist(timediffLaBr0entries_plot, bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L0', range=(mintime, 2), color = 'darkorange', linewidth = 2)
plt.hist(timediffLaBr1entries_plot, bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L1', range=(mintime, 2), color = 'darkgreen', linewidth = 2)
plt.hist(timediffLaBr2entries_plot, bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L2', range=(mintime, 2), color = 'magenta', linewidth = 2)
plt.hist(timediffLaBr3entries_plot, bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L3', range=(mintime, 2), color = 'violet', linewidth = 2)

plt.title('Time Difference Between Consecutive Entries For Each Detector', fontsize = 25)
plt.xlabel('$\\Delta$t ($\\mu$s)', fontsize = 30)
plt.ylim(0,225000 )
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.text(0.5, 0.5, 'Avg. $\\Delta$t POLARIS = ' + str(meantimePOLARIS) + ' $\\mu$s', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize = 30)
plt.text(0.5, 0.4, 'Avg. $\\Delta$t LaBr$_3$(Ce) L0 = ' + str(meantimeLaBr0)+ ' $\\mu$s', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize = 30)
plt.text(0.5, 0.3, 'Avg. $\\Delta$t LaBr$_3$(Ce) L1 = ' + str(meantimeLaBr1)+ ' $\\mu$s', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize = 30)
plt.text(0.5, 0.2, 'Avg. $\\Delta$t LaBr$_3$(Ce) L2 = ' + str(meantimeLaBr2)+ ' $\\mu$s', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize = 30)
plt.text(0.5, 0.1, 'Avg. $\\Delta$t LaBr$_3$(Ce) L3 = ' + str(meantimeLaBr3)+ ' $\\mu$s', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize = 30)
plt.ylabel('Counts', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.savefig(save_dir+ 'Time_Diff_Between_Consecutive_Entries_For_Each_Detector_AfterSync.png')

# plot the Timestamp_LaBr for TD_LaBr0, TD_LaBr1, TD_LaBr2, TD_LaBr3 and TD_POLARIS on the same axes 
mintime = 1e12
maxtime = TD_POLARIS['Timestamp_LaBr'].max()

plt.figure(figsize = (30,15))
plt.hist(TD_LaBr0['Timestamp_LaBr'], bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L0', range=(mintime, maxtime), color = 'darkorange', linewidth = 2)
plt.hist(TD_LaBr1['Timestamp_LaBr'], bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L1', range=(mintime, maxtime), color = 'darkgreen', linewidth = 2)
plt.hist(TD_LaBr2['Timestamp_LaBr'], bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L2', range=(mintime, maxtime), color = 'magenta', linewidth = 2)
plt.hist(TD_LaBr3['Timestamp_LaBr'], bins = 1000, histtype = 'step', label = 'LaBr$_3$(Ce) L3', range=(mintime, maxtime), color = 'violet', linewidth = 2)
plt.hist(TD_POLARIS['Timestamp_LaBr'], bins = 1000, histtype = 'step', label = 'POLARIS', range=(mintime, maxtime), color = 'royalblue', linewidth = 2)
plt.title('Timestamps For Each Detector', fontsize = 25)
plt.xlabel('Timestamp ($\\mu$s)', fontsize = 30)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.ylabel('Counts', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.savefig(save_dir+ 'Timestamps_For_Each_Detector_AfterSync.png')

#'' --------------------------------------------------------------------------------

'''file_nameCLF = 'Coincidence_df_CLF.txt'
coincidence = pd.read_csv(file_dir + '/' + file_nameCLF, sep = '\t', header = 0)
coincidence.columns = ['Index', 'Time', 'POLARIS_energy', 'x', 'y', 'z', 'LaBr_Energy', 'Window', 'TimeDiff', 'TimeDiff_CumSum', 'theta1_deg']

# DATA BEFORE CLF BUT AFTER WINDOW CUT
file_name_noCLF = 'BeforeCLF_df.txt'
coincidence_no_CLF = pd.read_csv(file_dir + '/' + file_name_noCLF, sep = '\t', header = 0)
coincidence_no_CLF.columns = ['E0', 'E1', 'theta1_deg']
coincidencePOLARISNOCLFplotE1 = coincidence_no_CLF[(coincidence_no_CLF['E1'] > 0) & (coincidence_no_CLF['E1'] < 10000)]
coincidencePOLARISNOCLFplotE1 = coincidencePOLARISNOCLFplotE1.groupby(['E1', 'theta1_deg']).size().reset_index(name = 'counts')'''
#scaler = MinMaxScaler()
#coincidencePOLARISNOCLFplotE1['counts'] = scaler.fit_transform(coincidencePOLARISNOCLFplotE1[['counts']])


'''plt.figure(figsize = (30, 15))
plt.scatter(coincidencePOLARISNOCLFplotE1['theta1_deg'], coincidencePOLARISNOCLFplotE1['E1'] , c = coincidencePOLARISNOCLFplotE1['counts'] , cmap = 'turbo',vmin=1, s=0.5)
plt.xlim(0, 180)
plt.ylim(0, 6000)
plt.ylabel('Energy (keV)', fontsize = 30)
plt.xlabel('$\\theta$ (deg)', fontsize = 30)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.title('E$_{1}$ vs $\\theta$', fontsize = 30)
plt.grid(axis='both', alpha=0.5, linestyle='--', linewidth=0.5)
cbar = plt.colorbar()
plt.savefig(save_dir + 'E1_vs_theta_NOCLF' + '.png')
print(' PLOT SAVED: ' + save_dir + 'E1_vs_theta_NOCLF' + '.png')'''


'''coincidencePOLARISplotE = coincidence[(coincidence['POLARIS_energy'] > 0) & (coincidence['POLARIS_energy'] < 10000)]
coincidencePOLARISplotE = coincidencePOLARISplotE.groupby(['POLARIS_energy', 'theta1_deg']).size().reset_index(name = 'counts')'''
# take sklearn min max scalar to scale the counts to a range of 0 to 1
#coincidencePOLARISplotE['counts'] = scaler.fit_transform(coincidencePOLARISplotE[['counts']])



'''plt.figure(figsize = (30, 15))
plt.scatter(coincidencePOLARISplotE['theta1_deg'], coincidencePOLARISplotE['POLARIS_energy'] , c = coincidencePOLARISplotE['counts'] , cmap = 'turbo',vmin=1, s=0.5)
plt.xlim(0, 180)
plt.ylim(0, 6000)
plt.xlabel('$\\theta$ (deg)', fontsize = 30)
plt.ylabel('Energy (keV)', fontsize = 30)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.title('E$_{1}$ vs $\\theta$', fontsize = 30)
plt.tick_params(axis='both', which='major', labelsize=20)
cbar = plt.colorbar()
plt.savefig(save_dir + 'POLARIS_Energy_vs_Theta.png')'''

'''# plot the above two graphs on the same figure with E1_vs_theta_NOCLF above and  POLARIS_Energy_vs_Theta below 
fig, ax = plt.subplots(2, 1, sharex=True, figsize = (30,15))
a = ax[0].scatter(coincidencePOLARISNOCLFplotE1['theta1_deg'], coincidencePOLARISNOCLFplotE1['E1'] , c =coincidencePOLARISNOCLFplotE1['counts'] ,vmin=1, vmax = 10, cmap = 'turbo', s=0.5)
ax[0].set_xlim(0, 180)
ax[0].set_ylim(0, 600)
ax[0].set_ylabel('Energy (keV)', fontsize = 30)
ax[0].set_title('E$_{1}$ vs $\\theta$', fontsize = 30)
ax[0].grid(axis='both', alpha=0.5, linestyle='--', linewidth=0.5)
ax[0].tick_params(axis='both', which='both', labelsize=20)
cbar = fig.colorbar(a, ax=ax[0])

b = ax[1].scatter(coincidencePOLARISplotE['theta1_deg'], coincidencePOLARISplotE['POLARIS_energy'] , c = (coincidencePOLARISplotE['counts']) ,vmin=1, vmax=2, cmap = 'turbo', s=0.5)
ax[1].set_xlim(0, 180)
ax[1].set_ylim(0, 600)
ax[1].set_title('E$_{1}$ vs $\\theta$ - CLF applied', fontsize = 30)
ax[0].tick_params(axis='x', which='both', labelbottom=False)
ax[1].set_ylabel('Energy (keV)', fontsize = 30)
ax[1].set_xlabel('$\\theta_{1}$ (deg)', fontsize = 30)
ax[1].grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
ax[1].tick_params(axis='both', which='both', labelsize=20)
cbar = fig.colorbar(b, ax=ax[1])
plt.savefig(save_dir + 'LowCBRange_E1_vs_theta_NOCLF' + '_and_' + 'POLARIS_Energy_vs_Theta' + '.png')

# plot the above two graphs on the same figure with E1_vs_theta_NOCLF above and  POLARIS_Energy_vs_Theta below 
fig, ax = plt.subplots(2, 1, sharex=True, figsize = (30,15))
a = ax[0].scatter(coincidencePOLARISNOCLFplotE1['theta1_deg'], coincidencePOLARISNOCLFplotE1['E1'] , c = coincidencePOLARISNOCLFplotE1['counts'] , cmap = 'turbo', s=0.5)
ax[0].set_xlim(0, 180)
ax[0].set_ylim(0, 600)
ax[0].set_ylabel('Energy (keV)', fontsize = 30)
ax[0].set_title('E$_{1}$ vs $\\theta$', fontsize = 30)
ax[0].grid(axis='both', alpha=0.5, linestyle='--', linewidth=0.5)
ax[0].tick_params(axis='both', which='both', labelsize=20)
cbar = fig.colorbar(a, ax=ax[0])
b = ax[1].scatter(coincidencePOLARISplotE['theta1_deg'], coincidencePOLARISplotE['POLARIS_energy'] , c = coincidencePOLARISplotE['counts'], cmap = 'turbo',vmin=1, s=0.5)
ax[1].set_xlim(0, 180)
ax[1].set_ylim(0, 600)
ax[1].set_title('E$_{1}$ vs $\\theta$ - CLF applied', fontsize = 30)
ax[0].tick_params(axis='x', which='both', labelbottom=False)
ax[1].set_ylabel('Energy (keV)', fontsize = 30)
ax[1].set_xlabel('$\\theta_{1}$ (deg)', fontsize = 30)
ax[1].grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
ax[1].tick_params(axis='both', which='both', labelsize=20)
cbar = fig.colorbar(b, ax=ax[1])
plt.savefig(save_dir + 'HighCBRange_E1_vs_theta_NOCLF' + '_and_' + 'POLARIS_Energy_vs_Theta' + '.png')

# plot E0 as a histogram from 0 to 5000 kev
plt.figure(figsize = (30, 15))
plt.hist(coincidence_no_CLF['E0'], bins = 5000, range = (0, 5000), color = 'royalblue', histtype='step', linewidth = 2)
plt.xlim(0, 600)
plt.xlabel('Reconstructed incident $\\gamma$ energy (keV)', fontsize = 30)
plt.ylabel('Counts', fontsize = 30)
plt.title('E$_{0}$', fontsize = 30)
plt.yscale('log')
plt.tick_params(axis='both', which='both', labelsize=20)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.savefig(save_dir + 'E0CCSpectrum' + '.png')


# plot E0 as a histogram from 0 to 5000 kev
plt.figure(figsize = (30, 15))
plt.hist(coincidence_no_CLF['E0'], bins = 5000, range = (0, 5000), color = 'royalblue', histtype='step', linewidth = 2)
plt.xlim(0, 5000)
plt.xlabel('Reconstructed incident $\\gamma$ energy (keV)', fontsize = 30)
plt.ylabel('Counts', fontsize = 30)
plt.title('E$_{0}$', fontsize = 30)
plt.tick_params(axis='both', which='both', labelsize=20)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.savefig(save_dir + 'E0CCSpectrum_no_log' + '.png')'''



'''coincidencePOLARISplotE = coincidence[(coincidence['POLARIS_energy'] > 0) & (coincidence['POLARIS_energy'] < 10000)]
plt.figure(figsize = (30, 15))
plt.scatter(coincidencePOLARISplotE['theta1_deg'], coincidencePOLARISplotE['POLARIS_energy'] ,  cmap = 'turbo',vmin=1, s=0.5)
plt.xlim(0, 180)
plt.ylim(0, 6000)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.xlabel('$\\theta$ (deg)', fontsize = 30)
plt.ylabel('Energy (keV)', fontsize = 30)
plt.title('E$_{1}$ vs Theta', fontsize = 30)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.savefig(save_dir + 'NoHeatBar_POLARIS_Energy_vs_Theta.png')
plt.close()'''


#--------------------------------------------------------------------------------


'''coincidencePOLARISNOCLFplotE0 = coincidence_no_CLF[(coincidence_no_CLF['E0'] > 0) & (coincidence_no_CLF['E0'] < 10000)]
coincidencePOLARISNOCLFplotE0 = coincidencePOLARISNOCLFplotE0.groupby(['E0', 'theta1_deg']).size().reset_index(name = 'counts')
plt.figure(figsize = (30, 15))
plt.scatter(coincidencePOLARISNOCLFplotE0['theta1_deg'], coincidencePOLARISNOCLFplotE0['E0'] , c = coincidencePOLARISNOCLFplotE0['counts'] , cmap = 'turbo',vmin=1, s=0.5)
plt.xlim(0, 180)
plt.ylim(0, 6000)
plt.ylabel('Energy  (keV)', fontsize = 30)
plt.xlabel('$\\theta$ (deg)', fontsize = 30)
plt.title('E$_{0}$ vs $\\theta$', fontsize = 30)
plt.tick_params(axis='both', which='major', labelsize=20)
cbar = plt.colorbar()
plt.grid(axis='both', alpha=0.5, linestyle='--', linewidth=0.5)
plt.savefig(save_dir + 'E0_vs_theta_NOCLF' + '.png')
print(' PLOT SAVED: ' + save_dir + 'E0_vs_theta_NOCLF' + '.png')'''

