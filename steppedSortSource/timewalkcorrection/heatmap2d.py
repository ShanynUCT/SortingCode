import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import ROOT as ROOT
import os
import seaborn as sns
import pandas as pd


plt.rcParams.update({'font.size': 18})
# make axis labels bold
plt.rcParams["axes.labelweight"] = "bold"
# label size
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)

# LaBr0 slow data 2 inch
""" risetime = [0.1,0.3,0.5,0.7,0.9,1.1,1.3,0.1,0.3,0.5,0.7,0.9,1.1,1.3,0.1,0.3,0.5,0.7,0.9,1.1,1.3,0.1,0.3,0.5,0.7,0.9,1.1,1.3,0.1,0.3,0.5,0.7,0.9,1.1,1.3,0.1,0.3,0.5,0.7,0.9,1.1,1.3,0.1,0.3,0.5,0.7,0.9,1.1,1.3]
flattop = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.3,1.3,1.3,1.3,1.3,1.3,1.3]
fwhm = [4.842072263,4.670767933,4.678504721,4.818509871,4.13570544,4.218403548,4.207315597,9.729257421,4.096109293,4.320639856,4.35512696,4.397651442,4.109715233,4.219867638,5.811796862,5.175578848,4.109363806,4.192081763,4.231223429,4.076850438,4.113794048,4.893864215,4.771275146,5.023747564,4.143504514,4.376423075,4.228794341,4.338986238,4.138489434,4.331958703,4.503138444,4.776084131,4.072459975,4.216067769,4.270145349,6.496681802,4.18298987,4.237316283,4.419273432,4.429455532,4.094702776,4.319307187,5.492929602,5.897336573,4.172581799,4.320114783,4.300154792,5.798990255,8.117832246]
 """
# LaBr0 fast data 2 inch
""" risetime = [0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
flattop = [0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
fwhm = [5.159916746,4.811093988,5.104379069,4.440307473,5.326172313,4.591222848,5.395107142,5.011440789,5.408073747,4.540947506,5.18731261,4.503821616,5.315133944,4.317939001,4.979643325,5.102338548,5.120671093,4.268220905,5.123475306,4.453958867,5.391531558,5.240001849,5.346815146,4.233552764,5.326498774,4.778041439,4.749719731,5.276438248,5.219700696,4.973029447,5.325418626,4.442850907,5.387447047,4.891174943,5.573296106,4.975371736,4.838432233,4.526518255,5.454565834,4.492508181,5.319706276,5.099871786,5.49119246,4.89828547,5.416075504,4.538128424,5.484313725,4.516728678,5.399482694]
 """
# LaBr Assembly 05 slow data
""" risetime = [0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7,0.5,0.7,0.9,1.1,1.3,1.5,1.7]
flattop = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.3,1.3,1.3,1.3,1.3,1.3,1.3,0.06,0.06,0.06,0.06,0.06,0.06,0.06]
fwhm = [3.37742277,3.536883679,4.007363353,4.012699764,4.062065771,4.127602976,4.343850786,4.229192356,4.432847533,4.220038734,4.450411636,4.224158165,4.203329168,4.479616307,4.580012453,4.547309616,4.525003114,4.390426943,4.409294083,4.629235942,3.975691473,4.524340793,4.648310035,4.702237449,4.705581759,4.896343745,4.598371315,4.703056079,4.930543468,4.764686373,4.781758443,4.293128359,4.473211022,5.020488292,3.861096572,4.944049631,4.960828057,4.442110233,4.116342618,4.307849273,4.015588347,3.757931239,4.725238295,4.733327577,4.904384618,4.464956928,3.820583711,4.249403204,4.735084238,3.302489909,3.905936636,3.90132106,4.198917584,4.246437541,4.176753061,4.130262599] 
 """

# LaBr Assembly 05 slow data triangular fit
""" risetime = [2.2,2,2.4,1.8,1.6,1.4,1.2,1,0.88,0.6]
peaksample = [2.2,2,2.4,1.8,1.6,1.4,1.2,1,0.88,0.6]
fwhm = [3.942,4.02,4.04,3.77,3.85,3.6,3.552,3.36,3.4,3.85]
 """

# SrI2 Assembly 04 slow data triangular fit
risetime = [10,10.16,12,14.08,20,12.96,15.04,19.04,17.92,15.84,16.32,16.96]
peaksample = [10,10.16,12,14.08,20,12.96,15.04,19.04,17.92,15.84,16.32,16.96]
fwhm = [4.444,4.27,4.03507,4.24,4.0506,4.012,3.923,4.47,3.8401,4.04,4,4.3]

data = pd.DataFrame({'Risetime': risetime, 'Peaksample': peaksample, 'FWHM': fwhm})
""" duplicates = data.duplicated(subset=['Risetime', 'Peaksample'], keep=False)
print(duplicates) """
# Pivot the DataFrame to get it in the right format for a heatmap, sorting 'Flattop' in ascending order
data_pivot = data.pivot('Risetime', 'Peaksample', 'FWHM').sort_index(ascending=True)

# Create the heatmap
plt.figure( figsize=(10,6))
sns.heatmap(data_pivot,  fmt=".2f", cmap='viridis', cbar_kws={'label': 'Energy Resolution Percentage at 662 keV'}, norm=PowerNorm(gamma=0.45))
plt.gca().invert_yaxis()
plt.grid( linestyle='--', linewidth=0.5, alpha=0.5)
plt.xlabel('Trapezoid rise time (µs)')
plt.ylabel('Trapezoid peak sample (µs)')
#make a box around the plot
plt.gca().add_patch(plt.Rectangle((0, 0), 12, 12, fill=None, edgecolor='black', lw=2))

plt.show()