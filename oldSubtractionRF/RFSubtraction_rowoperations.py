#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------
from cmath import nan
from ctypes import sizeof
from tokenize import Double
from unittest import skip
from venv import create
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
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
import multiprocessing as mp
from joblib import Parallel, delayed # for parallel processing
cores = 25 # number of cores to use

#------------------------------------------------------------------
data_dir = sys.argv[1]
ROOTFile = sys.argv[2] 

save_dir = data_dir + '/' + 'RFSubtraction/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))

f = uproot.open(data_dir + ROOTFile)
tree0 = f['LaBr0Data']
EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det0 = "TimingSlow_LaBr_Det0"
TimingSlow_RF = "TimingSlow_RF"
EnergySlow_RF = "EnergySlow_RF"
slowLaBr0calib1stFullRange = tree0[EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange].array()
slowtimeLaBr0 = tree0[TimingSlow_LaBr_Det0].array()
RFsignaltime0 = tree0[TimingSlow_RF].array()
RFsignalenergy0 = tree0[EnergySlow_RF].array()
df = pd.DataFrame({'slowLaBr0calib1stFullRange': slowLaBr0calib1stFullRange, 'slowtimeLaBr0': slowtimeLaBr0, 'RFsignaltime0': RFsignaltime0, 'RFsignalenergy0': RFsignalenergy0})

tree1 = f['LaBr1Data']
EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det1 = "TimingSlow_LaBr_Det1"
TimingSlow_RF = "TimingSlow_RF"
EnergySlow_RF = "EnergySlow_RF"
slowLaBr1calib1stFullRange = tree1[EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange].array()
slowtimeLaBr1 = tree1[TimingSlow_LaBr_Det1].array()
RFsignaltime1 = tree1[TimingSlow_RF].array()
RFsignalenergy1 = tree1[EnergySlow_RF].array()
df1 = pd.DataFrame({'slowLaBr1calib1stFullRange': slowLaBr1calib1stFullRange, 'slowtimeLaBr1': slowtimeLaBr1, 'RFsignaltime1': RFsignaltime1, 'RFsignalenergy1': RFsignalenergy1})

tree2 = f['LaBr2Data']
EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det2 = "TimingSlow_LaBr_Det2"
TimingSlow_RF = "TimingSlow_RF"
EnergySlow_RF = "EnergySlow_RF"
slowLaBr2calib1stFullRange = tree2[EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange].array()
slowtimeLaBr2 = tree2[TimingSlow_LaBr_Det2].array()
RFsignaltime2 = tree2[TimingSlow_RF].array()
RFsignalenergy2 = tree2[EnergySlow_RF].array()
df2 = pd.DataFrame({'slowLaBr2calib1stFullRange': slowLaBr2calib1stFullRange, 'slowtimeLaBr2': slowtimeLaBr2, 'RFsignaltime2': RFsignaltime2, 'RFsignalenergy2': RFsignalenergy2})

tree3 = f['LaBr3Data']
EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det3 = "TimingSlow_LaBr_Det3"
TimingSlow_RF = "TimingSlow_RF"
EnergySlow_RF = "EnergySlow_RF"
slowLaBr3calib1stFullRange = tree3[EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange].array()
slowtimeLaBr3 = tree3[TimingSlow_LaBr_Det3].array()
RFsignaltime3 = tree3[TimingSlow_RF].array()
RFsignalenergy3 = tree3[EnergySlow_RF].array()
df3 = pd.DataFrame({'slowLaBr3calib1stFullRange': slowLaBr3calib1stFullRange, 'slowtimeLaBr3': slowtimeLaBr3, 'RFsignaltime3': RFsignaltime3, 'RFsignalenergy3': RFsignalenergy3})

df['slowtimeLaBr0'] = df['slowtimeLaBr0']/10
df['RFsignaltime0'] = df['RFsignaltime0']/10
df = df.sort_values(by=['slowtimeLaBr0'])
df = df.reset_index(drop=True)
df1['slowtimeLaBr1'] = df1['slowtimeLaBr1']/10
df1['RFsignaltime1'] = df1['RFsignaltime1']/10
df1 = df1.sort_values(by=['slowtimeLaBr1'])
df1 = df1.reset_index(drop=True)
df2['slowtimeLaBr2'] = df2['slowtimeLaBr2']/10
df2['RFsignaltime2'] = df2['RFsignaltime2']/10
df2 = df2.sort_values(by=['slowtimeLaBr2'])
df2 = df2.reset_index(drop=True)
df3['slowtimeLaBr3'] = df3['slowtimeLaBr3']/10
df3['RFsignaltime3'] = df3['RFsignaltime3']/10
df3 = df3.sort_values(by=['slowtimeLaBr3'])
df3 = df3.reset_index(drop=True)


for i in range(len(df)):
    if df['slowtimeLaBr0'][i] > 1e13 and df['slowLaBr0calib1stFullRange'][i] != 0 and df['slowLaBr0calib1stFullRange'][i] != nan:
        firstindex = i
        break

# plot rf signal energy vs time as a heatmap. make the background white. Do this for the first 1min of data and bin in 10 ns bins

# first 1 min of data
df1 = df.loc[df['RFsignaltime0'] > 1e13, :]
df1 = df1.reset_index(drop=True)


#use root to plot a 2d histogram of rf energy vs time
'''c1 = root.TCanvas("c1","c1", 1000, 1000)
h1 = root.TH2F("h1", "h1", 60000, 0, 600000, 18000, 0, 18000)
for i in range(len(df1)):
    h1.Fill(df1['RFsignaltime0'][i], df1['RFsignalenergy0'][i])
h1.Draw("colz")
# label axes
h1.GetXaxis().SetTitle("RF signal time (ns)")
h1.GetYaxis().SetTitle("RF signal energy (keV)")
# set color palette
h1.SetContour(100)
c1.SetLogz()
c1.SaveAs(save_dir + "RFsignalenergy0vstime.png")
'''

timediffLaBr0 = df['slowtimeLaBr0'].diff()
timediffRF = df['RFsignaltime0'].diff()
#timediffLaBr0 = timediffLaBr0[firstindex:]
#timediffRF = timediffRF[firstindex:]
# find the mean timediffLaBr0 and timediffRF
timediffLaBr0mean = timediffLaBr0.mean()
timediffRFmean = timediffRF.mean()

plt.figure(figsize=(30,20))
mu, sigma = timediffLaBr0mean, timediffLaBr0.std()
n, bins, patches = plt.hist(timediffLaBr0, 80, density=True, facecolor='g', alpha=0.75)
y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
plt.plot(bins, y, '--')
plt.xlabel('Time difference between LaBr0 events (ns)')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ time\ difference\ between\ LaBr0\ events:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)
plt.xlim(0, 8)
plt.savefig(save_dir + 'timediffLaBr0.png')


# plot subplots of timediffLaBr0 and timediffRF and set their ranges to be the same
print ('Plotting time difference between consecutive RF events and time difference between consecutive slowtimeLaBr0 events')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30,20), sharex=True)
ax1.plot(timediffLaBr0, 'b-', label='LaBr$_3$')
ax2.plot(timediffRF, 'r-', label='RF')
ax1.legend( loc='upper right', fontsize=40)
ax2.legend( loc='upper right', fontsize=40)
ax1.set_xlim(0, 20)
ax2.set_xlim(0, 20)
ax2.set_xlabel('Time Difference (ns)', fontsize=40)
ax1.set_ylabel('Counts', fontsize=40)
ax2.set_ylabel('Counts', fontsize=40)
ax1.set_title('Time difference between consecutive slow time L0 events', fontsize=40)
ax2.set_title('Time difference between consecutive RF events', fontsize=40)
plt.savefig(save_dir + 'timediffLaBr0andRF.png')


print ('Plotting time difference between consecutive RF events and time difference between consecutive slowtimeLaBr0 events')
plt.figure(figsize=(30,20))
plt.hist(df['slowtimeLaBr0'], bins=60, range=(df['slowtimeLaBr0'][firstindex], df['slowtimeLaBr0'][firstindex]+600), histtype='step', label='slowtimeLaBr0', linewidth=2) 
plt.hist(df['RFsignaltime0'], bins=60, range=(df['slowtimeLaBr0'][firstindex], df['slowtimeLaBr0'][firstindex]+600), histtype='step', label='RFsignaltime0', linewidth=2)
plt.xlabel('Time (ns)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.legend( fontsize=40)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
plt.savefig(save_dir + 'slowtimeLaBr0_RFsignaltime0_600ns.png')
plt.close()


df['timedifflabrRF'] = df['slowtimeLaBr0']-df['RFsignaltime0']

# plot timedifflabrRF in 10 ns bins from 0 to 1000 ns
print ('Plotting time difference between RF and slowtimeLaBr0 events')
plt.figure(figsize=(30,20))
plt.hist(df['timedifflabrRF'], bins=np.arange(0, 1000, 10), color='blue', edgecolor='black', linewidth=1.2)
plt.xlabel('Time difference (ns)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.xticks(np.arange(0, 1000, 50), fontsize=20)
plt.title('Time difference (RF - LaBr$_3$) for events at same index', fontsize=40)
plt.savefig(save_dir + 'timedifflabrRFSameIndex_10nsbins.png')

df['boolinout'] = np.where((df['timedifflabrRF'] > -50) & (df['timedifflabrRF'] < 50), 1, 0) # 1 means in the range, 0 means out of the range

df['RFsignalshiftup1'] = df['RFsignaltime0'].shift(1)
df['timedifflabrRFshiftup1'] = df['slowtimeLaBr0']-df['RFsignalshiftup1']
df['boolinoutshiftup1'] = np.where((df['timedifflabrRFshiftup1'] > -50) & (df['timedifflabrRFshiftup1'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup2'] = df['RFsignaltime0'].shift(2)
df['timedifflabrRFshiftup2'] = df['slowtimeLaBr0']-df['RFsignalshiftup2']
df['boolinoutshiftup2'] = np.where((df['timedifflabrRFshiftup2'] > -50) & (df['timedifflabrRFshiftup2'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup3'] = df['RFsignaltime0'].shift(3)
df['timedifflabrRFshiftup3'] = df['slowtimeLaBr0']-df['RFsignalshiftup3']
df['boolinoutshiftup3'] = np.where((df['timedifflabrRFshiftup3'] > -50) & (df['timedifflabrRFshiftup3'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup4'] = df['RFsignaltime0'].shift(4)
df['timedifflabrRFshiftup4'] = df['slowtimeLaBr0']-df['RFsignalshiftup4']
df['boolinoutshiftup4'] = np.where((df['timedifflabrRFshiftup4'] > -50) & (df['timedifflabrRFshiftup4'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup5'] = df['RFsignaltime0'].shift(5)
df['timedifflabrRFshiftup5'] = df['slowtimeLaBr0']-df['RFsignalshiftup5']
df['boolinoutshiftup5'] = np.where((df['timedifflabrRFshiftup5'] > -50) & (df['timedifflabrRFshiftup5'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup6'] = df['RFsignaltime0'].shift(6)
df['timedifflabrRFshiftup6'] = df['slowtimeLaBr0']-df['RFsignalshiftup6']
df['boolinoutshiftup6'] = np.where((df['timedifflabrRFshiftup6'] > -50) & (df['timedifflabrRFshiftup6'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup7'] = df['RFsignaltime0'].shift(7)
df['timedifflabrRFshiftup7'] = df['slowtimeLaBr0']-df['RFsignalshiftup7']
df['boolinoutshiftup7'] = np.where((df['timedifflabrRFshiftup7'] > -50) & (df['timedifflabrRFshiftup7'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup8'] = df['RFsignaltime0'].shift(8)
df['timedifflabrRFshiftup8'] = df['slowtimeLaBr0']-df['RFsignalshiftup8']
df['boolinoutshiftup8'] = np.where((df['timedifflabrRFshiftup8'] > -50) & (df['timedifflabrRFshiftup8'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup9'] = df['RFsignaltime0'].shift(9)
df['timedifflabrRFshiftup9'] = df['slowtimeLaBr0']-df['RFsignalshiftup9']
df['boolinoutshiftup9'] = np.where((df['timedifflabrRFshiftup9'] > -50) & (df['timedifflabrRFshiftup9'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftup10'] = df['RFsignaltime0'].shift(10)
df['timedifflabrRFshiftup10'] = df['slowtimeLaBr0']-df['RFsignalshiftup10']
df['boolinoutshiftup10'] = np.where((df['timedifflabrRFshiftup10'] > -50) & (df['timedifflabrRFshiftup10'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
print('number of boolinoutshiftup10 == 1: ', len(df[df['boolinoutshiftup10'] == 1]))

df['RFsignalshiftdown1'] = df['RFsignaltime0'].shift(-1)
df['timedifflabrRFshiftdown1'] = df['slowtimeLaBr0']-df['RFsignalshiftdown1']
df['boolinoutshiftdown1'] = np.where((df['timedifflabrRFshiftdown1'] > -50) & (df['timedifflabrRFshiftdown1'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown2'] = df['RFsignaltime0'].shift(-2)
df['timedifflabrRFshiftdown2'] = df['slowtimeLaBr0']-df['RFsignalshiftdown2']
df['boolinoutshiftdown2'] = np.where((df['timedifflabrRFshiftdown2'] > -50) & (df['timedifflabrRFshiftdown2'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown3'] = df['RFsignaltime0'].shift(-3)
df['timedifflabrRFshiftdown3'] = df['slowtimeLaBr0']-df['RFsignalshiftdown3']
df['boolinoutshiftdown3'] = np.where((df['timedifflabrRFshiftdown3'] > -50) & (df['timedifflabrRFshiftdown3'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown4'] = df['RFsignaltime0'].shift(-4)
df['timedifflabrRFshiftdown4'] = df['slowtimeLaBr0']-df['RFsignalshiftdown4']
df['boolinoutshiftdown4'] = np.where((df['timedifflabrRFshiftdown4'] > -50) & (df['timedifflabrRFshiftdown4'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown5'] = df['RFsignaltime0'].shift(-5)
df['timedifflabrRFshiftdown5'] = df['slowtimeLaBr0']-df['RFsignalshiftdown5']
df['boolinoutshiftdown5'] = np.where((df['timedifflabrRFshiftdown5'] > -50) & (df['timedifflabrRFshiftdown5'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown6'] = df['RFsignaltime0'].shift(-6)
df['timedifflabrRFshiftdown6'] = df['slowtimeLaBr0']-df['RFsignalshiftdown6']
df['boolinoutshiftdown6'] = np.where((df['timedifflabrRFshiftdown6'] > -50) & (df['timedifflabrRFshiftdown6'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown7'] = df['RFsignaltime0'].shift(-7) 
df['timedifflabrRFshiftdown7'] = df['slowtimeLaBr0']-df['RFsignalshiftdown7']
df['boolinoutshiftdown7'] = np.where((df['timedifflabrRFshiftdown7'] > -50) & (df['timedifflabrRFshiftdown7'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown8'] = df['RFsignaltime0'].shift(-8)
df['timedifflabrRFshiftdown8'] = df['slowtimeLaBr0']-df['RFsignalshiftdown8']
df['boolinoutshiftdown8'] = np.where((df['timedifflabrRFshiftdown8'] > -50) & (df['timedifflabrRFshiftdown8'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown9'] = df['RFsignaltime0'].shift(-9)
df['timedifflabrRFshiftdown9'] = df['slowtimeLaBr0']-df['RFsignalshiftdown9']
df['boolinoutshiftdown9'] = np.where((df['timedifflabrRFshiftdown9'] > -50) & (df['timedifflabrRFshiftdown9'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
df['RFsignalshiftdown10'] = df['RFsignaltime0'].shift(-10)
df['timedifflabrRFshiftdown10'] = df['slowtimeLaBr0']-df['RFsignalshiftdown10']
df['boolinoutshiftdown10'] = np.where((df['timedifflabrRFshiftdown10'] > -50) & (df['timedifflabrRFshiftdown10'] < 50), 1, 0) # 1 means in the range, 0 means out of the range
print('number of boolinoutshiftdown10 == 1: ', len(df[df['boolinoutshiftdown10'] == 1]))

# do the same plot as timedifflabrRFSameIndex but with the shifted RF signal. plot the shiftup timediff histograms on the same plot with no alpha
plt.figure(figsize=(30, 20))
plt.hist(df['timedifflabrRFshiftup1'], bins=np.arange(0, 50, 10),  label='shifted up 1 row', color='red')
plt.hist(df['timedifflabrRFshiftup2'], bins=np.arange(0, 50, 10),  label='shifted up 2 rows', color='orange')
plt.hist(df['timedifflabrRFshiftup3'], bins=np.arange(0, 50, 10),  label='shifted up 3 rows', color='blue')
plt.hist(df['timedifflabrRFshiftup4'], bins=np.arange(0, 50, 10),  label='shifted up 4 rows', color='green')
plt.hist(df['timedifflabrRFshiftup5'], bins=np.arange(0, 50, 10),  label='shifted up 5 rows', color='purple')
plt.hist(df['timedifflabrRFshiftup6'], bins=np.arange(0, 50, 10),  label='shifted up 6 rows', color='pink')
plt.hist(df['timedifflabrRFshiftup7'], bins=np.arange(0, 50, 10),  label='shifted up 7 rows', color='black')
plt.hist(df['timedifflabrRFshiftup8'], bins=np.arange(0, 50, 10),  label='shifted up 8 rows', color='lightblue')
plt.hist(df['timedifflabrRFshiftup9'], bins=np.arange(0, 50, 10),  label='shifted up 9 rows', color='lightgreen')
plt.hist(df['timedifflabrRFshiftup10'], bins=np.arange(0, 50, 10),  label='shifted up 10 rows', color='brown')
plt.xlabel('Time difference (ns)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.xticks(np.arange(0, 50, 5))
plt.legend(fontsize=30)
plt.title('Time difference (RF - LaBr$_3$) for events at same index', fontsize=40)
plt.savefig(save_dir + 'timedifflabrRFShiftedUp_10nsbins.png')

# do the same plot as timedifflabrRFSameIndex but with the shifted RF signal. plot the shiftdown timediff histograms on the same plot with no alpha
plt.figure(figsize=(30, 20))
plt.hist(df['timedifflabrRFshiftdown1'], bins=np.arange(0, 50, 10),  label='shifted down 1 row' , color='red')
plt.hist(df['timedifflabrRFshiftdown2'], bins=np.arange(0, 50, 10),  label='shifted down 2 rows' , color='orange')
plt.hist(df['timedifflabrRFshiftdown3'], bins=np.arange(0, 50, 10),  label='shifted down 3 rows' , color='blue')
plt.hist(df['timedifflabrRFshiftdown4'], bins=np.arange(0, 50, 10),  label='shifted down 4 rows' , color='green')
plt.hist(df['timedifflabrRFshiftdown5'], bins=np.arange(0, 50, 10),  label='shifted down 5 rows' , color='purple')
plt.hist(df['timedifflabrRFshiftdown6'], bins=np.arange(0, 50, 10),  label='shifted down 6 rows' , color='pink')
plt.hist(df['timedifflabrRFshiftdown7'], bins=np.arange(0, 50, 10),  label='shifted down 7 rows' , color='black')
plt.hist(df['timedifflabrRFshiftdown8'], bins=np.arange(0, 50, 10),  label='shifted down 8 rows' , color='lightblue')
plt.hist(df['timedifflabrRFshiftdown9'], bins=np.arange(0, 50, 10),  label='shifted down 9 rows' , color='lightgreen')
plt.hist(df['timedifflabrRFshiftdown10'], bins=np.arange(0, 50, 10),  label='shifted down 10 rows' , color='brown')
plt.xlabel('Time difference (ns)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.legend(fontsize=30)
plt.xticks(np.arange(0, 50, 5))
plt.title('Time difference (RF - LaBr$_3$) for events at same index', fontsize=40)
plt.savefig(save_dir + 'timedifflabrRFShiftedDown_10nsbins.png')


# pull out all df['slowLaBr0calib1stFullRange'] that have a corresponding RF signal within 50 ns
insync = df['slowLaBr0calib1stFullRange'][(df['boolinout'] == 1) | (df['boolinoutshiftup1'] == 1) | (df['boolinoutshiftup2'] == 1) | (df['boolinoutshiftup3'] == 1) | (df['boolinoutshiftup4'] == 1) | (df['boolinoutshiftup5'] == 1) | (df['boolinoutshiftup6'] == 1) | (df['boolinoutshiftup7'] == 1) | (df['boolinoutshiftup8'] == 1) | (df['boolinoutshiftup9'] == 1) | (df['boolinoutshiftup10'] == 1) | (df['boolinoutshiftdown1'] == 1) | (df['boolinoutshiftdown2'] == 1) | (df['boolinoutshiftdown3'] == 1) | (df['boolinoutshiftdown4'] == 1) | (df['boolinoutshiftdown5'] == 1) | (df['boolinoutshiftdown6'] == 1) | (df['boolinoutshiftdown7'] == 1) | (df['boolinoutshiftdown8'] == 1) | (df['boolinoutshiftdown9'] == 1) | (df['boolinoutshiftdown10'] == 1)]

outofsync = df['slowLaBr0calib1stFullRange'][(df['boolinout'] == 0) | (df['boolinoutshiftup1'] == 0) | (df['boolinoutshiftup2'] == 0) | (df['boolinoutshiftup3'] == 0) | (df['boolinoutshiftup4'] == 0) | (df['boolinoutshiftup5'] == 0) | (df['boolinoutshiftup6'] == 0) | (df['boolinoutshiftup7'] == 0) | (df['boolinoutshiftup8'] == 0) | (df['boolinoutshiftup9'] == 0) | (df['boolinoutshiftup10'] == 0) | (df['boolinoutshiftdown1'] == 0) | (df['boolinoutshiftdown2'] == 0) | (df['boolinoutshiftdown3'] == 0) | (df['boolinoutshiftdown4'] == 0) | (df['boolinoutshiftdown5'] == 0) | (df['boolinoutshiftdown6'] == 0) | (df['boolinoutshiftdown7'] == 0) | (df['boolinoutshiftdown8'] == 0) | (df['boolinoutshiftdown9'] == 0) | (df['boolinoutshiftdown10'] == 0)]

print('RF signals in sync with LaBr0: \n', df['RFsignaltime0'][(df['boolinout'] == 1) |(df['boolinoutshiftup1'] == 1) | (df['boolinoutshiftup2'] == 1) | (df['boolinoutshiftup3'] == 1) | (df['boolinoutshiftup4'] == 1) | (df['boolinoutshiftup5'] == 1) | (df['boolinoutshiftup6'] == 1) | (df['boolinoutshiftup7'] == 1) | (df['boolinoutshiftup8'] == 1) | (df['boolinoutshiftup9'] == 1) | (df['boolinoutshiftup10'] == 1) | (df['boolinoutshiftdown1'] == 1) | (df['boolinoutshiftdown2'] == 1) | (df['boolinoutshiftdown3'] == 1) | (df['boolinoutshiftdown4'] == 1) | (df['boolinoutshiftdown5'] == 1) | (df['boolinoutshiftdown6'] == 1) | (df['boolinoutshiftdown7'] == 1) | (df['boolinoutshiftdown8'] == 1) | (df['boolinoutshiftdown9'] == 1) | (df['boolinoutshiftdown10'] == 1)])
print('RF signals out of sync with LaBr0: \n', df['RFsignaltime0'][(df['boolinout'] == 0) |(df['boolinoutshiftup1'] == 0) | (df['boolinoutshiftup2'] == 0) | (df['boolinoutshiftup3'] == 0) | (df['boolinoutshiftup4'] == 0) | (df['boolinoutshiftup5'] == 0) | (df['boolinoutshiftup6'] == 0) | (df['boolinoutshiftup7'] == 0) | (df['boolinoutshiftup8'] == 0) | (df['boolinoutshiftup9'] == 0) | (df['boolinoutshiftup10'] == 0) | (df['boolinoutshiftdown1'] == 0) | (df['boolinoutshiftdown2'] == 0) | (df['boolinoutshiftdown3'] == 0) | (df['boolinoutshiftdown4'] == 0) | (df['boolinoutshiftdown5'] == 0) | (df['boolinoutshiftdown6'] == 0) | (df['boolinoutshiftdown7'] == 0) | (df['boolinoutshiftdown8'] == 0) | (df['boolinoutshiftdown9'] == 0) | (df['boolinoutshiftdown10'] == 0)])

timeRFinsync = df['RFsignaltime0'][(df['boolinout'] == 1) | (df['boolinoutshiftup1'] == 1) | (df['boolinoutshiftup2'] == 1) | (df['boolinoutshiftup3'] == 1) | (df['boolinoutshiftup4'] == 1) | (df['boolinoutshiftup5'] == 1) | (df['boolinoutshiftup6'] == 1) | (df['boolinoutshiftup7'] == 1) | (df['boolinoutshiftup8'] == 1) | (df['boolinoutshiftup9'] == 1) | (df['boolinoutshiftup10'] == 1) | (df['boolinoutshiftdown1'] == 1) | (df['boolinoutshiftdown2'] == 1) | (df['boolinoutshiftdown3'] == 1) | (df['boolinoutshiftdown4'] == 1) | (df['boolinoutshiftdown5'] == 1) | (df['boolinoutshiftdown6'] == 1) | (df['boolinoutshiftdown7'] == 1) | (df['boolinoutshiftdown8'] == 1) | (df['boolinoutshiftdown9'] == 1) | (df['boolinoutshiftdown10'] == 1)]
timeRFoutofsync = df['RFsignaltime0'][(df['boolinout'] == 0) | (df['boolinoutshiftup1'] == 0) | (df['boolinout'] == 0) |(df['boolinoutshiftup2'] == 0) | (df['boolinoutshiftup3'] == 0) | (df['boolinoutshiftup4'] == 0) | (df['boolinoutshiftup5'] == 0) | (df['boolinoutshiftup6'] == 0) | (df['boolinoutshiftup7'] == 0) | (df['boolinoutshiftup8'] == 0) | (df['boolinoutshiftup9'] == 0) | (df['boolinoutshiftup10'] == 0) | (df['boolinoutshiftdown1'] == 0) | (df['boolinoutshiftdown2'] == 0) | (df['boolinoutshiftdown3'] == 0) | (df['boolinoutshiftdown4'] == 0) | (df['boolinoutshiftdown5'] == 0) | (df['boolinoutshiftdown6'] == 0) | (df['boolinoutshiftdown7'] == 0) | (df['boolinoutshiftdown8'] == 0) | (df['boolinoutshiftdown9'] == 0) | (df['boolinoutshiftdown10'] == 0)]
timeRFinsync = timeRFinsync.reset_index(drop=True)
timeRFoutofsync = timeRFoutofsync.reset_index(drop=True)

minRFinsync = timeRFinsync[timeRFinsync > 1e14].min()
maxRFinsync = max(timeRFinsync)
minRFoutofsync = timeRFoutofsync[timeRFoutofsync > 1e14].min() 
maxRFoutofsync = max(timeRFoutofsync)


b = int((maxRFinsync - minRFinsync)/100000)

print( 'minRF = ', minRFinsync, 'maxRF = ', maxRFinsync, 'bins = ', b)

print ('Plotting histogram of RF signals in sync with LaBr0')
plt.figure(figsize=(30,20))
m, bins = np.histogram(timeRFinsync, bins = b, range = (minRFinsync, maxRFinsync))
n, bins = np.histogram(timeRFoutofsync, bins = b, range = (minRFinsync, maxRFinsync))
plt.plot(bins[:-1], m, label = 'RF signals in sync with LaBr0')
plt.plot(bins[:-1], n, label = 'RF signals out of sync with LaBr0')
plt.xlabel('RF signal time (ns)', fontsize = 30)
plt.ylabel('Counts', fontsize = 30)
plt.legend(fontsize=40)
plt.title('RF signal time', fontsize = 30)
plt.savefig(save_dir + 'RFsignalsInOutSync_Time.png')


# plot the timeRFinsync vs timeRFoutofsync histograms between minRF and maxRF in bins
print ('Plotting histogram of RF signals in sync with LaBr0 vs RF signals out of sync with LaBr0')
plt.figure(figsize=(30,20))
plt.hist(timeRFinsync , bins = b, range = (minRFinsync, maxRFinsync), histtype='step', label = 'RF signals in sync with LaBr0')
plt.hist(timeRFoutofsync , bins = b, range = (minRFinsync, maxRFinsync), histtype='step', label = 'RF signals out of sync with LaBr0')
plt.legend(loc = 'upper right', fontsize = 30)
plt.xlabel('RF signal time (ns)', fontsize = 30)
plt.ylabel('Counts', fontsize = 30)
plt.title('RF signal time for events in sync and out of sync with LaBr$_3$ detector', fontsize = 30)
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.savefig(save_dir + 'InOutSync_10UpDown_RF.png')


# plot the timeRFinsync vs timeRFoutofsync histograms between minRF and maxRF in bins
print ('Plotting histogram of RF signals in sync with LaBr0 vs RF signals out of sync with LaBr0')
plt.figure(figsize=(30,20))
plt.hist(insync, bins=6000, range=(0, 6000), histtype='step', label='in sync')
plt.hist(outofsync, bins=6000, range=(0, 6000), histtype='step', label='out of sync')
plt.xlabel('Energy (keV)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.legend( fontsize=40)
plt.xticks(fontsize=40)
plt.yscale('log')
plt.yticks(fontsize=40)
plt.savefig(save_dir + 'pm50ns_insyncenergy_outsyncenergy_10rowsaboveandbelow.png')

print ('Plotting histogram of RF signals in sync with LaBr0 vs RF signals out of sync with LaBr0')
plt.figure(figsize=(30,20))
plt.hist(insync, bins=6000, range=(0, 6000), histtype='step', label='in sync')
plt.hist(outofsync, bins=6000, range=(0, 6000), histtype='step', label='out of sync')
plt.xlabel('Energy (keV)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.legend( fontsize=40)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
plt.savefig(save_dir + 'pm50ns_nolog_insyncenergy_outsyncenergy_10rowsaboveandbelow.png')

# use root to plot the df['slowLaBr0calib1stFullRange'] and insync on the same plot
c = root.TCanvas('c', 'c', 800, 600)
c.SetLogy()
h1 = root.TH1F('h1', 'h1', 6000, 0, 6000)
h2 = root.TH1F('h2', 'h2', 6000, 0, 6000)
h3 = root.TH1F('h3', 'h3', 6000, 0, 6000)
for i in range(len(insync)):
    h1.Fill(insync[i])
for i in range(len(df['slowLaBr0calib1stFullRange'])):
    h2.Fill(df['slowLaBr0calib1stFullRange'][i])
h1.SetLineColor(2)
h1.SetLineWidth(2)
h2.SetLineColor(4)
h2.SetLineWidth(2)
h1.Draw()
h2.Draw('same')
h3 = h2.Clone()
h3.SetLineColor(1)
h3.SetLineWidth(2)
h3.Add(h1, -1)
h3.Draw('same')
c.BuildLegend()
c.SaveAs(save_dir + 'pm50ns_insyncenergy_outsyncenergy_10rowsaboveandbelow.root')