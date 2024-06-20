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
'''import multiprocessing as mp
from joblib import Parallel, delayed # for parallel processing
cores = 25 # number of cores to use
'''
#------------------------------------------------------------------
data_dir = sys.argv[1]
ROOTFile = sys.argv[2] 

save_dir = data_dir + '/' + 'RFSubtraction/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))

f = uproot.open(data_dir + ROOTFile)
tree = f['LaBr1Data']
EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det0 = "TimingSlow_LaBr_Det1"
TimingSlow_RF = "TimingSlow_RF"
EnergySlow_RF = "EnergySlow_RF"
slowLaBr0calib1stFullRange = tree[EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange].array()
slowtimeLaBr0 = tree[TimingSlow_LaBr_Det0].array()
RFsignaltime = tree[TimingSlow_RF].array()
RFsignalenergy = tree[EnergySlow_RF].array()
df = pd.DataFrame({'slowLaBr0calib1stFullRange': slowLaBr0calib1stFullRange, 'slowtimeLaBr0': slowtimeLaBr0, 'RFsignaltime': RFsignaltime, 'RFsignalenergy': RFsignalenergy})

df['slowtimeLaBr0'] = df['slowtimeLaBr0']/10
df['RFsignaltime'] = df['RFsignaltime']/10
df = df.sort_values(by=['slowtimeLaBr0'])
df = df.reset_index(drop=True)

for i in range(len(df)):
    if df['slowtimeLaBr0'][i] > 1e13 and df['slowLaBr0calib1stFullRange'][i] != 0 and df['slowLaBr0calib1stFullRange'][i] != nan:
        firstindex = i
        time1 = df['slowtimeLaBr0'][i]
        break

for i in range(len(df)):
    if df['slowtimeLaBr0'][i] > 1e13 and df['slowLaBr0calib1stFullRange'][i] != 0 and df['slowLaBr0calib1stFullRange'][i] != nan:
        lastindex = i
        time2 = df['slowtimeLaBr0'][i]
        break
    
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30,20), sharex=True)
ax1.plot(df['slowtimeLaBr0'][firstindex:], 'b-', label='slowtimeLaBr0')
ax2.plot(df['RFsignaltime'][firstindex:], 'r-', label='RFsignaltime')
ax1.set_ylabel('Counts', fontsize=40)
ax2.set_ylabel('Counts', fontsize=40)
ax2.set_xlabel('Time (ns)', fontsize=40)
ax1.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', which='major', labelsize=30)
ax1.legend(fontsize=40)
ax2.legend(fontsize=40)
plt.savefig(save_dir + 'RFLaBrSignalTimes.png')
plt.close()

timediffLaBr0 = df['slowtimeLaBr0'].diff()
mintime = timediffLaBr0.min()

plt.figure(figsize = (30,15))
plt.hist(timediffLaBr0, bins = 1000, histtype = 'step', label = 'POLARIS', range=(mintime, 2), color = 'royalblue', linewidth = 2)
plt.xlabel('Time (ns)', fontsize = 40)
plt.ylabel('Counts', fontsize = 40)
plt.tick_params(axis='both', which='major', labelsize=30)
plt.legend(fontsize=40)
plt.savefig(save_dir + 'TimeDiffLaBr0.png')

plt.figure(figsize=(30,20))
plt.hist(df['slowtimeLaBr0'], bins=80000, range=(time1,time1+80000), histtype='step', color='b', label='slowtimeLaBr0')
plt.xlabel('Time (ns)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.title('Time of LaBr0 events', fontsize=40)
plt.tick_params(axis='both', which='major', labelsize=30)
plt.savefig(save_dir + 'LaBrTime.png')

plt.figure(figsize=(30,20))
plt.hist(df['RFsignaltime'], bins=800000, range=(time1,time1+800000), histtype='step', color='r', label='RFsignaltime')
plt.xlabel('Time (ns)', fontsize=40)
plt.ylabel('Counts', fontsize=40)
plt.title('Time of RF events', fontsize=40)
plt.tick_params(axis='both', which='major', labelsize=30)
plt.savefig(save_dir + 'RFTime.png')



df['TimeDiff'] = abs(df['RFsignaltime'] - df['slowtimeLaBr0']) 
insync = df.loc[(df['TimeDiff'] <=70, 'slowLaBr0calib1stFullRange')]

for i in range(20):
    df['RFsignalshiftup{}'.format(i+1)] = df['RFsignaltime'].shift(i+1)
    df['timedifflabrRFshiftup{}'.format(i+1)] = df['slowtimeLaBr0']-df['RFsignalshiftup{}'.format(i+1)].abs()
    df['RFsignalshiftdown{}'.format(i+1)] = df['RFsignaltime'].shift(-(i+1))
    df['timedifflabrRFshiftdown{}'.format(i+1)] = df['slowtimeLaBr0']-df['RFsignalshiftdown{}'.format(i+1)].abs()
    insync = insync.append(df.loc[(df['timedifflabrRFshiftup{}'.format(i+1)] <=70, 'slowLaBr0calib1stFullRange')])
    insync = insync.append(df.loc[(df['timedifflabrRFshiftdown{}'.format(i+1)] <=70, 'slowLaBr0calib1stFullRange')])


# plot the insync data as a histogram in root

h1 = root.TH1F("h1", "InSync", 8000, 0, 8000)
h2 = root.TH1F("h2", "All ", 8000, 0, 8000)
h3 = root.TH1F("h3", "OutSync", 8000, 0, 8000)
h4 = root.TH1F("h4", "Prompt", 8000, 0, 8000)

for i in range(len(insync)):
    h1.Fill(insync[i])
h1.Draw()
h1.Write()


for i in range(len(df['slowLaBr0calib1stFullRange'])):
    h2.Fill(df['slowLaBr0calib1stFullRange'][i])
h2.Draw()
h2.Write()
# plot c2 - c1 using Add() and Draw()
c3 = root.TCanvas("c3","c3",1000,800)
h3 = h2.Clone()
h3.Add(h1, -1)
h3.SetLineColor(root.kRed)
h3.Draw()
h2.SetLineColor(root.kBlue)
h2.Draw("same")
h1.SetLineColor(root.kBlack)
h1.Draw("same")
l = root.TLegend(0.7,0.7,0.9,0.9)
l.AddEntry(h1, "In Sync Data", "l")
l.AddEntry(h2, "All Data", "l")
l.AddEntry(h3, "Out of Sync Data", "l")
l.Draw()
c3.SetLogy()
c3.SaveAs(save_dir + 'SubtractedData.root')
c3.SaveAs(save_dir + 'SubtractedData.png')

h3.Write()
# plot h1 - h3 using Add() and Draw()
c4 = root.TCanvas("c4","c4",1000,800)
h4 = h1.Clone()
h4.Add(h3, -1)
h4.Draw()
h1.Draw("same")
h3.Draw("same")
h3.SetLineColor(root.kRed)
h1.SetLineColor(root.kBlack)
h4.SetLineColor(root.kMagenta)
l = root.TLegend(0.7,0.7,0.9,0.9)
l.AddEntry(h1, "In Sync Data", "l")
l.AddEntry(h3, "Out of Sync Data", "l")
l.AddEntry(h4, "In Sync Data - Out of Sync Data", "l")
l.Draw()
c4.SetLogy()
c4.SaveAs(save_dir + 'SubtractedData2.root')
c4.SaveAs(save_dir + 'SubtractedData2.png')