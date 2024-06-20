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

#------------------------------------------------------------------
data_dir = sys.argv[1]
ROOTFile = sys.argv[2] 

save_dir = data_dir + '/' + 'RFSubtraction/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))

f = uproot.open(data_dir + ROOTFile)
tree = f['LaBr0Data']
EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det0 = "TimingSlow_LaBr_Det0"
TimingSlow_RF = "TimingSlow_RF"
slowLaBr0calib1stFullRange = tree[EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange].array()
slowtimeLaBr0 = tree[TimingSlow_LaBr_Det0].array()
RFsignaltime = tree[TimingSlow_RF].array()
df = pd.DataFrame({'slowLaBr0calib1stFullRange': slowLaBr0calib1stFullRange, 'slowtimeLaBr0': slowtimeLaBr0, 'RFsignaltime': RFsignaltime})

df['slowtimeLaBr0'] = df['slowtimeLaBr0']/10
df['RFsignaltime'] = df['RFsignaltime']/10
df = df.sort_values(by=['slowtimeLaBr0'])
df = df.reset_index(drop=True)

for i in range(len(df)):
    if df['slowtimeLaBr0'][i] > 1e13 and df['slowLaBr0calib1stFullRange'][i] != 0 and df['slowLaBr0calib1stFullRange'][i] != nan:
        firstindex = i
        break

# for every RFsignaltimes, calculate the frequency of the RF signal
# the frequency of the RF signal is a list of RFsignaltimes and how often it occurs

for i in range(len(df['RFsignaltime'])):
    if i == 0: # first RF signal time
        RFsignaltimes = [df['RFsignaltime'][i]] # list of RF signal times
        RFsignalfrequency = [1] # list of RF signal frequencies
    else:
        if df['RFsignaltime'][i] == df['RFsignaltime'][i-1]: # if the RF signal time is the same as the previous RF signal time
            RFsignalfrequency[-1] += 1 # add 1 to the frequency of the last RF signal time
        # if the RF signal time is not the same as the previous RF signal time but matches a previous RF signal time, add 1 to the frequency of that RF signal time
        elif df['RFsignaltime'][i] in RFsignaltimes:
            RFsignalfrequency[RFsignaltimes.index(df['RFsignaltime'][i])] += 1
        else:
            RFsignaltimes.append(df['RFsignaltime'][i]) # add the RF signal time to the list of RF signal times
            RFsignalfrequency.append(1)

dfRF = pd.DataFrame({'RFsignaltimes': RFsignaltimes, 'RFsignalfrequency': RFsignalfrequency})
print(dfRF)

# plot a bar graph of the RF signal times and their frequencies
plt.bar(dfRF['RFsignaltimes'], dfRF['RFsignalfrequency'])
plt.xlabel('RF signal time (ns)')
plt.ylabel('Frequency')
plt.title('RF signal time vs. frequency')
plt.savefig(save_dir + 'RFsignaltimevsfrequency.png')
plt.show()

    


