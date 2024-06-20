
#------------------------------------------------------------------
# python import statements
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
import decimal
import time as time
from ydata_profiling import ProfileReport

# run this script in multiple threads
os.environ["NUMEXPR_MAX_THREADS"] = "16"

#------------------------------------------------------------------
print('________________________________________', '\n')
data_dir = sys.argv[1] # e.g. /home/physicist/Desktop/Research/Research_Data/2021_08_05_LaBr3_POLARIS/2021_08_05_LaBr3_POLARIS_1/
MyFile = uproot.open(data_dir)
print('MyFile:\n', MyFile)
print('________________________________________', '\n')
MyTree = MyFile['LaBrData']
insyncE_L0 = "insyncEL0"   
slowECalib_L0 = "slowECalibL0" 
fastTime_L0 = "fastTimeL0"
slowTime_L0 = "slowTimeL0"
insyncE_L1 = "insyncEL1"
slowECalib_L1 = "slowECalibL1"
fastTime_L1 = "fastTimeL1"
slowTime_L1 = "slowTimeL1"
insyncE_L2 = "insyncEL2"
slowECalib_L2 = "slowECalibL2"
fastTime_L2 = "fastTimeL2"
slowTime_L2 = "slowTimeL2"
insyncE_L3 = "insyncEL3"
slowECalib_L3 = "slowECalibL3"
fastTime_L3 = "fastTimeL3"
slowTime_L3 = "slowTimeL3"
slowE_POLARIS = "slowEPOLARIS"
slowTime_POLARIS = "slowTimePOLARIS"

#  form an array for each leaf
insyncE_L0 = MyTree[insyncE_L0].array()
slowECalib_L0 = MyTree[slowECalib_L0].array()
fastTime_L0 = MyTree[fastTime_L0].array()
slowTime_L0 = MyTree[slowTime_L0].array()
insyncE_L1 = MyTree[insyncE_L1].array()
slowECalib_L1 = MyTree[slowECalib_L1].array()
fastTime_L1 = MyTree[fastTime_L1].array()
slowTime_L1 = MyTree[slowTime_L1].array()
insyncE_L2 = MyTree[insyncE_L2].array()
slowECalib_L2 = MyTree[slowECalib_L2].array()
fastTime_L2 = MyTree[fastTime_L2].array()
slowTime_L2 = MyTree[slowTime_L2].array()
insyncE_L3 = MyTree[insyncE_L3].array()
slowECalib_L3 = MyTree[slowECalib_L3].array()
fastTime_L3 = MyTree[fastTime_L3].array()
slowTime_L3 = MyTree[slowTime_L3].array()
slowE_POLARIS = MyTree[slowE_POLARIS].array()
slowTime_POLARIS = MyTree[slowTime_POLARIS].array()

# convert the arrays into a dataframe
# dataLaBr = {'insyncE_L0': insyncE_L0, 'slowECalib_L0': slowECalib_L0, 'fastTime_L0': fastTime_L0, 'slowTime_L0': slowTime_L0, 'insyncE_L1': insyncE_L1, 'slowECalib_L1': slowECalib_L1, 'fastTime_L1': fastTime_L1, 'slowTime_L1': slowTime_L1, 'insyncE_L2': insyncE_L2, 'slowECalib_L2': slowECalib_L2, 'fastTime_L2': fastTime_L2, 'slowTime_L2': slowTime_L2, 'insyncE_L3': insyncE_L3, 'slowECalib_L3': slowECalib_L3, 'fastTime_L3': fastTime_L3, 'slowTime_L3': slowTime_L3, 'slowE_POLARIS': slowE_POLARIS, 'slowTime_POLARIS': slowTime_POLARIS}
dataLaBr = {'slowECalib_L0': slowECalib_L0, 'fastTime_L0': fastTime_L0, 'slowTime_L0': slowTime_L0, 'slowECalib_L1': slowECalib_L1, 'fastTime_L1': fastTime_L1, 'slowTime_L1': slowTime_L1, 'slowECalib_L2': slowECalib_L2, 'fastTime_L2': fastTime_L2, 'slowTime_L2': slowTime_L2, 'slowECalib_L3': slowECalib_L3, 'fastTime_L3': fastTime_L3, 'slowTime_L3': slowTime_L3, 'slowE_POLARIS': slowE_POLARIS, 'slowTime_POLARIS': slowTime_POLARIS}
dataLaBr_df = pd.DataFrame(dataLaBr)
dataLaBr_df.reset_index(drop = True, inplace = True)
print('dataLaBr_df:\n',dataLaBr_df)

# dir = str(data_dir[:-8])
# data_small = dataLaBr_df.head(1000)
# data_small.to_csv(dir+'LaBr_Data.txt', sep = '\t', index = False)

# profile = ProfileReport(dataLaBr_df, title="LaBr Data Profile Report")

# profile.to_file("LaBr_Data_Profile_Report.html")

# import subprocess

# html_file = 'your_report.html'
# # save the html file to the save_dir
# subprocess.call(['cp', html_file, dir])

print('________________________________________', '\n')