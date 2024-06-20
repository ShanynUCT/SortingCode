################################################################################
#                               to run:                                        #
# python3 PlotTwoCalibrations.py /directory/to/rootfile(R**_rawData.root)/  /directory/to/rootfile(R**_CalibratedData.root)/             #
__author__ = "Shanyn Hart"
__date__ = "2022-08-22"
__version__ = "1.0"
                                            #
################################################################################
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
#os.environ["OPENBLAS_NUM_THREADS"] = "200"
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()}, language_level=3)

data_dir = sys.argv[1]
EuCalibratedFileName = sys.argv[2] 
EuCalibratedFile = uproot.open(data_dir + EuCalibratedFileName)

save_dir = data_dir + '/' + 'timing/' 
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))


#************************ Eu Calibrated Data (_rawData.root) ************************
EuCaliBratedTreeL0 = EuCalibratedFile['LaBr0Data']
TimingSlow_LaBr_Det0 = "TimingSlow_LaBr_Det0"
Time_LaBrL0 = EuCaliBratedTreeL0[TimingSlow_LaBr_Det0].array()
L0 = {'Time_LaBrL0': Time_LaBrL0}
L0 = pd.DataFrame(L0)
# time = time * 10e-9 to be in seconds
L0['Time_LaBrL0'] = L0['Time_LaBrL0'] * 10e-9
L0.reset_index(drop = True, inplace = True)

EuCaliBratedTreeL1 = EuCalibratedFile['LaBr1Data']
TimingSlow_LaBr1 = "TimingSlow_LaBr1"
Time_LaBrL1 = EuCaliBratedTreeL1[TimingSlow_LaBr1].array()
L1 = {'Time_LaBrL1': Time_LaBrL1}
L1 = pd.DataFrame(L1)
# time = time * 10e-9 to be in seconds
L1['Time_LaBrL1'] = L1['Time_LaBrL1'] * 10e-9
L1.reset_index(drop = True, inplace = True)

EuCaliBratedTreeL2 = EuCalibratedFile['LaBr2Data']
TimingSlow_LaBr2 = "TimingSlow_LaBr2"
Time_LaBrL2 = EuCaliBratedTreeL2[TimingSlow_LaBr2].array()
L2 = {'Time_LaBrL2': Time_LaBrL2}
L2 = pd.DataFrame(L2)
# time = time * 10e-9 to be in seconds
L2['Time_LaBrL2'] = L2['Time_LaBrL2'] * 10e-9
L2.reset_index(drop = True, inplace = True)

EuCaliBratedTreeL3 = EuCalibratedFile['LaBr3Data']
TimingSlow_LaBr_Det3 = "TimingSlow_LaBr_Det3"
Time_LaBrL3 = EuCaliBratedTreeL3[TimingSlow_LaBr_Det3].array()
L3 = {'Time_LaBrL3': Time_LaBrL3}
L3 = pd.DataFrame(L3)
# time = time * 10e-9 to be in seconds
L3['Time_LaBrL3'] = L3['Time_LaBrL3'] * 10e-9
L3.reset_index(drop = True, inplace = True)

max = 47753457.26562487
bins = 47753457


#*****************************************************************************************
c1 = root.TCanvas("c1", "c1", 2000, 1500)
# 10 ns bins
h1 = root.TH1F("h1", "L0", bins, 0, max)
h2 = root.TH1F("h2", "L1", bins, 0, max)
h3 = root.TH1F("h3", "L2", bins, 0, max)
h4 = root.TH1F("h4", "L3", bins, 0, max)

for i in range(len(L0['Time_LaBrL0'])):
    h1.Fill(L0['Time_LaBrL0'][i])
for i in range(len(L1['Time_LaBrL1'])):
    h2.Fill(L1['Time_LaBrL1'][i])
for i in range(len(L2['Time_LaBrL2'])):
    h3.Fill(L2['Time_LaBrL2'][i])
for i in range(len(L3['Time_LaBrL3'])):
    h4.Fill(L3['Time_LaBrL3'][i])
h1.SetLineColor(2)
h2.SetLineColor(4)
h3.SetLineColor(6)
h4.SetLineColor(8)
h1.SetLineWidth(2)
h2.SetLineWidth(2)
h3.SetLineWidth(2)
h4.SetLineWidth(2)
h1.SetStats(0)
h2.SetStats(0)
h3.SetStats(0)
h4.SetStats(0)
h1.GetXaxis().SetTitle("Time (s)")
h1.GetYaxis().SetTitle("Counts")
h1.GetXaxis().SetTitleSize(0.05)
h1.GetYaxis().SetTitleSize(0.05)
h1.GetXaxis().SetTitleOffset(0.9)
h1.GetYaxis().SetTitleOffset(0.9)
h1.GetXaxis().SetLabelSize(0.05)
h1.GetYaxis().SetLabelSize(0.05)
h1.GetXaxis().SetRangeUser(0, max)
h1.Draw()
h2.Draw("same")
h3.Draw("same")
h4.Draw("same")
leg = root.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(h1, "L0", "l")
leg.AddEntry(h2, "L1", "l")
leg.AddEntry(h3, "L2", "l")
leg.AddEntry(h4, "L3", "l")
leg.Draw()
c1.SaveAs(save_dir + '/' + 'TimeNoBeam.root')
c1.Close()

