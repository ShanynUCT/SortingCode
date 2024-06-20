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




data_dir = sys.argv[1]
data_filename = 'mod51.txt'


save_dir = sys.argv[1] + 'TimeSyncResults/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))


mod_event_data = pd.read_csv(data_dir + data_filename, sep = '	', header = None)  
mod_event_data.columns = ['scatters', 'x', 'y', 'z', 'energy', 'time']   

mod_event_data = mod_event_data[ (mod_event_data['scatters'] == 1)]
mod_event_data_df = pd.DataFrame(mod_event_data)
mod_event_data_df['energy'] = mod_event_data_df['energy']/1e3

LaBrFileName = sys.argv[2] # just the name of the file, e.g. R00_rawData.root 
MyFile = uproot.open(data_dir + LaBrFileName)

MyTree3 = MyFile['LaBr0Data']
EnergySlow_LaBr3_1stOrderLowRange = "EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange"
EnergySlow_LaBr3_1stOrderLowRange = MyTree3[EnergySlow_LaBr3_1stOrderLowRange].array()
# convert the arrays into a dataframe
LaBrData_Det3 = {'EnergySlow_LaBr3_1stOrderLowRange': EnergySlow_LaBr3_1stOrderLowRange}
LaBrData_Det3_df = pd.DataFrame(LaBrData_Det3)
LaBrData_Det3_df.reset_index(drop = True, inplace = True)

c1 = root.TCanvas("c1","c1", 1000, 800)
c1.SetLogy()
h1 = root.TH1F("h1", "Energy spectra of LaBr_3:Ce detector vs. POLARIS detector for a Water target", 7000, 0, 7000)
h2 = root.TH1F("h2", "Energy spectra of LaBr_3:Ce detector vs. POLARIS detector for a Water target", 7000, 0, 7000)
for i in range(len(LaBrData_Det3_df['EnergySlow_LaBr3_1stOrderLowRange'])):
    h1.Fill(LaBrData_Det3_df['EnergySlow_LaBr3_1stOrderLowRange'].iloc[i])
for i in range(len(mod_event_data_df['energy'])):
    h2.Fill(mod_event_data_df['energy'].iloc[i])
h1.SetLineColor(2)
h1.SetLineWidth(2)
h2.SetLineColor(4)
h2.SetLineWidth(2)
h1.GetXaxis().SetTitle("Energy (keV/bin)")
h1.GetYaxis().SetTitle("Counts")
h1.Draw()
h2.Draw("same")
l = root.TLegend(0.7,0.7,0.9,0.9)
l.AddEntry(h1,"LaBr_3:Ce","l")
l.AddEntry(h2,"POLARIS","l")
l.Draw()
c1.SaveAs(save_dir + 'POLARIS_LaBr3_D0_Energy.root')