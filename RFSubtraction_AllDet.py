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

save_dir = data_dir + '/' + 'RFSubtraction/AllDet/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print ('Creating directory: {}'.format(save_dir))

f = uproot.open(data_dir + ROOTFile)
tree0 = f['LaBr0Data']
EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det0 = "TimingSlow_LaBr0"
TimingFast_LaBr_Det0 = "TimingFast_LaBr0"
TimingFast_RF = "TimingFast_RF"
slowLaBr0calib1stFullRange = tree0[EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange].array()
slowtimeLaBr0 = tree0[TimingSlow_LaBr_Det0].array()
fasttimeLaBr0 = tree0[TimingFast_LaBr_Det0].array()
RFsignaltime0 = tree0[TimingFast_RF].array()
df0 = pd.DataFrame({'slowLaBr0calib1stFullRange': slowLaBr0calib1stFullRange, 'slowtimeLaBr0': slowtimeLaBr0, 'fasttimeLaBr0': fasttimeLaBr0, 'RFsignaltime0': RFsignaltime0})

tree1 = f['LaBr1Data']
EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det1 = "TimingSlow_LaBr1"
TimingFast_LaBr_Det1 = "TimingFast_LaBr1"
TimingFast_RF = "TimingFast_RF"
slowLaBr1calib1stFullRange = tree1[EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange].array()
slowtimeLaBr1 = tree1[TimingSlow_LaBr_Det1].array()
fasttimeLaBr1 = tree1[TimingFast_LaBr_Det1].array()
RFsignaltime1 = tree1[TimingFast_RF].array()
df1 = pd.DataFrame({'slowLaBr1calib1stFullRange': slowLaBr1calib1stFullRange, 'slowtimeLaBr1': slowtimeLaBr1, 'fasttimeLaBr1': fasttimeLaBr1, 'RFsignaltime1': RFsignaltime1})

tree2 = f['LaBr2Data']
EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det2 = "TimingSlow_LaBr2"
TimingFast_LaBr_Det2 = "TimingFast_LaBr2"
TimingFast_RF = "TimingFast_RF"
slowLaBr2calib1stFullRange = tree2[EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange].array()
slowtimeLaBr2 = tree2[TimingSlow_LaBr_Det2].array()
fasttimeLaBr2 = tree2[TimingFast_LaBr_Det2].array()
RFsignaltime2 = tree2[TimingFast_RF].array()
df2 = pd.DataFrame({'slowLaBr2calib1stFullRange': slowLaBr2calib1stFullRange, 'slowtimeLaBr2': slowtimeLaBr2, 'fasttimeLaBr2': fasttimeLaBr2, 'RFsignaltime2': RFsignaltime2})

tree3 = f['LaBr3Data']
EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange = "EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange"
TimingSlow_LaBr_Det3 = "TimingSlow_LaBr3"
TimingFast_LaBr_Det3 = "TimingFast_LaBr3"
TimingFast_RF = "TimingFast_RF"
slowLaBr3calib1stFullRange = tree3[EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange].array()
slowtimeLaBr3 = tree3[TimingSlow_LaBr_Det3].array()
fasttimeLaBr3 = tree3[TimingFast_LaBr_Det3].array()
RFsignaltime3 = tree3[TimingFast_RF].array()
df3 = pd.DataFrame({'slowLaBr3calib1stFullRange': slowLaBr3calib1stFullRange, 'slowtimeLaBr3': slowtimeLaBr3, 'fasttimeLaBr3': fasttimeLaBr3, 'RFsignaltime3': RFsignaltime3})

df0['slowtimeLaBr0'] = df0['slowtimeLaBr0']/10
df0['fasttimeLaBr0'] = df0['fasttimeLaBr0']/10
df0['RFsignaltime0'] = df0['RFsignaltime0']/10
df0 = df0.sort_values(by=['fasttimeLaBr0'])
df0 = df0.reset_index(drop=True)
df1['slowtimeLaBr1'] = df1['slowtimeLaBr1']/10
df1['fasttimeLaBr1'] = df1['fasttimeLaBr1']/10
df1['RFsignaltime1'] = df1['RFsignaltime1']/10
df1 = df1.sort_values(by=['fasttimeLaBr1'])
df1 = df1.reset_index(drop=True)
df2['slowtimeLaBr2'] = df2['slowtimeLaBr2']/10
df2['fasttimeLaBr2'] = df2['fasttimeLaBr2']/10
df2['RFsignaltime2'] = df2['RFsignaltime2']/10
df2 = df2.sort_values(by=['fasttimeLaBr2'])
df2 = df2.reset_index(drop=True)
df3['slowtimeLaBr3'] = df3['slowtimeLaBr3']/10
df3['fasttimeLaBr3'] = df3['fasttimeLaBr3']/10
df3['RFsignaltime3'] = df3['RFsignaltime3']/10
df3 = df3.sort_values(by=['fasttimeLaBr3'])
df3 = df3.reset_index(drop=True)


print('L0 df0', df0)
print('L1 df1', df1)
print('L2 df2', df2)
print('L3 df3', df3)

#---------------------------------------------------------------


def RFSubtraction(df0, df1, df2, df3, save_dir, ROOTFile):

    '''        # plot fasttimeLaBr0 - fasttimeRF vs slowLaBr0calib1stFullRange.  Where the time difference is -300 - 300 ns, color the points red.
    # plot the histogram
    c11 =root.TCanvas("c11","c11", 800, 600)
    c11.SetGrid()
    h14 = root.TH2F("h14", "h14", 6000, 0, 6000, 1000, -300, 300) # 100 ns range in 1 ps bins
    for i in range(len(df0)):
        h14.Fill(df0['slowLaBr0calib1stFullRange'][i], df0['fasttimeLaBr0'][i] - df0['RFsignaltime0'][i])
    h14.SetLineColor(2)
    h14.SetLineWidth(2)
    h14.SetTitle("fasttimeLaBr0 - fasttimeRF vs slowLaBr0calib1stFullRange")
    h14.GetXaxis().SetTitle("slowLaBr0calib1stFullRange")
    h14.GetYaxis().SetTitle("fasttimeLaBr0 - fasttimeRF")
    h14.Draw("colz")
    c11.SaveAs(save_dir + ROOTFile + "_fasttimeLaBr0_fasttimeRF_vs_slowLaBr0calib1stFullRange.png")

    # plot fasttimeLaBr1 - fasttimeRF vs slowLaBr1calib1stFullRange.  Where the time difference is -300 - 300 ns, color the points red.
    # plot the histogram
    c12 =root.TCanvas("c12","c12", 800, 600)
    c12.SetGrid()
    h15 = root.TH2F("h15", "h15", 6000, 0, 6000, 1000, -300, 300) # 100 ns range in 1 ps bins
    for i in range(len(df1)):
        h15.Fill(df1['slowLaBr1calib1stFullRange'][i], df1['fasttimeLaBr1'][i] - df1['RFsignaltime1'][i])
    h15.SetLineColor(2)
    h15.SetLineWidth(2)
    h15.SetTitle("fasttimeLaBr1 - fasttimeRF vs slowLaBr1calib1stFullRange")
    h15.GetXaxis().SetTitle("slowLaBr1calib1stFullRange")
    h15.GetYaxis().SetTitle("fasttimeLaBr1 - fasttimeRF")
    h15.Draw("colz")
    c12.SaveAs(save_dir + ROOTFile + "_fasttimeLaBr1_fasttimeRF_vs_slowLaBr1calib1stFullRange.png")

    # plot fasttimeLaBr2 - fasttimeRF vs slowLaBr2calib1stFullRange.  Where the time difference is -300 - 300 ns, color the points red.
    # plot the histogram
    c13 =root.TCanvas("c13","c13", 800, 600)
    c13.SetGrid()
    h16 = root.TH2F("h16", "h16", 6000, 0, 6000, 1000, -300, 300) # 100 ns range in 1 ps bins
    for i in range(len(df2)):
        h16.Fill(df2['slowLaBr2calib1stFullRange'][i], df2['fasttimeLaBr2'][i] - df2['RFsignaltime2'][i])
    h16.SetLineColor(2)
    h16.SetLineWidth(2)
    h16.SetTitle("fasttimeLaBr2 - fasttimeRF vs slowLaBr2calib1stFullRange")
    h16.GetXaxis().SetTitle("slowLaBr2calib1stFullRange")
    h16.GetYaxis().SetTitle("fasttimeLaBr2 - fasttimeRF")
    h16.Draw("colz")
    c13.SaveAs(save_dir + ROOTFile + "_fasttimeLaBr2_fasttimeRF_vs_slowLaBr2calib1stFullRange.png")

    # plot fasttimeLaBr3 - fasttimeRF vs slowLaBr3calib1stFullRange.  Where the time difference is -300 - 300 ns, color the points red.
    # plot the histogram
    c14 =root.TCanvas("c14","c14", 800, 600)
    c14.SetGrid()
    h17 = root.TH2F("h17", "h17", 6000, 0, 6000, 1000, -300, 300) # 100 ns range in 1 ps bins
    for i in range(len(df3)):
        h17.Fill(df3['slowLaBr3calib1stFullRange'][i], df3['fasttimeLaBr3'][i] - df3['RFsignaltime3'][i])
    h17.SetLineColor(2)
    h17.SetLineWidth(2)
    h17.SetTitle("fasttimeLaBr3 - fasttimeRF vs slowLaBr3calib1stFullRange")
    h17.GetXaxis().SetTitle("slowLaBr3calib1stFullRange")
    h17.GetYaxis().SetTitle("fasttimeLaBr3 - fasttimeRF")
    h17.Draw("colz")
    c14.SaveAs(save_dir + ROOTFile + "_fasttimeLaBr3_fasttimeRF_vs_slowLaBr3calib1stFullRange.png")'''
    
    # find the time difference between the RF and LaBr signals and add it to the dataframe as a new column
    df0['RFSubtraction'] = df0['fasttimeLaBr0'] - df0['RFsignaltime0']
    df1['RFSubtraction'] = df1['fasttimeLaBr1'] - df1['RFsignaltime1']
    df2['RFSubtraction'] = df2['fasttimeLaBr2'] - df2['RFsignaltime2']
    df3['RFSubtraction'] = df3['fasttimeLaBr3'] - df3['RFsignaltime3']

    # when RFSubtraction is between -300 and 300 ns, plot the slowLaBr0calib1stFullRange etc. histograms alongside the total slowLaBr0calibFullRange etc. histograms. Then subtract the slowLaBr0calib1stFullRange etc. histograms from the total slowLaBr0calibFullRange etc. histograms and plot the difference histograms.
    # plot the histograms
    c15 =root.TCanvas("c15","c15", 800, 600)
    c15.SetGrid()
    h18 = root.TH1F("h18", "h18", 6000, 0, 6000)
    h19 = root.TH1F("h19", "h19", 6000, 0, 6000)
    h20 = root.TH1F("h20", "h20", 6000, 0, 6000)
    h21 = root.TH1F("h21", "h21", 6000, 0, 6000)
    h22 = root.TH1F("h22", "h22", 6000, 0, 6000)
    h23 = root.TH1F("h23", "h23", 6000, 0, 6000)
    h24 = root.TH1F("h24", "h24", 6000, 0, 6000)
    h25 = root.TH1F("h25", "h25", 6000, 0, 6000)
    h26 = root.TH1F("h26", "h26", 6000, 0, 6000)
    h27 = root.TH1F("h27", "h27", 6000, 0, 6000)
    h28 = root.TH1F("h28", "h28", 6000, 0, 6000)
    h29 = root.TH1F("h29", "h29", 6000, 0, 6000)
    
    for i in range(len(df0)):
        if (df0['RFSubtraction'][i] >= 0 and df0['RFSubtraction'][i] <= 100) | (df0['RFSubtraction'][i] >= -270 and df0['RFSubtraction'][i] <= -200):
            h18.Fill(df0['slowLaBr0calib1stFullRange'][i])
        h19.Fill(df0['slowLaBr0calib1stFullRange'][i])
    h18.SetLineColor(2)
    h18.SetLineWidth(2)
    h18.SetTitle("Energy Spectrum")
    h18.GetXaxis().SetTitle("Energy (keV)")
    h18.GetYaxis().SetTitle("Counts")
    h18.Draw()
    h19.SetLineColor(4)
    h19.SetLineWidth(2)
    h19.Draw("same")
    h20 = h19.Clone()
    h20.Add(h18, -1)
    h20.SetLineColor(6)
    h20.SetLineWidth(2)
    h20.Draw("same")
    l = root.TLegend(0.7, 0.7, 0.9, 0.9)
    l.AddEntry(h18, "In Sync 0 to 100 ns", "l")
    l.AddEntry(h19, "Full Spectrum", "l")
    l.AddEntry(h20, "Full Spectrum - In Sync 0 to 100 ns", "l")
    l.Draw()
    c15.SaveAs(save_dir + ROOTFile + "_slowLaBr0calibFullRange.root")
    '''
    for i in range(len(df1)):
        if df1['RFSubtraction'][i] >=0 and df1['RFSubtraction'][i] <= 100:
            h21.Fill(df1['slowLaBr1calib1stFullRange'][i])
        h22.Fill(df1['slowLaBr1calib1stFullRange'][i])
    h21.SetLineColor(2)
    h21.SetLineWidth(2)
    h21.SetTitle("Energy Spectrum")
    h21.GetXaxis().SetTitle("Energy (keV)")
    h21.GetYaxis().SetTitle("Counts")
    h21.Draw()
    h22.SetLineColor(4)
    h22.SetLineWidth(2)
    h22.Draw("same")
    h23 = h22.Clone()
    h23.Add(h21, -1)
    h23.SetLineColor(6)
    h23.SetLineWidth(2)
    h23.Draw("same")
    l = root.TLegend(0.7, 0.7, 0.9, 0.9)
    l.AddEntry(h21, "In Sync 0 to 100 ns", "l")
    l.AddEntry(h22, "Full Spectrum", "l")
    l.AddEntry(h23, "Full Spectrum - In Sync 0 to 100 ns", "l")
    l.Draw()
    c15.SaveAs(save_dir + ROOTFile + "_slowLaBr1calibFullRange.root")

    for i in range(len(df2)):
        if df2['RFSubtraction'][i] >=0 and df2['RFSubtraction'][i] <=100:
            h24.Fill(df2['slowLaBr2calib1stFullRange'][i])
        h25.Fill(df2['slowLaBr2calib1stFullRange'][i])
    h24.SetLineColor(2)
    h24.SetLineWidth(2)
    h24.SetTitle("Energy Spectrum")
    h24.GetXaxis().SetTitle("Energy (keV)")
    h24.GetYaxis().SetTitle("Counts")
    h24.Draw()
    h25.SetLineColor(4)
    h25.SetLineWidth(2)
    h25.Draw("same")
    h26 = h25.Clone()
    h26.Add(h24, -1)    
    h26.SetLineColor(6)
    h26.SetLineWidth(2)
    h26.Draw("same")
    l = root.TLegend(0.7, 0.7, 0.9, 0.9)
    l.AddEntry(h24, "In Sync 0 to 100 ns", "l")
    l.AddEntry(h25, "Full Spectrum", "l")
    l.AddEntry(h26, "Full Spectrum - In Sync 0 to 100 ns", "l")
    l.Draw()
    c15.SaveAs(save_dir + ROOTFile + "_slowLaBr2calibFullRange.root")

    for i in range(len(df3)):
        if df3['RFSubtraction'][i] >=0 and df3['RFSubtraction'][i] <=100:
            h27.Fill(df3['slowLaBr3calib1stFullRange'][i])
        h28.Fill(df3['slowLaBr3calib1stFullRange'][i])
    h27.SetLineColor(2)
    h27.SetLineWidth(2)
    h27.SetTitle("Energy Spectrum")
    h27.GetXaxis().SetTitle("Energy (keV)")
    h27.GetYaxis().SetTitle("Counts")
    h27.Draw()
    h28.SetLineColor(4)
    h28.SetLineWidth(2)
    h28.Draw("same")
    h29 = h28.Clone()
    h29.Add(h27, -1)
    h29.SetLineColor(6)
    h29.SetLineWidth(2)
    h29.Draw("same")
    l = root.TLegend(0.7, 0.7, 0.9, 0.9)
    l.AddEntry(h27, "In Sync 0 to 100 ns", "l")
    l.AddEntry(h28, "Full Spectrum", "l")
    l.AddEntry(h29, "Full Spectrum - In Sync 0 to 100 ns", "l")
    l.Draw()
    c15.SaveAs(save_dir + ROOTFile + "_slowLaBr3calibFullRange.root")
    '''
    #---------------------------------------------------------------
    first_valid_time_L0 = df0['fasttimeLaBr0'][df0['fasttimeLaBr0'] > 1e13].iloc[0]
    first_valid_index_L0 = df0.index[df0['fasttimeLaBr0'] == first_valid_time_L0].tolist()[0]
    first_valid_time_L1 = df1['fasttimeLaBr1'][df1['fasttimeLaBr1'] > 1e13].iloc[0]
    first_valid_index_L1 = df1.index[df1['fasttimeLaBr1'] == first_valid_time_L1].tolist()[0]
    first_valid_time_L2 = df2['fasttimeLaBr2'][df2['fasttimeLaBr2'] > 1e13].iloc[0]
    first_valid_index_L2 = df2.index[df2['fasttimeLaBr2'] == first_valid_time_L2].tolist()[0]
    first_valid_time_L3 = df3['fasttimeLaBr3'][df3['fasttimeLaBr3'] > 1e13].iloc[0]
    first_valid_index_L3 = df3.index[df3['fasttimeLaBr3'] == first_valid_time_L3].tolist()[0]

    timediffLaBr0 = df0['fasttimeLaBr0'].diff()
    timediffLaBr1 = df1['fasttimeLaBr1'].diff()
    timediffLaBr2 = df2['fasttimeLaBr2'].diff()
    timediffLaBr3 = df3['fasttimeLaBr3'].diff()
    timediffRF = df0['RFsignaltime0'].diff()


    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex=True, figsize=(30,20))
    ax1.plot(timediffLaBr0, 'k-', label='L0')
    ax2.plot(timediffLaBr1, 'k-', label='L1')
    ax3.plot(timediffLaBr2, 'k-', label='L2')
    ax4.plot(timediffLaBr3, 'k-', label='L3')
    ax5.plot(timediffRF, 'r-', label='RF')
    ax1.legend( loc='upper right', fontsize=40)
    ax2.legend( loc='upper right', fontsize=40)
    ax3.legend( loc='upper right', fontsize=40)
    ax4.legend( loc='upper right', fontsize=40)
    ax5.legend( loc='upper right', fontsize=40)
    ax1.set_xlim(0, 10)
    ax2.set_xlim(0, 10)
    ax3.set_xlim(0, 10)
    ax4.set_xlim(0, 10)
    ax5.set_xlim(0, 10)
    ax5.set_xlabel('Time Difference (ns)', fontsize=40)
    ax1.set_ylabel('Counts', fontsize=40)
    ax2.set_ylabel('Counts', fontsize=40)
    ax3.set_ylabel('Counts', fontsize=40)
    ax4.set_ylabel('Counts', fontsize=40)
    ax5.set_ylabel('Counts', fontsize=40)
    ax1.set_title('Time difference between consecutive fast time L0 events', fontsize=40)
    ax2.set_title('Time difference between consecutive fast time L1 events', fontsize=40)
    ax3.set_title('Time difference between consecutive fast time L2 events', fontsize=40)
    ax4.set_title('Time difference between consecutive fast time L3 events', fontsize=40)
    ax5.set_title('Time difference between consecutive fast time RF events', fontsize=40)
    ax1.tick_params(axis='both', which='major', labelsize=40)
    ax2.tick_params(axis='both', which='major', labelsize=40)
    ax3.tick_params(axis='both', which='major', labelsize=40)
    ax4.tick_params(axis='both', which='major', labelsize=40)
    ax5.tick_params(axis='both', which='major', labelsize=40)
    plt.savefig(save_dir + 'timediffLaBr0andRF.png')

    first_time = min(first_valid_time_L0, first_valid_time_L1, first_valid_time_L2, first_valid_time_L3)

    h1 = root.TH1F("h1", "h1", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - L0
    h2 = root.TH1F("h2", "h2", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - L1
    h3 = root.TH1F("h3", "h3", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - L2
    h4 = root.TH1F("h4", "h4", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - L3
    r1 = root.TH1F("r1", "r1", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - RF

    #h5 = root.TH1F("h5", "h5", 100000 , 0, 100) # 100 ns range in 1 ps bins  - L0 -L1
    #h6 = root.TH1F("h6", "h6", 100000 , 0, 100) # 100 ns range in 1 ps bins  - L0 -L2
    #h7 = root.TH1F("h7", "h7", 100000 , 0, 100) # 100 ns range in 1 ps bins  - L0 -L3
    #h8 = root.TH1F("h8", "h8", 100000 , 0, 100) # 100 ns range in 1 ps bins  - L0 - RF
    
    c =root.TCanvas("c","c", 800, 600)
    c.SetGrid()
    for i in range(len(df0)):
        h1.Fill(df0['fasttimeLaBr0'][i])
    for i in range(len(df1)):
        h2.Fill(df1['fasttimeLaBr1'][i])
    for i in range(len(df2)):
        h3.Fill(df2['fasttimeLaBr2'][i])
    for i in range(len(df3)):
        h4.Fill(df3['fasttimeLaBr3'][i])
    for i in range(len(df3)):
        r1.Fill(df3['RFsignaltime3'][i])
    h1.SetLineColor(2)
    h1.SetLineWidth(2)
    h1.SetTitle("L0")
    h1.GetXaxis().SetTitle("Fast Time (ns)")
    h1.GetYaxis().SetTitle("Counts")
    h1.Draw()
    h2.SetLineColor(3)
    h2.SetLineWidth(2)
    h2.SetTitle("L1")
    h2.Draw("same")
    h3.SetLineColor(4)
    h3.SetLineWidth(2)
    h3.SetTitle("L2")
    h3.Draw("same")
    h4.SetLineColor(6)
    h4.SetLineWidth(2)
    h4.SetTitle("L3")
    h4.Draw("same")
    r1.SetLineColor(1)
    r1.SetLineWidth(2)
    r1.SetTitle("RF")
    r1.Draw("same")
    c.BuildLegend()
    c.SaveAs(save_dir + ROOTFile + "_fast_L0_L1_L2_RF.png")
    c.SaveAs(save_dir + ROOTFile + "_fast_L0_L1_L2_RF.root")

    #---------------------------------------------------------------
    '''# draw h1 - h2 
    c1 =root.TCanvas("c1","c1", 800, 600)
    c1.SetGrid()
    h5 = h1.Clone()
    h5.Add(h2, -1)
    h5.SetLineColor(2)
    h5.SetLineWidth(2)
    h5.SetTitle("L0 - L1")
    h5.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h5.GetYaxis().SetTitle("Counts")
    h5.Draw()
    c1.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_L1.png")
    c1.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_L1.root")'''

    #---------------------------------------------------------------
    ''' # draw h1 - h3 (l0 - l2)
    c2 =root.TCanvas("c2","c2", 800, 600)
    c2.SetGrid()
    h6 = h1.Clone()
    h6.Add(h3, -1)
    h6.SetLineColor(2)
    h6.SetLineWidth(2)
    h6.SetTitle("L0 - L2")
    h6.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h6.GetYaxis().SetTitle("Counts")
    h6.Draw()
    c2.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_L2.png")
    c2.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_L2.root")
    '''

    #---------------------------------------------------------------
    # draw h1 - h4 (l0 - l3)
    '''c3 =root.TCanvas("c3","c3", 800, 600)
    c3.SetGrid()
    h7 = h1.Clone()
    h7.Add(h4, -1)
    h7.SetLineColor(2)
    h7.SetLineWidth(2)
    h7.SetTitle("L0 - L3")
    h7.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h7.GetYaxis().SetTitle("Counts")
    h7.Draw()
    c3.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_L3.png")
    c3.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_L3.root")'''

    #---------------------------------------------------------------
    # draw h1 - r1 (l0 - rf)
    '''c4 =root.TCanvas("c4","c4", 800, 600)
    c4.SetGrid()
    h8 = h1.Clone()
    h8.Add(r1, -1)
    h8.SetLineColor(2)
    h8.SetLineWidth(2)
    h8.SetTitle("L0 - RF")
    h8.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h8.GetYaxis().SetTitle("Counts")
    h8.Draw()
    c4.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_RF.png")
    c4.SaveAs(save_dir + ROOTFile + "_fast_L0_minus_RF.root")'''

    #---------------------------------------------------------------
    # the time difference between consecutive events in slowtimeLaBr0
    # plot the histogram
    c5 =root.TCanvas("c5","c5", 800, 600)
    c5.SetGrid()
    h9 = root.TH1F("h9", "h9", 100000 , 0, 100) # 100 ns range in 1 ps bins
    for i in range(len(df0)-1):
        h9.Fill(df0['fasttimeLaBr0'][i+1] - df0['fasttimeLaBr0'][i])
    h9.SetLineColor(2)
    h9.SetLineWidth(2)
    h9.SetTitle("Time difference between consecutive events in L0")
    h9.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h9.GetYaxis().SetTitle("Counts")
    h9.Draw()
    c5.SaveAs(save_dir + ROOTFile + "_fast_time_diff_L0.png")

    #---------------------------------------------------------------
    # the time difference between consecutive events in fasttimeLaBr1
    # plot the histogram
    c6 =root.TCanvas("c6","c6", 800, 600)
    c6.SetGrid()
    h10 = root.TH1F("h10", "h10", 100000 , 0, 100) # 100 ns range in 1 ps bins
    for i in range(len(df1)-1):
        h10.Fill(df1['fasttimeLaBr1'][i+1] - df1['fasttimeLaBr1'][i])
    h10.SetLineColor(2)
    h10.SetLineWidth(2)
    h10.SetTitle("Time difference between consecutive events in L1")
    h10.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h10.GetYaxis().SetTitle("Counts")
    h10.Draw()
    c6.SaveAs(save_dir + ROOTFile + "_fast_time_diff_L1.png")

    #---------------------------------------------------------------
    # the time difference between consecutive events in fasttimeLaBr2
    # plot the histogram
    c7 =root.TCanvas("c7","c7", 800, 600)
    c7.SetGrid()
    h11 = root.TH1F("h11", "h11", 100000 , 0, 100) # 100 ns range in 1 ps bins
    for i in range(len(df3)-1):
        h11.Fill(df2['fasttimeLaBr2'][i+1] - df2['fasttimeLaBr2'][i])
    h11.SetLineColor(2)
    h11.SetLineWidth(2)
    h11.SetTitle("Time difference between consecutive events in L2")
    h11.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h11.GetYaxis().SetTitle("Counts")
    h11.Draw()
    c7.SaveAs(save_dir + ROOTFile + "_fast_time_diff_L2.png")

    #---------------------------------------------------------------
    # the time difference between consecutive events in fasttimeLaBr3
    # plot the histogram
    c8 =root.TCanvas("c8","c8", 800, 600)
    c8.SetGrid()
    h12 = root.TH1F("h12", "h12", 100000 , 0, 100) # 100 ns range in 1 ps bins
    for i in range(len(df3)-1):
        h12.Fill(df3['fasttimeLaBr3'][i+1] - df3['fasttimeLaBr3'][i])
    h12.SetLineColor(2)
    h12.SetLineWidth(2)
    h12.SetTitle("Time difference between consecutive events in L3")
    h12.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h12.GetYaxis().SetTitle("Counts")
    h12.Draw()
    c8.SaveAs(save_dir + ROOTFile + "_fast_time_diff_L3.png")

    #---------------------------------------------------------------
    # the time difference between consecutive events in fasttimeRF
    # plot the histogram
    c9 =root.TCanvas("c9","c9", 800, 600)
    c9.SetGrid()
    h13 = root.TH1F("h13", "h13", 100000 , 0, 100) # 100 ns range in 1 ps bins
    for i in range(len(df3)-1):
        h13.Fill(df3['RFsignaltime3'][i+1] - df3['RFsignaltime3'][i])
    h13.SetLineColor(2)
    h13.SetLineWidth(2)
    h13.SetTitle("Time difference between consecutive events in RF")
    h13.GetXaxis().SetTitle("Fast Time Difference (ns)")
    h13.GetYaxis().SetTitle("Counts")
    h13.Draw()
    c9.SaveAs(save_dir + ROOTFile + "_fast_time_diff_RF.png")

    # Draw h9to h13 on the same canvas by dividing the canvas into 5 pads
    c10 = root.TCanvas("c10", "c10", 800, 600)
    c10.Divide(2,3)
    c10.cd(1)
    h9.Draw()
    c10.cd(2)
    h10.Draw()
    c10.cd(3)
    h11.Draw()
    c10.cd(4)
    h12.Draw()
    c10.cd(5)   
    h13.Draw()
    c10.SaveAs(save_dir + ROOTFile + "_fast_time_diff_all.png")
    c10.SaveAs(save_dir + ROOTFile + "_fast_time_diff_all.root")


    

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN CODE
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():
    RFSubtraction(df0, df1, df2, df3, save_dir, ROOTFile)

#---------------------------------------------------------------
if __name__ == "__main__":
    main()
