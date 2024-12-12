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
tree0 = f['LaBrData']
EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange = "slowECalibL0"
TimingSlow_LaBr_Det0 = "timeSL0"
TimingFast_LaBr_Det0 = "timeFL0"
TimingFast_RF = "timeRF"
slowLaBr0calib1stFullRange = tree0[EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange].array()
slowtimeLaBr0 = tree0[TimingSlow_LaBr_Det0].array()
fasttimeLaBr0 = tree0[TimingFast_LaBr_Det0].array()
RFsignaltime0 = tree0[TimingFast_RF].array()
df0 = pd.DataFrame({'slowLaBr0calib1stFullRange': slowLaBr0calib1stFullRange, 'slowtimeLaBr0': slowtimeLaBr0, 'fasttimeLaBr0': fasttimeLaBr0, 'RFsignaltime0': RFsignaltime0})


df0['slowtimeLaBr0'] = df0['slowtimeLaBr0']/10
df0['fasttimeLaBr0'] = df0['fasttimeLaBr0']/10
df0['RFsignaltime0'] = df0['RFsignaltime0']/10
df0 = df0.sort_values(by=['fasttimeLaBr0'])
df0 = df0.reset_index(drop=True)

print('L0 df0', df0)

#---------------------------------------------------------------


def RFSubtraction(df0, save_dir, ROOTFile):

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
    # when RFSubtraction is between -300 and 300 ns, plot the slowLaBr0calib1stFullRange etc. histograms alongside the total slowLaBr0calibFullRange etc. histograms. Then subtract the slowLaBr0calib1stFullRange etc. histograms from the total slowLaBr0calibFullRange etc. histograms and plot the difference histograms.
    # plot the histograms
    c15 =root.TCanvas("c15","c15", 800, 600)
    c15.SetGrid()
    h18 = root.TH1F("h18", "h18", 8000, 0, 8000)
    h19 = root.TH1F("h19", "h19", 8000, 0, 8000)
    h20 = root.TH1F("h20", "h20", 8000, 0, 8000)
    
    for i in range(len(df0)):
        if (df0['RFSubtraction'][i] >= 0 and df0['RFSubtraction'][i] <= 100) | (df0['RFSubtraction'][i] >= -270 and df0['RFSubtraction'][i] <= -200):
            h18.Fill(df0['slowLaBr0calib1stFullRange'][i])
        h19.Fill(df0['slowLaBr0calib1stFullRange'][i])
    
    for i in range(8000):
        h18.SetBinContent(i, h18.GetBinContent(i)*2)
    
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
    l.AddEntry(h18, "Out-of-sync", "l")
    l.AddEntry(h19, "In-sync", "l")
    l.AddEntry(h20, "Prompt", "l")
    l.Draw()
    c15.SaveAs(save_dir + ROOTFile + "_slowLaBr0calibFullRange.root")
    
    #---------------------------------------------------------------
    """ first_valid_time_L0 = df0['fasttimeLaBr0'][df0['fasttimeLaBr0'] > 1e13].iloc[0]
    first_valid_index_L0 = df0.index[df0['fasttimeLaBr0'] == first_valid_time_L0].tolist()[0]

    timediffLaBr0 = df0['fasttimeLaBr0'].diff()
    timediffRF = df0['RFsignaltime0'].diff()


    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex=True, figsize=(30,20))
    ax1.plot(timediffLaBr0, 'k-', label='L0')
    ax5.plot(timediffRF, 'r-', label='RF')
    ax1.legend( loc='upper right', fontsize=40)
    ax5.legend( loc='upper right', fontsize=40)
    ax1.set_xlim(0, 10)
    ax5.set_xlim(0, 10)
    ax5.set_xlabel('Time Difference (ns)', fontsize=40)
    ax1.set_ylabel('Counts', fontsize=40)
    ax5.set_ylabel('Counts', fontsize=40)
    ax1.set_title('Time difference between consecutive fast time L0 events', fontsize=40)
    ax5.set_title('Time difference between consecutive fast time RF events', fontsize=40)
    ax1.tick_params(axis='both', which='major', labelsize=40)
    ax5.tick_params(axis='both', which='major', labelsize=40)
    plt.savefig(save_dir + 'timediffLaBr0andRF.png')

    first_time = min(first_valid_time_L0)

    h1 = root.TH1F("h1", "h1", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - L0
    r1 = root.TH1F("r1", "r1", 900000 , first_time, first_time + 900) # 300 ns range in 100 ps bins - RF
    
    c =root.TCanvas("c","c", 800, 600)
    c.SetGrid()
    for i in range(len(df0)):
        h1.Fill(df0['fasttimeLaBr0'][i])
    h1.SetLineColor(2)
    h1.SetLineWidth(2)
    h1.SetTitle("L0")
    h1.GetXaxis().SetTitle("Fast Time (ns)")
    h1.GetYaxis().SetTitle("Counts")
    h1.Draw()
    r1.SetLineColor(1)
    r1.SetLineWidth(2)
    r1.SetTitle("RF")
    r1.Draw("same")
    c.BuildLegend()
    c.SaveAs(save_dir + ROOTFile + "_fast_L0_L1_L2_RF.png")
    c.SaveAs(save_dir + ROOTFile + "_fast_L0_L1_L2_RF.root")
 """


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
    c10.cd(5)   
    h13.Draw()
    c10.SaveAs(save_dir + ROOTFile + "_fast_time_diff_all.png")
    c10.SaveAs(save_dir + ROOTFile + "_fast_time_diff_all.root")


    

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN CODE
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():
    RFSubtraction(df0, save_dir, ROOTFile)

#---------------------------------------------------------------
if __name__ == "__main__":
    main()
