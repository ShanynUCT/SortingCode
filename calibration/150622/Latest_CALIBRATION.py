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
import math as math
import csv as csv
import array as array
import io as io
import sys, os
from math import sin, cos, pi, log, floor
import cProfile
import re
import pylab
import pyroot as pr
import decimal
import time as time

f = root.TFile.Open(sys.argv[1], "READ")
g = root.TFile.Open(sys.argv[2], "READ")
saveDirectory = sys.argv[1].split("/")[0:-1]
# if calibration_plots doesnt exist in saveDirectory, create it
plotsDirectory = "/".join(saveDirectory)+"/calibration_plots"
if not os.path.exists(plotsDirectory):
    os.makedirs(plotsDirectory)
    os.makedirs(plotsDirectory+"/gauss_fits")
    os.makedirs(plotsDirectory+"/energy_resolution")
    os.makedirs(plotsDirectory+"/energy_calibration")

# the first 11 channels correspond to the 11 peaks of 152Eu source, the last 3 peaks correspond to AmBe source
#         | 152 Eu                                                              | AmBe             |
energys = [121.8, 244.7, 344.3, 411.1, 444, 511, 778.9, 867.4, 964.1, 1408, 1470, 3416, 3927, 4438]

channels_L0 = [161, 325, 457, 547, 592, 682, 1037, 1152, 1282, 1872, 1940, 4494, 5170, 5820]
channels_L1 = [150, 304.6, 428, 505, 547, 629, 959, 1065, 1187, 1728, 1803, 4201, 4793, 5403]
channels_L2 = [237, 472, 664, 795, 863, 987, 1505, 1670, 1858, 2716, 2822, 6530, 7507, nan]
channels_L3 = [169, 339, 480, 567, 621, 712, 1088, 1206, 1345, 1960, 2042, 4742, 5441, 6115]


def getHistos(f,g):
    hists_152Eu = []
    hists_AmBe = []
    for i in range(0,4):
        hists_152Eu.append(f.Get("L"+str(i)+"_Uncalibrated"))
        hists_AmBe.append(g.Get("L"+str(i)+"_Uncalibrated"))
    return hists_152Eu, hists_AmBe


def getFitParams(hists_152Eu, hists_AmBe):
    fitParams = []
    for i in range(0,4):
        fitParams.append([])
        channels = eval("channels_L"+str(i))
        for j in range(0,10):
            #create a canvas to draw the histogram
            c = root.TCanvas("c", "c", 1000,800)
            fitParams[i].append([])
            fitParams[i][j].append(channels[j])
            flexibility = 30
            g1 = root.TF1("g1","gaus",channels[j]-flexibility,channels[j]+flexibility) # gauss peak
            # look on either side of the peak for the background
            pol1 = root.TF1("pol1","pol1",channels[j]-flexibility-60,channels[j]+flexibility+60) # linear background
            total = root.TF1( "total", "gaus(0) + pol1(3)", channels[j]-flexibility,channels[j]+flexibility) # total func
            h1 = hists_152Eu[i].Clone()
            h1.SetLineColor(root.kBlack)
            g1.SetLineColor(root.kBlue)
            pol1.SetLineColor(root.kGreen)
            total.SetLineColor(root.kRed)
            # Fit the histogram with gaussian
            h1.Fit(g1, "R")
            h1.Fit(pol1, "R+")
            h1.SetStats(000000000)
            
            
            p0 = g1.GetParameter(0) # amplitude of gauss peak
            p1 = g1.GetParameter(1) # mean value of gaussian
            # p2 is sigma of gaussian
            p2 = g1.GetParameter(2)
            print("amplitude value from gauss fit= ", p0 ,"\nmean value from gauss fit= ", p1, "\nsigma value from gauss fit= ", p2, "\nFWHM= ", 2.355*p2)
            
            p3 = pol1.GetParameter(0) # background value
            p4 = pol1.GetParameter(1) # slope of background
            
            e0 = g1.GetParError(0) # error on amplitude
            e1 = g1.GetParError(1) # error on mean
            e2 = g1.GetParError(2) # error on sigma

            e3 = pol1.GetParError(0) # error on background
            e4 = pol1.GetParError(1) # error on slope of background

            par = [p0, p1, p2, p3, p4] 
            h1.Fit(total, "R+")

            gaussParameters = [p0, p1, p2]
            gaussParameters = array.array('d', gaussParameters)
            MyConvarienceMatrix = [e0,0,0,0,e1,0,0,0,e2]
            MyConvarienceMatrix = array.array('d', MyConvarienceMatrix)


            HistoBinning = 1
            MyIntegral = g1.Integral( p1-3*p2, p1+3*p2 ) / HistoBinning #p4 is mean value and p5 is sigma value. Gauss range is about 8sigma
            IntegralError = g1.IntegralError(p1-3*p2, p1+3*p2, gaussParameters, MyConvarienceMatrix) / HistoBinning
            
            print("Fitting g1 integral and error: ", MyIntegral, "  ", IntegralError)
            print("Fitting pol1 integral and error: ", pol1.Integral(channels[j]-flexibility-50,channels[j]+flexibility+50))
            print("Fitting total integral and error: ", total.Integral( p1-3*p2, p1+3*p2 ))
            
            h1.GetXaxis().SetRangeUser(channels[j]-500,channels[j]+500)

            g1.Draw("same")
            pol1.Draw("same")
            total.Draw("same")

            xm1=0.6
            ym1=0.6
            xm2=0.8
            ym2=0.88
            legend1 = root.TLegend(xm1,ym1,xm2,ym2)
            legend1.SetHeader("^{152}Eu " + str(energys[j]) + " keV")
            legend1.AddEntry(g1,"Gauss Peak","l")
            legend1.AddEntry(pol1,"Polynomial Background","l")
            legend1.AddEntry(total,"Total Fit","l")
            legend1.SetBorderSize(0)
            legend1.SetTextSize(0.03)
            legend1.Draw()
            c.Update()
            c.SaveAs(plotsDirectory+"/gauss_fits/L" + str(i) + "_152Eu_" + str(energys[j]) + ".png")
            c.Close()
            fitParams[i][j].append(p1) #
            fitParams[i][j].append(p2)
            print("fitParams: ", fitParams)


        for j in range(10,13):
            c = root.TCanvas("c", "c", 1000,800)
            fitParams[i].append([])
            fitParams[i][j].append(channels[j])
            flexibility = 150
            g1 = root.TF1("g1","gaus",channels[j]-flexibility,channels[j]+flexibility) # gauss peak
            pol1 = root.TF1("pol1", "pol1", channels[j]-flexibility-300,channels[j]+flexibility+300) # pol backgrnd
            total = root.TF1( "total", "gaus(0) + pol1(3)", channels[j]-flexibility,channels[j]+flexibility) # total func
            h1 = hists_AmBe[i].Clone()
            h1.SetLineColor(root.kBlack)
            g1.SetLineColor(root.kBlue)
            pol1.SetLineColor(root.kGreen)
            total.SetLineColor(root.kRed)
            h1.Fit(g1, "R")
            h1.Fit(pol1, "R+")
            h1.SetStats(000000000)
            
            par = [p0, p1, p2, p3, p4]
            h1.Fit(total, "R+")

            p0 = g1.GetParameter(0) # amplitude of gauss peak
            p1 = g1.GetParameter(1) # mean value of gaussian
            # p2 is sigma of gaussian, width of gaussian=8sigma
            p2 = g1.GetParameter(2)
            print("amplitude value from gauss fit= ", p0 ,"\nmean value from gauss fit= ", p1, "\nsigma value from gauss fit= ", p2, "\nFWHM= ", 2.355*p2)
            
            p3 = pol1.GetParameter(0) # background value
            p4 = pol1.GetParameter(1) # slope of background
            
            e0 = g1.GetParError(0) # error on amplitude
            e1 = g1.GetParError(1) # error on mean
            e2 = g1.GetParError(2) # error on sigma

            e3 = pol1.GetParError(0) # error on background
            e4 = pol1.GetParError(1) # error on slope of background


            gaussParameters = [p0, p1, p2]
            gaussParameters = array.array('d', gaussParameters)
            MyConvarienceMatrix = [e0,0,0,0,e1,0,0,0,e2]
            MyConvarienceMatrix = array.array('d', MyConvarienceMatrix)

            h1.GetXaxis().SetRangeUser(channels[j]-500,channels[j]+500)

            HistoBinning = 1
            MyIntegral = g1.Integral( p1-3*p2, p1+3*p2 ) / HistoBinning #p4 is mean value and p5 is sigma value. Gauss range is about 8sigma
            IntegralError = g1.IntegralError(p1-3*p2, p1+3*p2, gaussParameters, MyConvarienceMatrix) / HistoBinning
            print("Fitting g1 integral and error: ", MyIntegral, "  ", IntegralError)
            print("Fitting pol1 integral and error: ", pol1.Integral(channels[j]-flexibility-300,channels[j]+flexibility+300))
            print("Fitting total integral and error: ", total.Integral( p1-3*p2, p1+3*p2 ))

            g1.Draw("same")
            pol1.Draw("same")
            total.Draw("same")

            xm1=0.6
            ym1=0.6
            xm2=0.8
            ym2=0.88
            legend1 = root.TLegend(xm1,ym1,xm2,ym2)
            legend1.SetHeader("^{152}Eu " + str(energys[j]) + " keV")
            legend1.AddEntry(g1,"Gauss Peak","l")
            legend1.AddEntry(pol1,"Polynomial Background","l")
            legend1.AddEntry(total,"Total Fit","l")
            legend1.SetBorderSize(0)
            legend1.SetTextSize(0.03)
            legend1.Draw()
            c.Update()
            c.SaveAs(plotsDirectory+"/gauss_fits/L" + str(i) + "_AmBe_" + str(energys[j]) + ".png")
            c.Close()
            fitParams[i][j].append(p1)
            fitParams[i][j].append(p2)
            print("fitParams: ", fitParams)
            
    return fitParams

def getEnergyResolution(fitParams):
    energyResolutions = []
    for i in range((0,4)):
        energyResolutions.append([])
        x = []
        y = []
        for j in range(0,13):
            energyResolutions[i].append([])
            energyResolutions[i][j].append((fitParams[i][j][1]))
            energyResolutions[i][j].append(((2.35*fitParams[i][j][2])/fitParams[i][j][1])*100) # 2.355 * sigma / mean 
            x.append((fitParams[i][j][1]))
            y.append(((2.35*fitParams[i][j][2])/fitParams[i][j][1])*100)
        c = root.TCanvas("c", "c", 1000,800)
        gr = root.TGraph(len(x), array.array('d', x), array.array('d', y))
        gr.SetTitle("Energy Resolution of Detector L" + str(i))
        gr.GetXaxis().SetTitle("Energy (keV)")
        gr.GetYaxis().SetTitle("Energy Resolution (%)")
        gr.Draw("AP")
        #increase the point size
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(1)
        gr.SetMarkerColor(root.kBlue)

        #legend1 = root.TLegend()
        #legend1.AddEntry(gr,"1st order polynomial fit","l")
        ##legend1.SetBorderSize(0)
        #legend1.SetTextSize(0.03)
        #legend1.Draw()
        c.Update()
        c.SaveAs(plotsDirectory+"/energy_resolution/L" + str(i) + ".png")
        c.Close()
    return energyResolutions


def getCalibration(fitParams, energyResolutions):
    # channel vs energy calibration for each detector where the error on channel total error on channel (x)
    # is error on background + sigma
    calibration = []

    return calibrations, residuals
    
def main():
    hists_152Eu, hists_AmBe = getHistos(f,g)
    fitParams = getFitParams(hists_152Eu, hists_AmBe)
    energyResolutions = getEnergyResolution(fitParams)
    calibrations, residuals = getCalibration(fitParams, energyResolutions)

# run the main function
if __name__ == "__main__":
    main()
