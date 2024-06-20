#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------
from cmath import nan
from ctypes import sizeof
from ftplib import parse150
from statistics import mean
from tokenize import Double
from unittest import skip
from venv import create
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 32})
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
from scipy.optimize import curve_fit

f = root.TFile.Open(sys.argv[1], "READ")
g = root.TFile.Open(sys.argv[2], "READ")
saveDirectory = sys.argv[1].split("/")[0:-1]
# if calibration_plots doesnt exist in saveDirectory, create it
#plotsDirectory = "/".join(saveDirectory)+"/calibration_plots"
plotsDirectory = "/home/shan/Documents/PhD/ExperimentResults/2022/Calibration/calibration_plots"
if not os.path.exists(plotsDirectory): # 
    os.makedirs(plotsDirectory)
    os.makedirs(plotsDirectory+"/gauss_fits")
    os.makedirs(plotsDirectory+"/energy_resolution")
    os.makedirs(plotsDirectory+"/energy_calibration")
# ----------------------------------------------------------------------------------------
# PEAK ENERGY
# ----------------------------------------------------------------------------------------
Energy = [121.81, 244.7, 344.3, 411.1, 444, 511, 778.9, 867.4, 964.1, 1408, 1470 , 3416, 3927, 4438, 2223.245, 2614.522]
EnergyL2 = [121.81, 244.7, 344.3, 411.1, 444, 511, 778.9, 867.4, 964.1, 1408, 1470, 2223.245, 2614.522, 3416, 3927,4438]

XRange = np.arange(0, 4500, 1)
CRange = np.arange(0, 8500, 1)
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L0
# ----------------------------------------------------------------------------------------
gaussAmpL0 = [38768 , 101066, 112799, 25432.7, 24914.2, 29304.3, 23437.3, 15163.5, 19678.5, 12617.2, 13882.5, 212.429 , 284.478 , 226.4640, 1034.82, 155.039]
gaussMeanL0 = [160.608 , 323.108 , 457.634 , 545.25, 589.187 , 680.668 , 1037.59, 1146.52, 1279.03, 1875.27, 1935.48, 4507.07, 5165.7, 5813.56, 2932.62, 3446.14]
gaussSigmaL0 = [9.46817 , 18.6455, 16.0223, 35.1651, 30.3395, 23.5746, 28.8243, 44.2992, 28.0501, 33.3803, 35.4086, 131.545 , 118.084 , 86.3007, 39.3711, 95.3512 ]
gaussMeanErrL0 = [0.00432861, 0.0186127 , 0.00897532, 0.0844736 , 0.072568, 0.0395102 , 0.0447314 , 0.17315 , 0.056759, 0.15068 , 0.120097, 1.65867 , 0.938417, 0.517557,0.163427, 1.29540 ]
gaussSigmaErrL0 = [0.00555433, 0.0310481 , 0.0103989 , 0.236424, 0.165381, 0.0723074 , 0.0802608 , 0.448676, 0.0968833 , 0.191186, 0.210336, 2.92745 , 1.83134 , 0.550225, 0.195005, 2.79274]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L0
# ----------------------------------------------------------------------------------------
mL0 = [310805.0 , 198182, 157102966, 14777.7, 92853.2, 28515.7, -11599.4, 40546, 56014.4, -104224.000 , 259079,-471.366 , -575.276 , 1364.17, 347.744, 367.124]
cL0 = [-943.881 , -406.153 , -145.737 , 15.0647, -119.763 , -13.7095, 25.2703, -24.0164, -35.4554, 61.4678, -127.498 , 0.140699, 0.150104, -0.220050, 0.0375189, -0.0673863]
mErrL0 = [189.473 , 275.541 , 141.228 , 567.236 , 585.807 , 258.504 , 201.383 , 234.906 , 157.483 , 1070.77, 891.762 ,20.8539, 25.8904, 14.4619, 79.562, 59.0168]
cErrL0 = [1.05475 , 0.844916, 0.297692, 1.04639 , 0.973431, 0.38171 , 0.195522, 0.203529, 0.121908, 0.572608, 0.455579, 0.00464896, 0.00505573, 0.00244120, 0.0271421, 0.0170952]


# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L1
# ----------------------------------------------------------------------------------------
gaussAmpL1 = [323430, 103807, 110965, 25363.4, 25114.5, 26455.6, 22140.2, 14997.90000000, 18691.1, 11653.7, 14603.2, 390.95, 524.866 , 424.8950, 1856.37, 251.876 ]
gaussMeanL1 = [149.599 , 302.029 , 425.452 , 507.292 , 547.778 , 632.156 , 962.781 , 1.06144e+03, 1185.38, 1739.17, 1794.15, 4191.13, 4791.47, 5389.37, 2721.87, 3191.39]
gaussSigmaL1 = [11.8851, 17.7681, 14.3319, 30.8899, 27.1413, 31.242 , 28.3204, 45.08500000, 26.8528, 32.4115, 33.1145, 151.343 , 115.221 , 83.6745, 40.4547, 113.764 ]
gaussMeanErrL1 = [0.00595092, 0.0228935 , 0.0103771 , 0.108486, 0.102783, 0.0573967 , 0.0537226 , 2.23883e-01, 0.0551473 , 0.209769, 0.119576, 1.21638 , 0.671244, 0.372742, 0.110805, 1.53115 ]
gaussSigmaErrL1 = [0.00824569, 0.0401211 , 0.013516, 0.361274, 0.359511, 0.126749, 0.0853805 , 6.64104e-01, 0.107345, 0.32069 , 0.300464, 1.81419 , 1.25165 , 0.444967, 0.123894, 4.26946 ]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L1
# ----------------------------------------------------------------------------------------
mL1 = [258671, 204758, 157301, 10446.300 ,94966.6, 44175.4, -14203.2, 52155.7, 88489.2, -187249, 123039, -1914.66, -462.57, 835.8730, 113.13, 317.968 ]
cL1 = [-653.937 , -446.7, -263.505 , 25.3397, -131.701 , -38.2299, 29.9405, -36.4701, -64.4047, 113.672 , -61.8742, 0.538685, 0.173 , -0.121443, 0.282326, -0.0269795]
mErrL1 = [248.704 ,332.264 , 346.04, 572.106 , 677.233 , 286.987 , 260.202 , 288.077 , 275.898 , 798.436 , 950.423 , 46.4773, 36.1925, 26.3084, 117.001, 82.2105 ]
cErrL1 = [1.63941, 1.08555, 0.794903, 1.14313, 1.21289, 0.448317, 0.269195, 0.267073, 0.22768, 0.463348, 0.527991, 0.0112324, 0.00756592, 0.00487826, 0.0429584, 0.0257422]

# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L2
# ----------------------------------------------------------------------------------------
gaussAmpL2 = [241422 ,65337.2, 77412.6, 16441.6, 16214.4, 19584.3, 14980.5, 9436.57, 12142.2, 7896.6, 9559.21, 396.327, 73.4288,81.2167,107.094,89.2483]
gaussMeanL2 = [234.116 , 471.254 , 666.923 , 794.766 ,859.126 ,988.201 , 1509.3, 1670.66, 1858.94, 2723.28, 2814.38, 4257.20, 4992.60, 6561.41,7499.25,8439.72]
gaussSigmaL2 = [14.38210000 ,25.5802, 19.6496, 44.4215, 39.8437, 33.1927, 36.796 , 60.5302, 40.5241, 44.2237 ,44.3424, 42.5720, 128.598, 222.939,182.446,124.555]
gaussMeanErrL2 = [0.00758056, 0.0308865 , 0.0150532 , 0.136702, 0.123888, 0.0602735 , 0.0728618 , 0.268571, 0.0746658 , 0.332815, 0.174987, 0.301183, 1.88251,2.90916,1.84703,0.888188]
gaussSigmaErrL2 = [0.0107252 , 0.0552928 , 0.0202322 , 0.41612 , 0.306796, 0.108268, 0.130866, 0.668506, 0.132893, 0.456536, 0.274263, 0.349926, 1.88251, 2.90916,1.84703,0.902082]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L2
# ----------------------------------------------------------------------------------------
mL2 = [207298, 149163, 108738, 12248.5 ,52276.4, 34406.7, -12791.1, 24967.1, 42065, -57350.7, 51831.8, 309.718, 193.537, -173.456, -28.9374, 220.113]
cL2 = [-435.545, -222.889, -121.892, 3.12652, -44.3865, -22.5422, 14.2698, -10.1473, -18.9637, 23.6345, -15.5607, -0.0488323, -0.0284277, 0.0354905, -28.9374, 220.113]
mErrL2 = [164.545, 208.271, 188.652, 373.187, 392.127, 146.156, 124.199, 141.686, 133.158, 764.557, 826.59, 8.45636, 9.33002, 5.53396, 1.77603, 0.896688]
cErrL2 = [0.662718, 0.426544, 0.272908, 0.471908, 0.450031, 0.145507, 0.0832166, 0.0831426, 0.0705606, 0.280658, 0.293641, 8.45636, 0.00185804, 0.000859117, 1.77603, 0.896688]

# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L3
# ----------------------------------------------------------------------------------------
gaussAmpL3 = [289739, 88764.3, 97416.5, 23943.8, 23305.4, 29827.9, 19854,13615.1, 16535.2, 10324.3, 12027.3, 58.8438, 77.125 , 60.2518, 27.6556, 72.5900 ]
gaussMeanL3 = [168.412 , 339.197 , 482.045 , 574.9, 619.126 , 714.041 , 1091.28, 1205.01, 1343, 1970.99, 2034.45, 4748.84, 5439.59, 6122.54, 3086.72, 3628.36]
gaussSigmaL3 = [14.8908, 22.8161, 17.7707, 41.7924, 38.0786, 28.7329, 33.7265, 50.8544, 34.1021, 38.1912, 39.2843, 137.839 , 133.008, 97.2427, 49.7551, 85.678]
gaussMeanErrL3 = [0.00821219, 0.0388648 , 0.0149587 , 0.159996, 0.168495, 0.0494714 , 0.10038 , 0.43267 , 0.0832548 , 0.295917, 0.201815, 3.10867 , 1.77808 , 1.1153, 0.379352, 1.58747]
gaussSigmaErrL3 = [0.0131817 , 0.0765902 , 0.02330370 ,0.648552, 0.501309, 0.10235600 ,0.174357, 1.53738 , 0.19439 ,0.493528, 0.343385, 5.83268 , 3.07196 ,1.349550, 0.485720, 2.96208]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L3
# ----------------------------------------------------------------------------------------
mL3 = [31953 , 186501, 107659, 14237.4, 78259.2, 57329.9, -32542.4, 41242.9, 59829.9, -145741, 45587, -186.838, -383.555, 335.205, -633.888, 290.833]
cL3 = [-944.519, -362.04, -143.738, 14.5566, -90.9263, -51.5461, 42.9951, -23.9671, -36.5185, 78.7147, -17.2289, 0.0490386, 0.0819119, -0.0503129, 0.262212, -0.0633561]
mErrL3 = [256.741, 182.493, 185.097, 662.072, 720.087,283.729, 244.234, 349.821, 204.132, 819.741, 969.665, 10.7665, 21.4641, 8.32757, 63.6361, 40.4554]
cErrL3 = [1.40382, 0.524332, 0.375088, 1.15936, 1.14492, 0.392508, 0.227195, 0.288529, 0.149852, 0.419187, 0.476942, 0.00228389, 0.00398461, 0.00133828, 0.0206613, 0.0111115]

# ----------------------------------------------------------------------------------------
# READ IN ROOT HISTOGRAMS FOR BACKGROUND SIGMA CALCULATION
# ----------------------------------------------------------------------------------------

def getHistos(f,g):
    hists_152Eu = []
    hists_AmBe = []
    for i in range(0,4):
        hists_152Eu.append(f.Get("L"+str(i)+"_Uncalibrated"))
        hists_AmBe.append(g.Get("L"+str(i)+"_Uncalibrated"))
    return hists_152Eu, hists_AmBe

def getBackgroundErrFitParams(hists_152Eu, hists_AmBe):
    for i in range(0,4):
            means = eval("gaussMeanL"+str(i))
            #print("MEANS: ",means)
            amp = eval("gaussAmpL"+str(i))
            sigma = eval("gaussSigmaL"+str(i))
            meanserr = eval("gaussMeanErrL"+str(i))
            sigmaerr = eval("gaussSigmaErrL"+str(i))
            polM = eval("mL"+str(i))
            polC = eval("cL"+str(i))
            merr = eval("mErrL"+str(i))
            cerr = eval("cErrL"+str(i))
            
            # 152Eu Part #
            for j in range(0,11):
                c = root.TCanvas("c", "c", 1000 ,800)          
                # draw a gaussian fit at every means value using the sigma and amp values from the fit
                # and the mean and sigma errors from the fit
                # the mean and sigma errors are means[j] +/- meanserr[j] and sigma[j] +/- sigmaerr[j]
                # the amp is the amp[j] value from the fit
                g1 = root.TF1("g1","gaus",means[j]-1*sigma[j],means[j]+1*sigma[j])
                g1.SetParameters(amp[j],means[j],sigma[j])
                pol1 = root.TF1("pol1","pol1",means[j]-1*sigma[j],means[j]+1*sigma[j]) # linear background
                pol1.SetParameters(polC[j],polM[j])
                total = root.TF1( "total", "gaus(0) + pol1(3)", means[j]-1*sigma[j],means[j]+1*sigma[j]) # total func
                h1 = hists_152Eu[i].Clone()
                h1.SetLineColor(root.kBlack)
                g1.SetLineColor(root.kBlue)
                pol1.SetLineColor(root.kGreen)
                total.SetLineColor(root.kRed)

                h1.Fit(g1, "R")
                h1.Fit(pol1, "R+")
                h1.SetStats(000000000)

                p0 = amp[j]
                p1 = means[j]
                p2 = sigma[j]

                print("amplitude value from gauss fit= ", p0 ,"\nmean value from gauss fit= ", p1, "\nsigma value from gauss fit= ", p2, "\nFWHM= ", 2.355*p2)

                p3 = polM[j]
                p4 = polC[j]
                
                e0 = g1.GetParError(0) # error on amplitude
                e1 = meanserr[j] # error on mean
                e2 = sigmaerr[j] # error on sigma
                e3 = merr[j] # error on background
                e4 = cerr[j] # error on background

                par = [p0 , p1, p2, p3, p4] 
                h1.Fit(total, "R+")

                gaussParameters = [p0 , p1, p2]
                gaussParameters = array.array('d', gaussParameters)
                MyConvarienceMatrix = [e0 ,0 ,0 ,0 ,e1, 0,0 ,0 ,e2]
                MyConvarienceMatrix = array.array('d', MyConvarienceMatrix)

                HistoBinning = 1
                MyIntegral = g1.Integral( p1-1*p2, p1+1*p2 ) / HistoBinning #p4 is mean value and p5 is sigma value. Gauss range is about 8sigma
                IntegralError = g1.IntegralError(p1-1*p2, p1+1*p2, gaussParameters, MyConvarienceMatrix) / HistoBinning
                
                #print("Fitting g1 integral and error: ", MyIntegral, "  ", IntegralError)
                #print("Fitting pol1 integral and error: ", pol1.Integral(means[j]-1*sigma[j],means[j]+1*sigma[j]))
                #print("Fitting total integral and error: ", total.Integral( p1-1*p2, p1+1*p2 ))
                
                h1.GetXaxis().SetRangeUser(means[j]-500,means[j]+500)

                g1.Draw("same")
                pol1.Draw("same")
                total.Draw("same")

                xm1=0.6
                ym1=0.6
                xm2=0.8
                ym2=0.88
                legend1 = root.TLegend(xm1,ym1,xm2,ym2)
                legend1.SetHeader("^{152}Eu " + str(Energy[j]) + " keV")
                legend1.AddEntry(g1,"Gauss Peak","l")
                legend1.AddEntry(pol1,"Polynomial Background","l")
                legend1.AddEntry(total,"Total Fit","l")
                legend1.SetBorderSize(0)
                legend1.SetTextSize(0.03)
                legend1.Draw()
                c.Update()
                c.SaveAs(plotsDirectory+"/gauss_fits/L" + str(i) + "_152Eu_" + str(Energy[j]) + ".png")
                c.Close()

            # AmBe Part #
            for j in range(11,16):                             
                c = root.TCanvas("c", "c", 1000 ,800)
                g1 = root.TF1("g1","gaus",means[j]-1*sigma[j],means[j]+1*sigma[j])
                g1.SetParameters(amp[j],means[j],sigma[j])
                pol1 = root.TF1("pol1","pol1",means[j]-1*sigma[j],means[j]+1*sigma[j]) # linear background
                pol1.SetParameters(polC[j],polM[j])
                total = root.TF1( "total", "gaus(0) + pol1(3)", means[j]-1*sigma[j],means[j]+1*sigma[j]) # total func
                h1 = hists_AmBe[i].Clone()
                h1.SetLineColor(root.kBlack)
                g1.SetLineColor(root.kBlue)
                pol1.SetLineColor(root.kGreen)
                total.SetLineColor(root.kRed)
                h1.Fit(g1, "R")
                h1.Fit(pol1, "R+")
                h1.SetStats(000000000)
                
                par = [p0 , p1, p2, p3, p4]
                h1.Fit(total, "R+")

                p0 = amp[j]
                p1 = means[j]
                p2 = sigma[j] # p2 is sigma of gaussian, width of gaussian=8sigma
                #print("amplitude value from gauss fit= ", p0 ,"\nmean value from gauss fit= ", p1, "\nsigma value from gauss fit= ", p2, "\nFWHM= ", 2.355*p2)
                
                p3 = polM[j]
                p4 = polC[j]
                
                e0 = g1.GetParError(0) # error on amplitude
                e1 = meanserr[j] # error on mean
                e2 = sigmaerr[j] # error on sigma

                e3 = merr[j] # error on background
                e4 = cerr[j] # error on slope


                gaussParameters = [p0 , p1, p2]
                gaussParameters = array.array('d', gaussParameters)
                MyConvarienceMatrix = [e0 ,0 , 0, 0,e1, 0,0 ,0 ,e2]
                MyConvarienceMatrix = array.array('d', MyConvarienceMatrix)

                h1.GetXaxis().SetRangeUser(means[j]-500,means[j]+500)

                HistoBinning = 1
                MyIntegral = g1.Integral( p1-1*p2, p1+1*p2 ) / HistoBinning #p4 is mean value and p5 is sigma value. Gauss range is about 8sigma
                IntegralError = g1.IntegralError(p1-1*p2, p1+1*p2, gaussParameters, MyConvarienceMatrix) / HistoBinning
                #print("Fitting g1 integral and error: ", MyIntegral, "  ", IntegralError)
                #print("Fitting pol1 integral and error: ", pol1.Integral(means[j]-1*sigma[j],means[j]+1*sigma[j]))
                #print("Fitting total integral and error: \n", total.Integral( p1-1*p2, p1+1*p2 ) , " +- " , total.IntegralError( p1-1*p2, p1+1*p2, gaussParameters, MyConvarienceMatrix))

                # save the total.IntegralError to a txt file called "IntegralError_L" + str(i) + "_AmBe_" + str(Energy[j]) + ".txt"


                g1.Draw("same")
                pol1.Draw("same")
                total.Draw("same")

                xm1=0.6
                ym1=0.6
                xm2=0.8
                ym2=0.88
                legend1 = root.TLegend(xm1,ym1,xm2,ym2)
                legend1.SetHeader("^{152}Eu " + str(Energy[j]) + " keV")
                legend1.AddEntry(g1,"Gauss Peak","l")
                legend1.AddEntry(pol1,"Polynomial Background","l")
                legend1.AddEntry(total,"Total Fit","l")
                legend1.SetBorderSize(0)
                legend1.SetTextSize(0.03)
                legend1.Draw()
                c.Update()
                c.SaveAs(plotsDirectory+"/gauss_fits/L" + str(i) + "_AmBe_" + str(Energy[j]) + ".png")
                c.Close()

def func(x, a, b):
    return np.dot(a, np.array(x)) + b

def func2(x, a, b, c):
    return np.dot(a, np.array(np.square(x))) + np.dot(b, np.array(x)) + c

def func3(x, a, b, c, d):
    return np.dot(a, np.array(np.power(x,3))) + np.dot(b, np.array(np.power(x,2))) + np.dot(c, np.array(x)) + d

def calibrationPlots():
    p1ParamsL0 = []
    p2ParamsL0 = []
    p1ParamsL1 = []
    p2ParamsL1 = []
    p1ParamsL3 = []
    p2ParamsL3 = []
    # for i in range 0,1 and 2,4
    for i in range(0,2):
        mean = eval("gaussMeanL" + str(i))
        sigmaerr = eval("gaussSigmaErrL" + str(i))
        
        plt.figure(figsize=(30,20))
        plot1 = plt.errorbar(mean, Energy, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func, mean, Energy, sigma=sigmaerr)
        if i == 0:
            p1ParamsL0.append(popt)
        elif i == 1:
            p1ParamsL1.append(popt)
        plt.plot(CRange, func(CRange, *popt), 'b-', label='fit: a=%10.8f, b=%10.8f' % tuple(popt)) # CRange is channel range

        plt.title("Energy Calibration for L" + str(i), fontsize=39)
        plt.ylabel("Energy (keV)", fontsize=32)
        plt.grid( linestyle='-', linewidth='0.5', color='gray')
        plt.xlabel("Channel Number", fontsize=32)
        plt.legend()

        plt.xticks(fontsize=32)
        plt.yticks(fontsize=32)
        plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + ".png")
        plt.close()


        # redo the above but with func2
        plt.figure(figsize=(30,20))
        plot2 = plt.errorbar(mean, Energy, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func2, mean, Energy, sigma=sigmaerr)
        if i == 0:
            p2ParamsL0.append(popt)
        elif i == 1:
            p2ParamsL1.append(popt)
        plt.plot(CRange, func2(CRange, *popt), 'b-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f' % tuple(popt)) 
        plt.title("Energy Calibration for L" + str(i), fontsize=39)
        plt.ylabel("Energy (keV)", fontsize=32)
        plt.grid( linestyle='-', linewidth='0.5', color='gray')
        plt.xlabel("Channel Number", fontsize=32)
        plt.legend()
        plt.xticks(fontsize=32)
        plt.yticks(fontsize=32)
        plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrder.png")
        plt.close()


    for i in range(3,4):
        mean = eval("gaussMeanL" + str(i))
        sigmaerr = eval("gaussSigmaErrL" + str(i))
        
        plt.figure(figsize=(30,20))
        plot1 = plt.errorbar(mean, Energy, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func, mean, Energy, sigma=sigmaerr)
        if i == 3:
            p1ParamsL3.append(popt)
        plt.plot(CRange, func(CRange, *popt), 'b-', label='fit: a=%10.8f, b=%10.8f' % tuple(popt))
        plt.title("Energy Calibration for L" + str(i), fontsize=39)
        plt.ylabel("Energy (keV)", fontsize=32)
        plt.grid( linestyle='-', linewidth='0.5', color='gray')
        plt.xlabel("Channel Number", fontsize=32)
        plt.legend()
        plt.xticks(fontsize=32)
        plt.yticks(fontsize=32)
        plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + ".png")
        plt.close()
    

        # redo the above but with func2
        plt.figure(figsize=(30,20))
        plot2 = plt.errorbar(mean, Energy, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func2, mean, Energy, sigma=sigmaerr)
        if i == 3:
            p2ParamsL3.append(popt)
        plt.plot(CRange, func2(CRange, *popt), 'b-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f' % tuple(popt))
        plt.grid( linestyle='-', linewidth='0.5', color='gray')
        plt.title("Energy Calibration for L" + str(i), fontsize=39)
        plt.ylabel("Energy (keV)", fontsize=32)
        plt.xlabel("Channel Number", fontsize=32)
        plt.legend()
        plt.xticks(fontsize=32)
        plt.yticks(fontsize=32)
        plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrder.png")
        plt.close()
        print("\np1ParamsL0: ", p1ParamsL0, "\np2ParamsL0: ", p2ParamsL0, "\np1ParamsL1: ", p1ParamsL1, "\np2ParamsL1: ", p2ParamsL1, "\np1ParamsL3: ", p1ParamsL3, "\np2ParamsL3: ", p2ParamsL3)
    return  p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3

def calibrationPlotsL2():
    p1ParamsL2 = []
    p2ParamsL2 = []
    # for i ==2 only
    for i in range(2,3):
        mean = eval("gaussMeanL" + str(i))
        sigmaerr = eval("gaussSigmaErrL" + str(i))
        
        plt.figure(figsize=(30,20))
        plot1 = plt.errorbar(mean, EnergyL2, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func, mean, EnergyL2, sigma=sigmaerr)
        p1ParamsL2.append(popt)
        plt.plot(CRange, func(CRange, *popt ), 'b-', label='fit: a=%10.8f, b=%10.8f' % tuple(popt))

        plt.title("Energy Calibration for L" + str(i), fontsize=39)
        plt.ylabel("Energy (keV)", fontsize=32)
        plt.grid( linestyle='-', linewidth='0.5', color='gray')
        plt.xlabel("Channel Number", fontsize=32)
        plt.legend()
        plt.xticks(fontsize=32)
        plt.yticks(fontsize=32)
        plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + ".png")
        plt.close()

        # redo the above but with func2
        plt.figure(figsize=(30,20))
        plot2 = plt.errorbar(mean, EnergyL2, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func2, mean, EnergyL2, sigma=sigmaerr)
        p2ParamsL2.append(popt)
        plt.plot(CRange, func2(CRange, *popt), 'b-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f' % tuple(popt))
        plt.grid( linestyle='-', linewidth='0.5', color='gray')
        plt.title("Energy Calibration for L" + str(i), fontsize=39)
        plt.ylabel("Energy (keV)", fontsize=32)
        plt.xlabel("Channel Number", fontsize=32)
        plt.legend()
        plt.xticks(fontsize=32)
        plt.yticks(fontsize=32)
        plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrder.png")
        plt.close()
        print("\np1ParamsL2: ", p1ParamsL2, "\np2ParamsL2: ", p2ParamsL2)
    return  p1ParamsL2, p2ParamsL2

def calculateResiduals(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3):
        for i in range(0,2):
            # ensure that the correct i value is used for the correct L value of p1Params and p2Params, i.e. L0 for i=0 and L1 for i=1
            mean = eval("gaussMeanL" + str(i))
            sigmaerr = eval("gaussSigmaErrL" + str(i))
            p1ParamsL = eval("p1ParamsL" + str(i))
            p2ParamsL = eval("p2ParamsL" + str(i))
            residuals = Energy - func(mean, *p1ParamsL[0])
            plt.figure(figsize=(30,20))
            plt.errorbar(Energy, residuals, yerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3,  capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Residuals for L" + str(i) + " 1$^{st}$ Order Polynomial Fit", fontsize=39)
            plt.ylabel("Residuals (keV)", fontsize=32)
            plt.xlabel("Energy (keV)", fontsize=32)
            plt.xticks(fontsize=32)
            plt.yticks(fontsize=32)
            plt.ylim(-10,70)
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_1stOrderResiduals.png")
            plt.close()

            residuals = Energy - func2(mean, *p2ParamsL[0])
            plt.figure(figsize=(30,20))
            plt.errorbar(Energy, residuals, yerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Residuals for L" + str(i) + " 2$^{nd}$ Order Polynomial Fit", fontsize=39)
            plt.ylabel("Residuals (keV)", fontsize=32)
            plt.xlabel("Energy (keV)", fontsize=32)
            plt.xticks(fontsize=32)
            plt.yticks(fontsize=32)
            plt.ylim(-30,15)
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrderResiduals.png")
            plt.close()


        for i in range(3,4):
            # for the first Order Polynomial Fit
            mean = eval("gaussMeanL" + str(i))
            sigmaerr = eval("gaussSigmaErrL" + str(i))
            p1ParamsL = eval("p1ParamsL" + str(i))
            p2ParamsL = eval("p2ParamsL" + str(i))

            # calculate residuals for 1st Order Polynomial Fit
            residuals = Energy - func(mean, *p1ParamsL[0])
            plt.figure(figsize=(30,20))
            plt.errorbar(Energy, residuals, yerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Residuals for L" + str(i) + " 1$^{st}$ Order Polynomial Fit", fontsize=39)
            plt.ylabel("Residuals (keV)", fontsize=32)
            plt.xlabel("Energy (keV)", fontsize=32)
            plt.xticks(fontsize=32)
            plt.ylim(-10,70)
            plt.yticks(fontsize=32)
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_1stOrderResiduals.png")
            plt.close()

            # for the second Order Polynomial Fit
            residuals = Energy - func2(mean, *p2ParamsL[0])
            plt.figure(figsize=(30,20))
            plt.errorbar(Energy, residuals, yerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Residuals for L" + str(i) + " 2$^{nd}$ Order Polynomial Fit", fontsize=39)
            plt.ylabel("Residuals (keV)", fontsize=32)
            plt.xlabel("Energy (keV)", fontsize=32)
            plt.xticks(fontsize=32)
            plt.yticks(fontsize=32)
            plt.ylim(-30,15)
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrderResiduals.png")
            plt.close()

def calculateResidualsL2(p1ParamsL2, p2ParamsL2):
        for i in range(2,3):
            mean = eval("gaussMeanL" + str(i))
            sigmaerr = eval("gaussSigmaErrL" + str(i))
            p1ParamsL = eval("p1ParamsL" + str(i))
            p2ParamsL = eval("p2ParamsL" + str(i))
            residuals = EnergyL2 - func(mean, *p1ParamsL[0])
            plt.figure(figsize=(30,20))
            plt.errorbar(EnergyL2, residuals, yerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Residuals for L" + str(i) + " 1$^{st}$ Order Polynomial Fit", fontsize=39)
            plt.ylabel("Residuals (keV)", fontsize=32)
            plt.xlabel("Energy (keV)", fontsize=32)
            plt.xticks(fontsize=32)
            plt.ylim(-10,70)
            plt.yticks(fontsize=32)
            
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_1stOrderResiduals.png")
            plt.close()

            # for the second Order Polynomial Fit
            residuals = EnergyL2 - func2(mean, *p2ParamsL[0])
            plt.figure(figsize=(30,20))
            plt.errorbar(EnergyL2, residuals, yerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Residuals for L" + str(i) + " 2$^{nd}$ Order Polynomial Fit", fontsize=39)
            plt.ylabel("Residuals (keV)", fontsize=32)
            plt.xlabel("Energy (keV)", fontsize=32) 
            plt.xticks(fontsize=32)
            plt.ylim(-30,15)
            plt.yticks(fontsize=32)
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrderResiduals.png")
            plt.close()

def energyResolution(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2):
        for i in range(0,4):
            mean = eval("gaussMeanL" + str(i))
            meanerr = eval("gaussMeanErrL" + str(i))
            p1ParamsL = eval("p1ParamsL" + str(i))
            energy = func(mean, *p1ParamsL[0])
            sigma = eval("gaussSigmaL" + str(i))
            sigmaerr = eval("gaussSigmaErrL" + str(i))
            energy = func (mean, *p1ParamsL[0])
            energyResolution = ((np.dot(2.355, sigma))/energy)*100
            energyResolutionErr = energyResolution*np.sqrt((2*np.array(sigmaerr)/np.array(sigma))**2 + (np.array(meanerr)/np.array(mean))**2)


            # plot energy resolution vs. energy
            plt.figure(figsize=(30,20))
            plt.errorbar(energy, energyResolution, yerr=energyResolutionErr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Energy Resolution for L" + str(i), fontsize=39)
            plt.ylabel("Energy Resolution [%]", fontsize=32)
            plt.xlabel("Energy [keV]", fontsize=32)
            plt.xticks(fontsize=32)
            plt.yticks(fontsize=32)
            popt, pcov = curve_fit(func2, energy, energyResolution, sigma=energyResolutionErr, absolute_sigma=True)
            plt.plot(XRange, func2(XRange, *popt), 'r-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f' % tuple(popt))
            plt.legend(loc='best')
            energyResolution662 = func2(662, *popt)
            plt.text(662, energyResolution662, "Energy Resolution at 662 keV = " + str(round(energyResolution662, 2)) + "%", fontsize=32)
            plt.savefig(plotsDirectory+"/energy_resolution/L" + str(i) + "_energyResolution_2ndOrder.png")
            plt.close()

            # do the same again but with func3
            plt.figure(figsize=(30,20))
            plt.errorbar(energy, energyResolution, yerr=energyResolutionErr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            plt.grid( linestyle='-', linewidth='0.5', color='gray')
            plt.title("Energy Resolution for L" + str(i), fontsize=39)
            plt.ylabel("Energy Resolution [%]", fontsize=32)
            plt.xlabel("Energy [keV]", fontsize=32)
            plt.xticks(fontsize=32)
            plt.yticks(fontsize=32)
            popt, pcov = curve_fit(func3, energy, energyResolution, sigma=energyResolutionErr, absolute_sigma=True)
            plt.plot(XRange, func3(XRange, *popt), 'r-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f, d=%10.8f' % tuple(popt))
            plt.legend(loc='best')
            energyResolution662 = func3(662, *popt)
            plt.text(662, energyResolution662, "Energy Resolution at 662 keV = " + str(round(energyResolution662, 2)) + "%", fontsize=32)
            plt.savefig(plotsDirectory+"/energy_resolution/L" + str(i) + "_energyResolution_3rdOrder.png")
            plt.close()


            
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN CODE
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():
    #hists_152Eu, hists_AmBe = getHistos(f,g)
    #getBackgroundErrFitParams(hists_152Eu, hists_AmBe)
    p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3 = calibrationPlots()
    p1ParamsL2, p2ParamsL2 = calibrationPlotsL2()

    calculateResiduals(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3)
    calculateResidualsL2(p1ParamsL2, p2ParamsL2)

    energyResolution(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2)

    # make a 2,2 plot of the 1stOrderResiduals for L0, L1, L2, L3
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(30,20))
    image0 = plt.imread(plotsDirectory+"/energy_calibration/L0_1stOrderResiduals.png")
    axs[0,0].imshow(image0)
    axs[0,0].axis('off')
    image1 = plt.imread(plotsDirectory+"/energy_calibration/L1_1stOrderResiduals.png")
    axs[0,1].imshow(image1)
    axs[0,1].axis('off')
    image2 = plt.imread(plotsDirectory+"/energy_calibration/L2_1stOrderResiduals.png")
    axs[1,0].imshow(image2)
    axs[1,0].axis('off')
    image3 = plt.imread(plotsDirectory+"/energy_calibration/L3_1stOrderResiduals.png")
    axs[1,1].imshow(image3)
    axs[1,1].axis('off')
    axs[0,0].set_ylabel("Residuals [keV]", fontsize=32)
    axs[1,0].set_ylabel("Residuals [keV]", fontsize=32)
    axs[1,0].set_xlabel("Energy [keV]", fontsize=32)
    axs[1,1].set_xlabel("Energy [keV]", fontsize=32)
    plt.savefig(plotsDirectory+"/energy_calibration/GRID_1stOrderResiduals.png")

    # make a 2,2 plot of the 2ndOrderResiduals for L0, L1, L2, L3
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(30,20))
    image0 = plt.imread(plotsDirectory+"/energy_calibration/L0_2ndOrderResiduals.png")
    axs[0,0].imshow(image0)
    axs[0,0].axis('off')
    image1 = plt.imread(plotsDirectory+"/energy_calibration/L1_2ndOrderResiduals.png")
    axs[0,1].imshow(image1)
    axs[0,1].axis('off')
    image2 = plt.imread(plotsDirectory+"/energy_calibration/L2_2ndOrderResiduals.png")
    axs[1,0].imshow(image2)
    axs[1,0].axis('off')
    image3 = plt.imread(plotsDirectory+"/energy_calibration/L3_2ndOrderResiduals.png")
    axs[1,1].imshow(image3)
    axs[1,1].axis('off')
    axs[0,0].set_ylabel("Residuals [keV]", fontsize=32)
    axs[1,0].set_ylabel("Residuals [keV]", fontsize=32)
    axs[1,0].set_xlabel("Energy [keV]", fontsize=32)
    axs[1,1].set_xlabel("Energy [keV]", fontsize=32)
    plt.savefig(plotsDirectory+"/energy_calibration/GRID_2ndOrderResiduals.png")

    # make a 2,2 plot of the Energy Resolution for L0, L1, L2, L3
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(30,20))
    image0 = plt.imread(plotsDirectory+"/energy_resolution/L0_energyResolution_2ndOrder.png")
    axs[0,0].imshow(image0)
    axs[0,0].axis('off')
    image1 = plt.imread(plotsDirectory+"/energy_resolution/L1_energyResolution_2ndOrder.png")
    axs[0,1].imshow(image1)
    axs[0,1].axis('off')
    image2 = plt.imread(plotsDirectory+"/energy_resolution/L2_energyResolution_2ndOrder.png")
    axs[1,0].imshow(image2)
    axs[1,0].axis('off')
    image3 = plt.imread(plotsDirectory+"/energy_resolution/L3_energyResolution_2ndOrder.png")
    axs[1,1].imshow(image3)
    axs[1,1].axis('off')
    axs[0,0].set_ylabel("Energy Resolution [keV]", fontsize=32)
    axs[1,0].set_ylabel("Energy Resolution [keV]", fontsize=32)
    axs[1,0].set_xlabel("Energy [keV]", fontsize=32)
    axs[1,1].set_xlabel("Energy [keV]", fontsize=32)
    plt.savefig(plotsDirectory+"/energy_resolution/GRID_energyResolution_2ndOrder.png")

    # make a 2,2 plot of the Energy Calibration for L0, L1, L2, L3
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(30,20))
    image0 = plt.imread(plotsDirectory+"/energy_calibration/L0.png")
    axs[0,0].imshow(image0)
    axs[0,0].axis('off')
    image1 = plt.imread(plotsDirectory+"/energy_calibration/L1.png")
    axs[0,1].imshow(image1)
    axs[0,1].axis('off')
    image2 = plt.imread(plotsDirectory+"/energy_calibration/L2.png")
    axs[1,0].imshow(image2)
    axs[1,0].axis('off')
    image3 = plt.imread(plotsDirectory+"/energy_calibration/L3.png")
    axs[1,1].imshow(image3)
    axs[1,1].axis('off')
    axs[0,0].set_ylabel("Energy [keV]", fontsize=32)
    axs[1,0].set_ylabel("Energy [keV]", fontsize=32)
    axs[1,0].set_xlabel("Channel", fontsize=32)
    axs[1,1].set_xlabel("Channel", fontsize=32)
    plt.savefig(plotsDirectory+"/energy_calibration/GRID_EnergyCalib1stOrder.png")

    # make a 2,2 plot of the Energy Calibration for L0, L1, L2, L3
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(30,20))
    image0 = plt.imread(plotsDirectory+"/energy_calibration/L0_2ndOrder.png")
    axs[0,0].imshow(image0)
    axs[0,0].axis('off')
    image1 = plt.imread(plotsDirectory+"/energy_calibration/L1_2ndOrder.png")
    axs[0,1].imshow(image1)
    axs[0,1].axis('off')
    image2 = plt.imread(plotsDirectory+"/energy_calibration/L2_2ndOrder.png")
    axs[1,0].imshow(image2)
    axs[1,0].axis('off')
    image3 = plt.imread(plotsDirectory+"/energy_calibration/L3_2ndOrder.png")
    axs[1,1].imshow(image3)
    axs[1,1].axis('off')
    axs[0,0].set_ylabel("Energy [keV]", fontsize=32)
    axs[1,0].set_ylabel("Energy [keV]", fontsize=32)
    axs[1,0].set_xlabel("Channel", fontsize=32)
    axs[1,1].set_xlabel("Channel", fontsize=32)
    plt.savefig(plotsDirectory+"/energy_calibration/GRID_EnergyCalib2ndOrder.png")




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
if __name__ == "__main__":
    main()