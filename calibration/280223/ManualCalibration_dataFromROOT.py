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
saveDirectory = sys.argv[1].split("/")[0:-1]
# if calibration_plots doesnt exist in saveDirectory, create it
#plotsDirectory = "/".join(saveDirectory)+"/calibration_plots"
plotsDirectory = "/home/shan/Documents/PhD/ExperimentResults/2023/230228/calibrationPlots"
if not os.path.exists(plotsDirectory): # 
    os.makedirs(plotsDirectory)
    os.makedirs(plotsDirectory+"/gauss_fits")
    os.makedirs(plotsDirectory+"/energy_resolution")
    os.makedirs(plotsDirectory+"/energy_calibration")
# ----------------------------------------------------------------------------------------
# PEAK ENERGY
# ----------------------------------------------------------------------------------------
 # origin: 22Na, 137Cs, 60Co, 22Na, 60Co, LaBr3
Energy = [511, 661.7, 1173.2, 1274.5, 1461] #keV

XRange = np.arange(0, 4500, 1)
CRange = np.arange(0, 8500, 1)
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L0
# ----------------------------------------------------------------------------------------
gaussAmpL0 = [1.20331e+04, 2.76116e+03, 9.29030e+02, 2.86056e+03,  4.54538e+02]
gaussMeanL0 =  [2.81122e+02, 3.63878e+02, 6.41776e+02, 6.96597e+02,  7.98512e+02]
gaussSigmaL0 = [2.81122e+02, 9.62880e+00, 1.28127e+01, 1.02426e+01,  1.52607e+01]
gaussMeanErrL0 = [1.70120e-02, 5.09876e-02, 1.09979e-01, 4.59728e-02,  1.28397e-01]
gaussSigmaErrL0 = [1.74488e-02, 6.48328e-02, 1.49285e-01, 4.63478e-02,  1.24923e-01]
# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L0
# ----------------------------------------------------------------------------------------
mL0 = [7821.55, 1688.11, 1010.83, -4453.34, 1370.88]
cL0 = [-20.3754, -1.74934, -0.982229, 7.47891, -1.56424]
mErrL0 = [35.7981, 5.77263, 1.93404, 87.7996,  27.6071]
cErrL0 = [0.125689, 0.00980303, 0.00286872, 0.126593, 0.0336021]


# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L1
# ----------------------------------------------------------------------------------------
gaussAmpL1 =  [4.98801e+04, 7.42290e+03, 1.97363e+03, 8.01436e+03, 6.04554e+02]
gaussMeanL1 = [1.94005e+02, 2.50845e+02, 4.40042e+02, 4.78339e+02, 5.47946e+02]
gaussSigmaL1 = [5.75739e+00, 7.83220e+00, 9.23967e+00, 8.30749e+00, 1.22772e+01]
gaussMeanErrL1 = [6.91850e-03, 2.40293e-02, 5.70580e-02, 2.40948e-02, 1.06273e-01]
gaussSigmaErrL1 = [6.28373e-03, 2.67651e-02, 7.01266e-02, 2.29503e-02, 1.09683e-01]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L1
# ----------------------------------------------------------------------------------------
mL1 = [6686.16, 5977.79, 1706.47, -213.241, 625.069]
cL1 = [-13.7765, -12.1168, -2.08807, 5.37668, -0.851133]
mErrL1 = [8.45893, 8.30936, 3.069, 0.970203, 1.61802]
cErrL1 = [0.0251954, 0.0268025, 0.00691405, 0.00789602, 0.00295028]

# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L2
# ----------------------------------------------------------------------------------------
gaussAmpL2 = [2.56490e+04, 3.84317e+03, 1.01858e+03, 4.02875e+03, 5.62676e+02]
gaussMeanL2 = [2.32282e+02, 2.99698e+02, 5.27794e+02, 5.73367e+02, 6.55365e+02]
gaussSigmaL2 = [5.87611e+00,  8.47660e+00, 1.09542e+01, 9.21392e+00, 1.22691e+01]
gaussMeanErrL2 = [9.79582e-03, 3.45061e-02, 8.49401e-02, 3.32470e-02, 1.08469e-01]
gaussSigmaErrL2 = [8.86909e-03, 3.96406e-02, 1.03170e-01, 3.25550e-02, 1.02654e-01]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L2
# ----------------------------------------------------------------------------------------
mL2 = [3854.64, 996.999, 1402.18, 266.758, 109.242]
cL2 = [-8.17376, -0.126573, -1.82821, 1.06244, -0.00903326]
mErrL2 = [7.74665, 2.97892, 2.39432, 1.51653, 0.646186]
cErrL2 = [0.0303118, 0.00882619, 0.00408295, 0.00462727, 0.00126455]

# ----------------------------------------------------------------------------------------
# GAUSS FIT DATA - L3
# ----------------------------------------------------------------------------------------
gaussAmpL3 = [4.05747e+04, 5.76246e+03, 1.49775e+03, 6.28109e+03, 4.70674e+02]
gaussMeanL3 = [2.48857e+02, 3.22064e+02, 5.67484e+02, 6.17619e+02, 7.07714e+02]
gaussSigmaL3 = [6.46421e+00, 9.61005e+00, 1.13187e+01, 1.00259e+01, 1.61900e+01]
gaussMeanErrL3 = [8.44852e-03, 2.93532e-02, 7.08496e-02, 3.02591e-02, 1.33809e-01]
gaussSigmaErrL3 = [7.63940e-03, 3.30563e-02, 8.69842e-02, 2.90052e-02, 1.39903e-01]

# ----------------------------------------------------------------------------------------
# BACKGROUND FIT DATA - L3
# ----------------------------------------------------------------------------------------
mL3 = [5474.2, 1670.07, 688.358, -12221, 1914.16]
cL3 = [-10.9778, -0.742393, -0.317576,22.2056,  -2.42654]
mErrL3 = [8.66758, 3.68925, 1.70245, 107.005, 44.4323]
cErrL3 = [0.030301, 0.00998107, 0.00320637, 0.175851, 0.0615133]

# ----------------------------------------------------------------------------------------
# READ IN ROOT HISTOGRAMS FOR BACKGROUND SIGMA CALCULATION
# ----------------------------------------------------------------------------------------
def getHistos(f):
    hists = []
    for i in range(0,4):
        hists.append(f.Get("Uncalib_Slow_Energy_ "+str(i)))
    return hists

def getBackgroundErrFitParams(hists):
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
        
        for j in range(0,5):
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
            h1 = hists[i].Clone()
            h1.SetLineColor(root.kBlack)
            g1.SetLineColor(root.kBlue)
            pol1.SetLineColor(root.kRed)
            total.SetLineColor(root.kGreen)

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
            h1.Fit(pol1, "R+")

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
            c.SetLogy()
            g1.Draw("same")
            pol1.Draw("same")
            total.Draw("same")

            xm1=0.6
            ym1=0.6
            xm2=0.8
            ym2=0.88
            legend1 = root.TLegend(xm1,ym1,xm2,ym2)
            legend1.SetHeader("R8 - " + str(Energy[j]) + " keV")
            legend1.AddEntry(g1,"Gauss Peak","l")
            legend1.AddEntry(pol1,"Polynomial Background","l")
            legend1.AddEntry(total,"Total Fit","l")
            legend1.SetBorderSize(0)
            legend1.SetTextSize(0.03)
            legend1.Draw()
            c.Update()
            c.SaveAs(plotsDirectory+"/gauss_fits/L" + str(i) + "_R8_" + str(Energy[j]) + ".png")
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
    p1ParamsL2 = []
    p2ParamsL2 = []
    p1ParamsL1 = []
    p2ParamsL1 = []
    p1ParamsL3 = []
    p2ParamsL3 = []
    # for i in range 0,1 and 2,4
    for i in range(0,4):
        mean = eval("gaussMeanL" + str(i))
        sigmaerr = eval("gaussSigmaErrL" + str(i))
        
        plt.figure(figsize=(30,20))
        plot1 = plt.errorbar(mean, Energy, xerr=sigmaerr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
        popt, pcov = curve_fit(func, mean, Energy, sigma=sigmaerr) 
        if i == 0:
            p1ParamsL0.append(popt)
        elif i == 1:
            p1ParamsL1.append(popt)
        elif i == 2:
            p1ParamsL2.append(popt)
        elif i == 3:
            p1ParamsL3.append(popt)
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
        elif i == 2:
            p2ParamsL2.append(popt)
        elif i == 3:
            p2ParamsL3.append(popt)
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

    return  p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2


def calculateResiduals(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2):
        for i in range(0,4):
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
            plt.ylim(-2,3)
            # increase the number of ticks on the y axis
            plt.yticks(np.arange(-2, 3, 0.5))
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
            plt.ylim(-4,4)
            plt.yticks(np.arange(-4, 4, 0.5))
            plt.savefig(plotsDirectory+"/energy_calibration/L" + str(i) + "_2ndOrderResiduals.png")
            plt.close()


        
def energyResolution(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2):
        for i in range(0,4):
            mean = eval("gaussMeanL" + str(i))
            meanerr = eval("gaussMeanErrL" + str(i))
            p2ParamsL = eval("p2ParamsL" + str(i))
           #print("Energy for L" + str(i) + " is " + str(energy))
            sigma = eval("gaussSigmaL" + str(i))
            sigmaerr = eval("gaussSigmaErrL" + str(i))
            energy = func2(mean, *p2ParamsL[0])
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
            popt, pcov = curve_fit(func2, energy, energyResolution, sigma=energyResolutionErr, absolute_sigma=False)
            plt.plot(XRange, func2(XRange, *popt), 'r-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f' % tuple(popt))
            plt.legend(loc='best')
            energyResolution662 = func2(662, *popt)
            plt.text(662, energyResolution662, "Energy Resolution at 662 keV = " + str(round(energyResolution662, 2)) + "%", fontsize=32)
            plt.savefig(plotsDirectory+"/energy_resolution/L" + str(i) + "_energyResolution_2ndOrder.png")
            plt.close()

            # # do the same again but with func3
            # plt.figure(figsize=(30,20))
            # plt.errorbar(energy, energyResolution, yerr=energyResolutionErr, fmt='o', color='black', ecolor='red', elinewidth=3, capsize=1, capthick=1, markersize=13)
            # plt.grid( linestyle='-', linewidth='0.5', color='gray')
            # plt.title("Energy Resolution for L" + str(i), fontsize=39)
            # plt.ylabel("Energy Resolution [%]", fontsize=32)
            # plt.xlabel("Energy [keV]", fontsize=32)
            # plt.xticks(fontsize=32)
            # plt.yticks(fontsize=32)
            # popt, pcov = curve_fit(func3, energy, energyResolution, sigma=energyResolutionErr, absolute_sigma=False)
            # plt.plot(XRange, func3(XRange, *popt), 'r-', label='fit: a=%10.8f, b=%10.8f, c=%10.8f, d=%10.8f' % tuple(popt))
            # plt.legend(loc='best')
            # energyResolution662 = func3(662, *popt)
            # plt.text(662, energyResolution662, "Energy Resolution at 662 keV = " + str(round(energyResolution662, 2)) + "%", fontsize=32)
            # plt.savefig(plotsDirectory+"/energy_resolution/L" + str(i) + "_energyResolution_3rdOrder.png")
            # plt.close()


            
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN CODE
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():

    hists = getHistos(f)

    #getBackgroundErrFitParams(hists)

    p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2 = calibrationPlots()

    calculateResiduals(p1ParamsL0, p2ParamsL0, p1ParamsL1, p2ParamsL1, p1ParamsL3, p2ParamsL3, p1ParamsL2, p2ParamsL2)

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