################################################################################
#                               to run:                                        #
# python3 Calibration.py /directory/to/rootfile(R**_rawData.root)/ Target_name #
################################################################################

from venv import create
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math
import pandas as pd
import uproot as up
import csv as csv

filename = sys.argv[1] #open a root file where the filename is given as an argument
target = sys.argv[2]  #create a string called target to store the name of the target as an argument
f = root.TFile(filename, "UPDATE")

s0none = f.Get("L0_Uncalibrated")
s1none = f.Get("L1_Uncalibrated")
s2none = f.Get("L2_Uncalibrated")
s3none = f.Get("L3_Uncalibrated")

timeSlowSig0 =f.Get("L0_TimeSlow")
timeSlowSig1 =f.Get("L1_TimeSlow")
timeSlowSig2 =f.Get("L2_TimeSlow")
timeSlowSig3 =f.Get("L3_TimeSlow")

######################### Second Order Polynomial Fit #########################
#L0
a20 = 0.00050627
b20 = -0.041062
c20 = 190.9778 
#L1
a21 = 0.00058954
b21 = -0.044864
c21 = 190.833 
#L2
a22 = 0.00024151
b22 = -0.033294
c22 = 192.6998 
#L3
a23 = 0.00045733
b23 = -0.03787
c23 = 190.7451 
##############################################################################

s0calib2 = root.TH1D("L0_Calibrated","LaBr3 Det 0 Slow Signal Energy Spectrum - Calibrated 2nd Order", 6000, 0, 6000) #create a new TH1D histogram called s0calib2 with the same binning as s0none
s1calib2 = root.TH1D("L1_Calibrated","LaBr3 Det 1 Slow Signal Energy Spectrum - Calibrated 2nd Order", 6000, 0, 6000) #create a new TH1D histogram called s1calib2 with the same binning as s1none
s2calib2 = root.TH1D("L2_Calibrated","LaBr3 Det 2 Slow Signal Energy Spectrum - Calibrated 2nd Order", 6000, 0, 6000) #create a new TH1D histogram called s2calib2 with the same binning as s2none
s3calib2 = root.TH1D("L3_Calibrated","LaBr3 Det 3 Slow Signal Energy Spectrum - Calibrated 2nd Order", 6000, 0, 6000) #create a new TH1D histogram called s3calib2 with the same binning as s3none

energyVtime0 = root.TH2D("L0_EnergyVTime","LaBr3 Det 0 Slow Signal Energy vs Time",6000,0,1e15,6000,0,6000);
energyVtime1 = root.TH2D("L1_EnergyVTime","LaBr3 Det 1 Slow Signal Energy vs Time",6000,0,1e15,6000,0,6000);
energyVtime2 = root.TH2D("L2_EnergyVTime","LaBr3 Det 2 Slow Signal Energy vs Time",6000,0,1e15,6000,0,6000);
energyVtime3 = root.TH2D("L3_EnergyVTime","LaBr3 Det 3 Slow Signal Energy vs Time",6000,0,1e15,6000,0,6000);


r = root.TRandom3()# create a TRandom3 object to generate random numbers

##############################################################################

N0 = s0none.GetXaxis().GetNbins()# get the number of bins in the histogram
binWidth0 = s0none.GetBinWidth(1)# get the bin width of the histogram
N1 = s1none.GetXaxis().GetNbins()# get the number of bins in the histogram
binWidth1 = s1none.GetBinWidth(1)# get the bin width of the histogram
N2 = s2none.GetXaxis().GetNbins()# get the number of bins in the histogram
binWidth2 = s2none.GetBinWidth(1)# get the bin width of the histogram
N3 = s3none.GetXaxis().GetNbins()# get the number of bins in the histogram
binWidth3 = s3none.GetBinWidth(1)# get the bin width of the histogram

##############################################################################

for i in range(N0):# loop over the bins of the histogram
    n0 = s0none.GetBinContent(i)# get the content of the bin
    n0 = int(n0) #make n an integer 
    ch0 = s0none.GetBinCenter(i)# get the center of the bin

    for j in range(n0):
        rndUnif=r.Rndm()*binWidth0-(binWidth0/2)# generate a random number between -binWidth/2 and +binWidth/2
        rndCh0 = ch0+rndUnif# add the random number to the center of the bin
        enerCalib0 = a20*rndCh0**2+b20*rndCh0+c20# calculate the calibrated energy of the bin
        s0calib2.Fill(enerCalib0)# fill the s0calib2 histogram with the calibrated energy of the bin
print("... calibrated L0") 

##############################################################################

for i in range(N1): 
    n1 = s1none.GetBinContent(i)
    n1 = int(n1) #make n an integer
    ch1 = s1none.GetBinCenter(i) # get the center of the bin

    for j in range(n1):
        rndUnif=r.Rndm()*binWidth1-(binWidth1/2)
        rndCh1 = ch1+rndUnif
        enerCalib1 = a21*ch1**2+b21*ch1+c21 
        s1calib2.Fill(enerCalib1)
print("... calibrated L1")

##############################################################################

for i in range(N2):
    n2 = s2none.GetBinContent(i)
    n2 = int(n2)
    ch2 = s2none.GetBinCenter(i)

    for j in range(n2):
        rndUnif=r.Rndm()*binWidth2-(binWidth2/2)
        rndCh2 = ch2+rndUnif
        enerCalib2 = a22*ch2**2+b22*ch2+c22
        s2calib2.Fill(enerCalib2)
print("... calibrated L2")

##############################################################################

for i in range(N3):
    n3 = s3none.GetBinContent(i)
    n3 = int(n3)
    ch3 = s3none.GetBinCenter(i)

    for j in range(n3):
        rndUnif=r.Rndm()*binWidth3-(binWidth3/2)
        rndCh3 = ch3+rndUnif
        enerCalib3 = a23*rndCh3**2+b23*rndCh3+c23
        s3calib2.Fill(enerCalib3)
print("... calibrated L3")

##############################################################################

f.WriteObject(s0calib2, "L0_Calibrated")# write the s0calib2 histogram to the root file
f.WriteObject(s1calib2, "L1_Calibrated")# write the s1calib2 histogram to the root file
f.WriteObject(s2calib2, "L2_Calibrated")# write the s2calib2 histogram to the root file
f.WriteObject(s3calib2, "L3_Calibrated")# write the s3calib2 histogram to the root file

print("Writing calibrated "+filename+"...")

##############################################################################

fileName = filename[filename.rfind("/")+1:]
fileName = fileName[:fileName.rfind("_rawData.root")]

saveDir = filename[:filename.rfind("/")]

# remove the last 16 characters from the filename
path = fileName[:-16]

# use os.path to join the saveDir and path
saveDir = os.path.join(saveDir, path)

df = pd.DataFrame(list())
df.to_csv(saveDir+"L0_Calibrated_"+fileName+"_"+target+".csv", index=False, header=False)
df.to_csv(saveDir+"L1_Calibrated_"+fileName+"_"+target+".csv", index=False, header=False)
df.to_csv(saveDir+"L2_Calibrated_"+fileName+"_"+target+".csv", index=False, header=False)
df.to_csv(saveDir+"L3_Calibrated_"+fileName+"_"+target+".csv", index=False, header=False)
 
# create a csv file  of the time and counts from the timeSlowSig0 histogram and the energy and counts from the s0calib2 histogram
with open(saveDir+"L0_Calibrated_"+fileName+"_"+target+".csv", "w") as g:       
    writer = csv.writer(g)
    writer.writerow(["Calibrated Energy [keV]", "Energy Counts", "TimeStamp Slow Signal [s]", "Time Counts"])
    for i in range(s0calib2.GetXaxis().GetNbins()):
        writer.writerow([s0calib2.GetXaxis().GetBinCenter(i), s0calib2.GetBinContent(i), timeSlowSig0.GetXaxis().GetBinCenter(i), timeSlowSig0.GetBinContent(i)])
with open(saveDir+"L1_Calibrated_"+fileName+"_"+target+".csv", "w") as g:
    writer = csv.writer(g)
    writer.writerow(["Calibrated Energy [keV]", "Energy Counts", "TimeStamp Slow Signal [s]", "Time Counts"])
    for i in range(s1calib2.GetXaxis().GetNbins()):
        writer.writerow([s1calib2.GetXaxis().GetBinCenter(i), s1calib2.GetBinContent(i), timeSlowSig1.GetXaxis().GetBinCenter(i), timeSlowSig1.GetBinContent(i)])
with open(saveDir+"L2_Calibrated_"+fileName+"_"+target+".csv", "w") as g:
    writer = csv.writer(g)
    writer.writerow(["Calibrated Energy [keV]", "Energy Counts", "TimeStamp Slow Signal [s]", "Time Counts"])
    for i in range(s2calib2.GetXaxis().GetNbins()):
        writer.writerow([s2calib2.GetXaxis().GetBinCenter(i), s2calib2.GetBinContent(i), timeSlowSig2.GetXaxis().GetBinCenter(i), timeSlowSig2.GetBinContent(i)])
with open(saveDir+"L3_Calibrated_"+fileName+"_"+target+".csv", "w") as g:
    writer = csv.writer(g)
    writer.writerow(["Calibrated Energy [keV]", "Energy Counts", "TimeStamp Slow Signal [s]", "Time Counts"])
    for i in range(s3calib2.GetXaxis().GetNbins()):
        writer.writerow([s3calib2.GetXaxis().GetBinCenter(i), s3calib2.GetBinContent(i), timeSlowSig3.GetXaxis().GetBinCenter(i), timeSlowSig3.GetBinContent(i)])

print("Writing CSV files to "+saveDir+"...")



##############################################################################

f.Close()# close the root file


