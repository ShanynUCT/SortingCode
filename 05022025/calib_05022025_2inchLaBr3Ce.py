import csv
import pandas as pd
import sys
import ROOT as ROOT
from pathlib import Path
from matplotlib import pyplot as plt
import random
import numpy as np
import os

########################
#HOW TO RUN: ~/miniconda/envs/my_root_env/bin/python calib_05020205_2inchLaBr3Ce.py /Users/shanyn/Documents/PhD/Exp/2025/05022025/RXX.root
# ______________________________________________________________________________
# Read in data
dir_path = sys.argv[1]
dir_path = dir_path[:-len(dir_path.split('/')[-1])]
calib_dir = Path(dir_path + 'calibrated/')
residuals_dir = Path(dir_path + 'residuals/')

if calib_dir.exists():
    print('\nThe directory {} already exists.'.format(calib_dir))
else:
    calib_dir.mkdir(parents=True, exist_ok=False)
    print('\nThe directory {} was created.'.format(calib_dir))

if residuals_dir.exists():
    print('\nThe directory {} already exists.'.format(residuals_dir))
else:
    residuals_dir.mkdir(parents=True, exist_ok=False)
    print('\nThe directory {} was created.'.format(residuals_dir))

file_path = Path(sys.argv[1])
name = file_path.stem
rnd = random.Random()  
r = ROOT.TRandom3(1)

detectorname = file_path.parts[-4]
detector = file_path.parts[-3]
date = file_path.parts[-2]
filename = file_path.parts[-1]
Rvalue = file_path.parts[-1].split('.')[0]

# print the folder name
print(' ')
print(' ')
print('Detector: {}'.format(detectorname) + ' ' + detector)
print('Date: {}'.format(date))
print('Filename: {}'.format(filename))
print(' ')
print(' ')

dfslowEL0 = []
dffastEL0 = []
dfslowEL1 = []
dffastEL1 = []

file = dir_path + filename
f = ROOT.TFile(file)
t = f.Get("LaBrData")
for i in range(t.GetEntries()):
    t.GetEntry(i)
    dfslowEL0.append(t.slowEL0)
    dffastEL0.append(t.fastEL0)
    dfslowEL1.append(t.slowEL1)
    dffastEL1.append(t.fastEL1)
f.Close()


# create a dataframe with the data
df = pd.DataFrame({'slowEL0': dfslowEL0, 'fastEL0': dffastEL0, 'slowEL1': dfslowEL1, 'fastEL1': dffastEL1})
df = df.dropna()
df = df[(df['slowEL0'] > 0) | (df['fastEL0'] > 0) | (df['slowEL1'] > 0) | (df['fastEL1'] > 0)]
df = df[(df['slowEL0'] < 16350) | (df['fastEL0'] < 16350) | (df['slowEL1'] < 16350) | (df['fastEL1'] < 16350)]
df = df.reset_index(drop=True)

# ______________________________________________________________________________ plot ROOT histogram for channels
h_fastEL0 = ROOT.TH1D("h_fastEL0", "h_fastEL0", 16350, 0, 16350)
for i in range(len(df)):
    h_fastEL0.Fill(df['fastEL0'][i])
h_fastEL0.Draw()
h_fastEL0.GetXaxis().SetTitle("Channel")
h_fastEL0.GetYaxis().SetTitle("Counts (a.u.)")
h_fastEL0.SetLineColor(1) # red
h_fastEL0.SetLineWidth(1)
h_fastEL0.SetStats(0)
h_fastEL0.SaveAs(str(dir_path + 'h_fastEL0.root'))

h_slowEL0 = ROOT.TH1D("h_slowEL0", "h_slowEL0", 16350, 0, 16350)
for i in range(len(df)):
    h_slowEL0.Fill(df['slowEL0'][i])
h_slowEL0.Draw()
h_slowEL0.GetXaxis().SetTitle("Channel")
h_slowEL0.GetYaxis().SetTitle("Counts (a.u.)")
h_slowEL0.SetLineColor(1) # red
h_slowEL0.SetLineWidth(1)
h_slowEL0.SetStats(0)
h_slowEL0.SaveAs(str(dir_path + 'h_slowEL0.root'))

h_fastEL1 = ROOT.TH1D("h_fastEL1", "h_fastEL1", 16350, 0, 16350)
for i in range(len(df)):
    h_fastEL1.Fill(df['fastEL1'][i])
h_fastEL1.Draw()
h_fastEL1.GetXaxis().SetTitle("Channel")
h_fastEL1.GetYaxis().SetTitle("Counts (a.u.)")
h_fastEL1.SetLineColor(1) # red
h_fastEL1.SetLineWidth(1)
h_fastEL1.SetStats(0)
h_fastEL1.SaveAs(str(dir_path + 'h_fastEL1.root'))

h_slowEL1 = ROOT.TH1D("h_slowEL1", "h_slowEL1", 16350, 0, 16350)
for i in range(len(df)):
    h_slowEL1.Fill(df['slowEL1'][i])
h_slowEL1.Draw()
h_slowEL1.GetXaxis().SetTitle("Channel")
h_slowEL1.GetYaxis().SetTitle("Counts (a.u.)")
h_slowEL1.SetLineColor(1) # red
h_slowEL1.SetLineWidth(1)
h_slowEL1.SetStats(0)
h_slowEL1.SaveAs(str(dir_path + 'h_slowEL1.root'))

# ______________________________________________________________________________ peaks
    
""" energy = [0, 511.0, 1224.5, 1173.2, 1332.5]
ch_slowEL0 = [0,1018.89, 2523.68, 2318.18, 2630.58]
ch_fastEL0 = [0, 1608.96, 3997.23, 3665.27, 4161.23]
ch_slowEL1 = [0, 1164.77, 2886.66, 2643.18, 2998.88]
ch_fastEL1 = [0, 1714.76, 4260.72, 3892.81, 4409.04] """

energy = [0, 1173.2, 1332.5]

# for runs <= R5 (2inch-2inch LaBr3Ce coincidence runs)
""" ch_slowEL0 = [0, 2318.18, 2630.58]
ch_fastEL0 = [0,3665.27, 4161.23]
ch_slowEL1 = [0, 2643.18, 2998.88]
ch_fastEL1 = [0, 3892.81, 4409.04]  """

# for runs <= R5 (2inch-DA LaBr3Ce coincidence runs)
ch_slowEL0 = [0, 10596.9, 11568.8]
ch_fastEL0 = [0, 9016.79, 9940.73]
ch_slowEL1 = [0, 2643.18, 2998.88]
ch_fastEL1 = [0, 3892.81, 4409.04] 


scaleslowEL0 = ch_slowEL0[1]/energy[1]
scalefastEL0 = ch_fastEL0[1]/energy[1]
scaleslowEL1 = ch_slowEL1[1]/energy[1]
scalefastEL1 = ch_fastEL1[1]/energy[1]

# ______________________________________________________________________________ residuals ch_slowEL0
residuals1_slowEL0 = []
x = ROOT.TCanvas("c", "c", 800, 600)
x.SetGrid()
pol1 = ROOT.TF1("pol1", "pol1", 0, 2000)
#plot energy vs ch_slowEL0. Fit the data with a first order polynomial, and plot the residuals (energy - fit)
h = ROOT.TH2D("h", "h", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h.Fill(ch_slowEL0[i], energy[i])
h.Draw()
h.GetXaxis().SetTitle("Channel (a.u.)")
h.GetYaxis().SetTitle("Energy (keV)")
h.Fit("pol1")
h.GetFunction("pol1").SetLineColor(2)
h.GetFunction("pol1").SetLineWidth(1)
h.GetFunction("pol1").Draw("same")
h.SetStats(0)
h.SetMarkerStyle(20)
h.SetMarkerSize(2)
h.SetMarkerColor(1)
x.Update()
x.SaveAs(str(residuals_dir / (name + '_pol1_slowEL0_fit.root')))
x.SaveAs(str(residuals_dir / (name + '_pol1_slowEL0_fit.png')))
x.Close()

pol1_params_slowEL0 = []
for i in range(2):
    pol1_params_slowEL0.append(h.GetFunction("pol1").GetParameter(i))
for i in range(len(energy)):
    residuals1_slowEL0.append([energy[i] - (pol1_params_slowEL0[0] + pol1_params_slowEL0[1]*ch_slowEL0[i])])
with open(str(residuals_dir / (name + '_pol1_slowEL0_residuals.txt')), 'w') as f:
    for item in residuals1_slowEL0:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol1_slowEL0_residuals.txt')))

x1 = ROOT.TCanvas("c1", "c1", 800, 600)
x1.SetGrid()
h1 = ROOT.TH2D("h1", "h1", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h1.Fill(energy[i], residuals1_slowEL0[i][0])
h1.Draw()
h1.GetXaxis().SetTitle("Energy (keV)")
h1.GetYaxis().SetTitle("Residuals (keV)")
h1.SetStats(0)
h1.SetMarkerStyle(20)
h1.SetMarkerSize(2)
x1.Update()
x1.SaveAs(str(residuals_dir / (name + '_pol1_slowEL0_residuals.root')))
x1.SaveAs(str(residuals_dir / (name + '_pol1_slowEL0_residuals.png')))
x1.Close()

residuals2_slowEL0 = []
x2 = ROOT.TCanvas("c2", "c2", 800, 600)
x2.SetGrid()
pol2 = ROOT.TF1("pol2", "pol2", 0, 2000)
#plot energy vs ch_slowEL0. Fit the data with a second order polynomial, and plot the residuals (energy - fit)
h2 = ROOT.TH2D("h2", "h2", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h2.Fill(ch_slowEL0[i], energy[i])
h2.Draw()
h2.GetXaxis().SetTitle("Channel (a.u.)")
h2.GetYaxis().SetTitle("Energy (keV)")
h2.Fit("pol2")
h2.GetFunction("pol2").SetLineColor(2)
h2.GetFunction("pol2").SetLineWidth(1)
h2.GetFunction("pol2").Draw("same")
h2.SetStats(0)
h2.SetMarkerStyle(20)
h2.SetMarkerSize(2)
h2.SetMarkerColor(1)
x2.Update()
x2.SaveAs(str(residuals_dir / (name + '_pol2_slowEL0_fit.root')))
x2.SaveAs(str(residuals_dir / (name + '_pol2_slowEL0_fit.png')))
x2.Close()

pol2_params_slowEL0 = []
for i in range(3):
    pol2_params_slowEL0.append(h2.GetFunction("pol2").GetParameter(i))
for i in range(len(energy)):
    residuals2_slowEL0.append([energy[i] - (pol2_params_slowEL0[0] + pol2_params_slowEL0[1]*ch_slowEL0[i] + pol2_params_slowEL0[2]*ch_slowEL0[i]**2)])
with open(str(residuals_dir / (name + '_pol2_slowEL0_residuals.txt')), 'w') as f:
    for item in residuals2_slowEL0:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol2_slowEL0_residuals.txt')))

x3 = ROOT.TCanvas("c3", "c3", 800, 600)
x3.SetGrid()
h3 = ROOT.TH2D("h3", "h3", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h3.Fill(energy[i], residuals2_slowEL0[i][0])
h3.Draw()
h3.GetXaxis().SetTitle("Energy (keV)")
h3.GetYaxis().SetTitle("Residuals (keV)")
h3.SetStats(0)
h3.SetMarkerStyle(20)
h3.SetMarkerSize(2)
x3.Update()
x3.SaveAs(str(residuals_dir / (name + '_pol2_slowEL0_residuals.root')))
x3.SaveAs(str(residuals_dir / (name + '_pol2_slowEL0_residuals.png')))
x3.Close()

residuals3_slowEL0 = []
x4 = ROOT.TCanvas("c4", "c4", 800, 600)
x4.SetGrid()
pol3 = ROOT.TF1("pol3", "pol3", 0, 2000)
#plot energy vs ch_slowEL0. Fit the data with a third order polynomial, and plot the residuals (energy - fit)
h4 = ROOT.TH2D("h4", "h4", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h4.Fill(ch_slowEL0[i], energy[i])
h4.Draw()
h4.GetXaxis().SetTitle("Channel (a.u.)")
h4.GetYaxis().SetTitle("Energy (keV)")
h4.Fit("pol3")
h4.GetFunction("pol3").SetLineColor(2)
h4.GetFunction("pol3").SetLineWidth(1)
h4.GetFunction("pol3").Draw("same")
h4.SetStats(0)
h4.SetMarkerStyle(20)
h4.SetMarkerSize(2)
h4.SetMarkerColor(1)
x4.Update()
x4.SaveAs(str(residuals_dir / (name + '_pol3_slowEL0_fit.root')))
x4.SaveAs(str(residuals_dir / (name + '_pol3_slowEL0_fit.png')))
x4.Close()

pol3_params_slowEL0 = []
for i in range(4):
    pol3_params_slowEL0.append(h4.GetFunction("pol3").GetParameter(i))
for i in range(len(energy)):
    residuals3_slowEL0.append([energy[i] - (pol3_params_slowEL0[0] + pol3_params_slowEL0[1]*ch_slowEL0[i] + pol3_params_slowEL0[2]*ch_slowEL0[i]**2 + pol3_params_slowEL0[3]*ch_slowEL0[i]**3)])
with open(str(residuals_dir / (name + '_pol3_slowEL0_residuals.txt')), 'w') as f:
    for item in residuals3_slowEL0:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol3_slowEL0_residuals.txt')))

x5 = ROOT.TCanvas("c5", "c5", 800, 600)
x5.SetGrid()
h5 = ROOT.TH2D("h5", "h5", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h5.Fill(energy[i], residuals3_slowEL0[i][0])
h5.Draw()
h5.GetXaxis().SetTitle("Energy (keV)")
h5.GetYaxis().SetTitle("Residuals (keV)")
h5.SetStats(0)
h5.SetMarkerStyle(20)
h5.SetMarkerSize(2)
x5.Update()
x5.SaveAs(str(residuals_dir / (name + '_pol3_slowEL0_residuals.root')))
x5.SaveAs(str(residuals_dir / (name + '_pol3_slowEL0_residuals.png')))
x5.Close()

# ______________________________________________________________________________ best fit
# Get chi-square values for each fit
chi2_1 = h.GetFunction("pol1").GetChisquare()
chi2_2 = h2.GetFunction("pol2").GetChisquare()
chi2_3 = h4.GetFunction("pol3").GetChisquare()

# Handle zero chi-square values
if chi2_1 == 0:
    residuals1_slowEL0 = [np.nan]
if chi2_2 == 0:
    residuals2_slowEL0 = [np.nan]
if chi2_3 == 0:
    residuals3_slowEL0 = [np.nan]

# Determine which polynomial fit has the smallest residuals
fits_slowEL0 = [(residuals1_slowEL0, "pol1"), (residuals2_slowEL0, "pol2"), (residuals3_slowEL0, "pol3")]
min_fit_slowEL0 = min(fits_slowEL0, key=lambda x: sum(abs(val) if isinstance(val, (int, float)) else sum(abs(inner_val) for inner_val in val) for val in x[0]))


# Assign the best fit function name
bestfit = min_fit_slowEL0[1]
if min_fit_slowEL0[1] == 'pol1':
    p0_slowEL0 = pol1_params_slowEL0[0]
    p1_slowEL0 = pol1_params_slowEL0[1]
    p2_slowEL0 = 0
    p3_slowEL0 = 0
elif min_fit_slowEL0[1] == 'pol2':
    p0_slowEL0 = pol2_params_slowEL0[0]
    p1_slowEL0 = pol2_params_slowEL0[1]
    p2_slowEL0 = pol2_params_slowEL0[2]
    p3_slowEL0 = 0
elif min_fit_slowEL0[1] == 'pol3':
    p0_slowEL0 = pol3_params_slowEL0[0]
    p1_slowEL0 = pol3_params_slowEL0[1]
    p2_slowEL0 = pol3_params_slowEL0[2]
    p3_slowEL0 = pol3_params_slowEL0[3]

# ______________________________________________________________________________ residuals ch_fastEL0
residuals1_fastEL0 = []
x = ROOT.TCanvas("c", "c", 800, 600)
x.SetGrid()
pol1 = ROOT.TF1("pol1", "pol1", 0, 2000)
pol2 = ROOT.TF1("pol2", "pol2", 0, 2000)
pol3 = ROOT.TF1("pol3", "pol3", 0, 2000)

pol1.SetParameter(0, 0)
pol2.SetParameter(0, 0)
pol3.SetParameter(0, 0)

h = ROOT.TH2D("h", "h", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h.Fill(ch_fastEL0[i], energy[i])
h.Draw()
h.GetXaxis().SetTitle("Channel (a.u.)")
h.GetYaxis().SetTitle("Energy (keV)")
h.Fit("pol1")
h.GetFunction("pol1").SetLineColor(2)
h.GetFunction("pol1").SetLineWidth(1)
h.GetFunction("pol1").Draw("same")
h.SetStats(0)
h.SetMarkerStyle(20)
h.SetMarkerSize(2)
h.SetMarkerColor(1)
x.Update()
x.SaveAs(str(residuals_dir / (name + '_pol1_fastEL0_fit.root')))
x.SaveAs(str(residuals_dir / (name + '_pol1_fastEL0_fit.png')))
x.Close()

pol1_params_fastEL0 = []
for i in range(2):
    pol1_params_fastEL0.append(h.GetFunction("pol1").GetParameter(i))
for i in range(len(energy)):
    residuals1_fastEL0.append([energy[i] - (pol1_params_fastEL0[0] + pol1_params_fastEL0[1]*ch_fastEL0[i])])
with open(str(residuals_dir / (name + '_pol1_fastEL0_residuals.txt')), 'w') as f:
    for item in residuals1_fastEL0:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol1_fastEL0_residuals.txt')))

x1 = ROOT.TCanvas("c1", "c1", 800, 600)
x1.SetGrid()
h1 = ROOT.TH2D("h1", "h1", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h1.Fill(energy[i], residuals1_fastEL0[i][0])
h1.Draw()
h1.GetXaxis().SetTitle("Energy (keV)")
h1.GetYaxis().SetTitle("Residuals (keV)")
h1.SetStats(0)
h1.SetMarkerStyle(20)
h1.SetMarkerSize(2)
x1.Update()
x1.SaveAs(str(residuals_dir / (name + '_pol1_fastEL0_residuals.root')))
x1.SaveAs(str(residuals_dir / (name + '_pol1_fastEL0_residuals.png')))
x1.Close()

residuals2_fastEL0 = []
x2 = ROOT.TCanvas("c2", "c2", 800, 600)
x2.SetGrid()
h2 = ROOT.TH2D("h2", "h2", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h2.Fill(ch_fastEL0[i], energy[i])
h2.Draw()
h2.GetXaxis().SetTitle("Channel (a.u.)")
h2.GetYaxis().SetTitle("Energy (keV)")
h2.Fit("pol2")
h2.GetFunction("pol2").SetLineColor(2)
h2.GetFunction("pol2").SetLineWidth(1)
h2.GetFunction("pol2").Draw("same")
h2.SetStats(0)
h2.SetMarkerStyle(20)
h2.SetMarkerSize(2)
h2.SetMarkerColor(1)
x2.Update()
x2.SaveAs(str(residuals_dir / (name + '_pol2_fastEL0_fit.root')))
x2.SaveAs(str(residuals_dir / (name + '_pol2_fastEL0_fit.png')))
x2.Close()

pol2_params_fastEL0 = []
for i in range(3):
    pol2_params_fastEL0.append(h2.GetFunction("pol2").GetParameter(i))
for i in range(len(energy)):
    residuals2_fastEL0.append([energy[i] - (pol2_params_fastEL0[0] + pol2_params_fastEL0[1]*ch_fastEL0[i] + pol2_params_fastEL0[2]*ch_fastEL0[i]**2)])
with open(str(residuals_dir / (name + '_pol2_fastEL0_residuals.txt')), 'w') as f:
    for item in residuals2_fastEL0:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol2_fastEL0_residuals.txt')))

x3 = ROOT.TCanvas("c3", "c3", 800, 600)
x3.SetGrid()
h3 = ROOT.TH2D("h3", "h3", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h3.Fill(energy[i], residuals2_fastEL0[i][0])
h3.Draw()
h3.GetXaxis().SetTitle("Energy (keV)")
h3.GetYaxis().SetTitle("Residuals (keV)")
h3.SetStats(0)
h3.SetMarkerStyle(20)
h3.SetMarkerSize(2)
x3.Update()
x3.SaveAs(str(residuals_dir / (name + '_pol2_fastEL0_residuals.root')))
x3.SaveAs(str(residuals_dir / (name + '_pol2_fastEL0_residuals.png')))
x3.Close()

residuals3_fastEL0 = []
x4 = ROOT.TCanvas("c4", "c4", 800, 600)
x4.SetGrid()
#plot energy vs ch_fastEL0. Fit the data with a third order polynomial, and plot the residuals (energy - fit)
h4 = ROOT.TH2D("h4", "h4", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h4.Fill(ch_fastEL0[i], energy[i])
h4.Draw()
h4.GetXaxis().SetTitle("Channel (a.u.)")
h4.GetYaxis().SetTitle("Energy (keV)")
h4.Fit("pol3")
h4.GetFunction("pol3").SetLineColor(2)
h4.GetFunction("pol3").SetLineWidth(1)
h4.GetFunction("pol3").Draw("same")
h4.SetStats(0)
h4.SetMarkerStyle(20)
h4.SetMarkerSize(2)
h4.SetMarkerColor(1)
x4.Update()
x4.SaveAs(str(residuals_dir / (name + '_pol3_fastEL0_fit.root')))
x4.SaveAs(str(residuals_dir / ( name + '_pol3_fastEL0_fit.png')))
x4.Close()

pol3_params_fastEL0 = []
for i in range(4):
    pol3_params_fastEL0.append(h4.GetFunction("pol3").GetParameter(i))
for i in range(len(energy)):
    residuals3_fastEL0.append([energy[i] - (pol3_params_fastEL0[0] + pol3_params_fastEL0[1]*ch_fastEL0[i] + pol3_params_fastEL0[2]*ch_fastEL0[i]**2 + pol3_params_fastEL0[3]*ch_fastEL0[i]**3)])
with open(str(residuals_dir / (name + '_pol3_fastEL0_residuals.txt')), 'w') as f:
    for item in residuals3_fastEL0:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol3_fastEL0_residuals.txt')))

x5 = ROOT.TCanvas("c5", "c5", 800, 600)
x5.SetGrid()
h5 = ROOT.TH2D("h5", "h5", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h5.Fill(energy[i], residuals3_fastEL0[i][0])
h5.Draw()
h5.GetXaxis().SetTitle("Energy (keV)")
h5.GetYaxis().SetTitle("Residuals (keV)")
h5.SetStats(0)
h5.SetMarkerStyle(20)
h5.SetMarkerSize(2)
x5.Update()
x5.SaveAs(str(residuals_dir / (name + '_pol3_fastEL0_residuals.root')))
x5.SaveAs(str(residuals_dir / (name + '_pol3_fastEL0_residuals.png')))
x5.Close()

# ______________________________________________________________________________ best fit
# Get chi-square values for each fit
chi2_1 = h.GetFunction("pol1").GetChisquare()
chi2_2 = h2.GetFunction("pol2").GetChisquare()
chi2_3 = h4.GetFunction("pol3").GetChisquare()

# Handle zero chi-square values
if chi2_1 == 0:
    residuals1_fastEL0 = [np.nan]
if chi2_2 == 0:
    residuals2_fastEL0 = [np.nan]
if chi2_3 == 0:
    residuals3_fastEL0 = [np.nan]

# Determine which polynomial fit has the smallest residuals
fits_fastEL0 = [(residuals1_fastEL0, "pol1"), (residuals2_fastEL0, "pol2"), (residuals3_fastEL0, "pol3")]
min_fit_fastEL0 = min(fits_fastEL0, key=lambda x: sum(abs(val) if isinstance(val, (int, float)) else sum(abs(inner_val) for inner_val in val) for val in x[0]))


# Assign the best fit function name
bestfit = min_fit_fastEL0[1]
if min_fit_fastEL0[1] == 'pol1':
    p0_fastEL0 = pol1_params_fastEL0[0]
    p1_fastEL0 = pol1_params_fastEL0[1]
    p2_fastEL0 = 0
    p3_fastEL0 = 0
elif min_fit_fastEL0[1] == 'pol2':
    p0_fastEL0 = pol2_params_fastEL0[0]
    p1_fastEL0 = pol2_params_fastEL0[1]
    p2_fastEL0 = pol2_params_fastEL0[2]
    p3_fastEL0 = 0
elif min_fit_fastEL0[1] == 'pol3':
    p0_fastEL0 = pol3_params_fastEL0[0]
    p1_fastEL0 = pol3_params_fastEL0[1]
    p2_fastEL0 = pol3_params_fastEL0[2]
    p3_fastEL0 = pol3_params_fastEL0[3]

# ______________________________________________________________________________ residuals ch_slowEL1
residuals1_slowEL1 = []
x = ROOT.TCanvas("c", "c", 800, 600)
x.SetGrid()
pol1 = ROOT.TF1("pol1", "pol1", 0, 2000)
pol2 = ROOT.TF1("pol2", "pol2", 0, 2000)
pol3 = ROOT.TF1("pol3", "pol3", 0, 2000)

pol1.SetParameter(0, 0)
pol2.SetParameter(0, 0)
pol3.SetParameter(0, 0)

h = ROOT.TH2D("h", "h", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h.Fill(ch_slowEL1[i], energy[i])
h.Draw()
h.GetXaxis().SetTitle("Channel (a.u.)")
h.GetYaxis().SetTitle("Energy (keV)")
h.Fit("pol1")
h.GetFunction("pol1").SetLineColor(2)
h.GetFunction("pol1").SetLineWidth(1)
h.GetFunction("pol1").Draw("same")
h.SetStats(0)
h.SetMarkerStyle(20)
h.SetMarkerSize(2)
h.SetMarkerColor(1)
x.Update()
x.SaveAs(str(residuals_dir / (name + '_pol1_slowEL1_fit.root')))
x.SaveAs(str(residuals_dir / (name + '_pol1_slowEL1_fit.png')))
x.Close()

pol1_params_slowEL1 = []
for i in range(2):
    pol1_params_slowEL1.append(h.GetFunction("pol1").GetParameter(i))
for i in range(len(energy)):
    residuals1_slowEL1.append([energy[i] - (pol1_params_slowEL1[0] + pol1_params_slowEL1[1]*ch_slowEL1[i])])
with open(str(residuals_dir / (name + '_pol1_slowEL1_residuals.txt')), 'w') as f:
    for item in residuals1_slowEL1:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol1_slowEL1_residuals.txt')))

x1 = ROOT.TCanvas("c1", "c1", 800, 600)
x1.SetGrid()
h1 = ROOT.TH2D("h1", "h1", 2000, 0, 2000, 50 ,-25,25)
for i in range(len(energy)):
    h1.Fill(energy[i], residuals1_slowEL1[i][0])
h1.Draw()
h1.GetXaxis().SetTitle("Energy (keV)")
h1.GetYaxis().SetTitle("Residuals (keV)")
h1.SetStats(0)
h1.SetMarkerStyle(20)
h1.SetMarkerSize(2)
x1.Update()
x1.SaveAs(str(residuals_dir / (name + '_pol1_slowEL1_residuals.root')))
x1.SaveAs(str(residuals_dir / (name + '_pol1_slowEL1_residuals.png')))
x1.Close()

residuals2_slowEL1 = []
x2 = ROOT.TCanvas("c2", "c2", 800, 600)
x2.SetGrid()
h2 = ROOT.TH2D("h2", "h2", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h2.Fill(ch_slowEL1[i], energy[i])
h2.Draw()
h2.GetXaxis().SetTitle("Channel (a.u.)")
h2.GetYaxis().SetTitle("Energy (keV)")
h2.Fit("pol2")
h2.GetFunction("pol2").SetLineColor(2)
h2.GetFunction("pol2").SetLineWidth(1)
h2.GetFunction("pol2").Draw("same")
h2.SetStats(0)
h2.SetMarkerStyle(20)
h2.SetMarkerSize(2)
h2.SetMarkerColor(1)
x2.Update()
x2.SaveAs(str(residuals_dir / (name + '_pol2_slowEL1_fit.root')))
x2.SaveAs(str(residuals_dir / (name + '_pol2_slowEL1_fit.png')))
x2.Close()

pol2_params_slowEL1 = []
for i in range(3):
    pol2_params_slowEL1.append(h2.GetFunction("pol2").GetParameter(i))
for i in range(len(energy)):
    residuals2_slowEL1.append([energy[i] - (pol2_params_slowEL1[0] + pol2_params_slowEL1[1]*ch_slowEL1[i] + pol2_params_slowEL1[2]*ch_slowEL1[i]**2)])
with open(str(residuals_dir / (name + '_pol2_slowEL1_residuals.txt')), 'w') as f:
    for item in residuals2_slowEL1:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol2_slowEL1_residuals.txt')))

x3 = ROOT.TCanvas("c3", "c3", 800, 600)
x3.SetGrid()
h3 = ROOT.TH2D("h3", "h3", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h3.Fill(energy[i], residuals2_slowEL1[i][0])
h3.Draw()
h3.GetXaxis().SetTitle("Energy (keV)")
h3.GetYaxis().SetTitle("Residuals (keV)")
h3.SetStats(0)
h3.SetMarkerStyle(20)
h3.SetMarkerSize(2)
x3.Update()
x3.SaveAs(str(residuals_dir / (name + '_pol2_slowEL1_residuals.root')))
x3.SaveAs(str(residuals_dir / (name + '_pol2_slowEL1_residuals.png')))
x3.Close()

residuals3_slowEL1 = []
x4 = ROOT.TCanvas("c4", "c4", 800, 600)
x4.SetGrid()
h4 = ROOT.TH2D("h4", "h4", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h4.Fill(ch_slowEL1[i], energy[i])
h4.Draw()
h4.GetXaxis().SetTitle("Channel (a.u.)")
h4.GetYaxis().SetTitle("Energy (keV)")
h4.Fit("pol3")
h4.GetFunction("pol3").SetLineColor(2)
h4.GetFunction("pol3").SetLineWidth(1)
h4.GetFunction("pol3").Draw("same")
h4.SetStats(0)
h4.SetMarkerStyle(20)
h4.SetMarkerSize(2)
h4.SetMarkerColor(1)
x4.Update()
x4.SaveAs(str(residuals_dir / (name + '_pol3_slowEL1_fit.root')))
x4.SaveAs(str(residuals_dir / (name + '_pol3_slowEL1_fit.png')))
x4.Close()

pol3_params_slowEL1 = []
for i in range(4):
    pol3_params_slowEL1.append(h4.GetFunction("pol3").GetParameter(i))
for i in range(len(energy)):
    residuals3_slowEL1.append([energy[i] - (pol3_params_slowEL1[0] + pol3_params_slowEL1[1]*ch_slowEL1[i] + pol3_params_slowEL1[2]*ch_slowEL1[i]**2 + pol3_params_slowEL1[3]*ch_slowEL1[i]**3)])
with open(str(residuals_dir / (name + '_pol3_slowEL1_residuals.txt')), 'w') as f:
    for item in residuals3_slowEL1:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol3_slowEL1_residuals.txt')))

x5 = ROOT.TCanvas("c5", "c5", 800, 600)
x5.SetGrid()
h5 = ROOT.TH2D("h5", "h5", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h5.Fill(energy[i], residuals3_slowEL1[i][0])
h5.Draw()
h5.GetXaxis().SetTitle("Energy (keV)")
h5.GetYaxis().SetTitle("Residuals (keV)")
h5.SetStats(0)
h5.SetMarkerStyle(20)
h5.SetMarkerSize(2)
x5.Update()
x5.SaveAs(str(residuals_dir / (name + '_pol3_slowEL1_residuals.root')))
x5.SaveAs(str(residuals_dir / (name + '_pol3_slowEL1_residuals.png')))
x5.Close()

# ______________________________________________________________________________ best fit
# Get chi-square values for each fit
chi2_1 = h.GetFunction("pol1").GetChisquare()
chi2_2 = h2.GetFunction("pol2").GetChisquare()
chi2_3 = h4.GetFunction("pol3").GetChisquare()

# Handle zero chi-square values
if chi2_1 == 0:
    residuals1_slowEL1 = [np.nan]
if chi2_2 == 0:
    residuals2_slowEL1 = [np.nan]
if chi2_3 == 0:
    residuals3_slowEL1 = [np.nan]

# Determine which polynomial fit has the smallest residuals
fits_slowEL1 = [(residuals1_slowEL1, "pol1"), (residuals2_slowEL1, "pol2"), (residuals3_slowEL1, "pol3")]
min_fit_slowEL1 = min(fits_slowEL1, key=lambda x: sum(abs(val) if isinstance(val, (int, float)) else sum(abs(inner_val) for inner_val in val) for val in x[0]))

# Assign the best fit function name
bestfit = min_fit_slowEL1[1]
if min_fit_slowEL1[1] == 'pol1':
    p0_slowEL1 = pol1_params_slowEL1[0]
    p1_slowEL1 = pol1_params_slowEL1[1]
    p2_slowEL1 = 0
    p3_slowEL1 = 0
elif min_fit_slowEL1[1] == 'pol2':
    p0_slowEL1 = pol2_params_slowEL1[0]
    p1_slowEL1 = pol2_params_slowEL1[1]
    p2_slowEL1 = pol2_params_slowEL1[2]
    p3_slowEL1 = 0
elif min_fit_slowEL1[1] == 'pol3':
    p0_slowEL1 = pol3_params_slowEL1[0]
    p1_slowEL1 = pol3_params_slowEL1[1]
    p2_slowEL1 = pol3_params_slowEL1[2]
    p3_slowEL1 = pol3_params_slowEL1[3]


# ______________________________________________________________________________ residuals ch_fastEL1
residuals1_fastEL1 = []
x = ROOT.TCanvas("c", "c", 800, 600)
x.SetGrid()
pol1 = ROOT.TF1("pol1", "pol1", 0, 2000)
pol2 = ROOT.TF1("pol2", "pol2", 0, 2000)
pol3 = ROOT.TF1("pol3", "pol3", 0, 2000)

pol1.SetParameter(0, 0)
pol2.SetParameter(0, 0)
pol3.SetParameter(0, 0)

h = ROOT.TH2D("h", "h", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h.Fill(ch_fastEL1[i], energy[i])
h.Draw()
h.GetXaxis().SetTitle("Channel (a.u.)")
h.GetYaxis().SetTitle("Energy (keV)")
h.Fit("pol1")
h.GetFunction("pol1").SetLineColor(2)
h.GetFunction("pol1").SetLineWidth(1)
h.GetFunction("pol1").Draw("same")
h.SetStats(0)
h.SetMarkerStyle(20)
h.SetMarkerSize(2)
h.SetMarkerColor(1)
x.Update()
x.SaveAs(str(residuals_dir / (name + '_pol1_fastEL1_fit.root')))
x.SaveAs(str(residuals_dir / (name + '_pol1_fastEL1_fit.png')))
x.Close()

pol1_params_fastEL1 = []
for i in range(2):
    pol1_params_fastEL1.append(h.GetFunction("pol1").GetParameter(i))
for i in range(len(energy)):
    residuals1_fastEL1.append([energy[i] - (pol1_params_fastEL1[0] + pol1_params_fastEL1[1]*ch_fastEL1[i])])
with open(str(residuals_dir / (name + '_pol1_fastEL1_residuals.txt')), 'w') as f:
    for item in residuals1_fastEL1:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol1_fastEL1_residuals.txt')))

x1 = ROOT.TCanvas("c1", "c1", 800, 600)
x1.SetGrid()

h1 = ROOT.TH2D("h1", "h1", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h1.Fill(energy[i], residuals1_fastEL1[i][0])
h1.Draw()
h1.GetXaxis().SetTitle("Energy (keV)")
h1.GetYaxis().SetTitle("Residuals (keV)")
h1.SetStats(0)
h1.SetMarkerStyle(20)
h1.SetMarkerSize(2)
x1.Update()
x1.SaveAs(str(residuals_dir / (name + '_pol1_fastEL1_residuals.root')))
x1.SaveAs(str(residuals_dir / (name + '_pol1_fastEL1_residuals.png')))
x1.Close()

residuals2_fastEL1 = []
x2 = ROOT.TCanvas("c2", "c2", 800, 600)
x2.SetGrid()
h2 = ROOT.TH2D("h2", "h2", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h2.Fill(ch_fastEL1[i], energy[i])
h2.Draw()
h2.GetXaxis().SetTitle("Channel (a.u.)")
h2.GetYaxis().SetTitle("Energy (keV)")
h2.Fit("pol2")
h2.GetFunction("pol2").SetLineColor(2)
h2.GetFunction("pol2").SetLineWidth(1)
h2.GetFunction("pol2").Draw("same")
h2.SetStats(0)
h2.SetMarkerStyle(20)
h2.SetMarkerSize(2)
h2.SetMarkerColor(1)
x2.Update()
x2.SaveAs(str(residuals_dir / (name + '_pol2_fastEL1_fit.root')))
x2.SaveAs(str(residuals_dir / (name + '_pol2_fastEL1_fit.png')))
x2.Close()

pol2_params_fastEL1 = []
for i in range(3):
    pol2_params_fastEL1.append(h2.GetFunction("pol2").GetParameter(i))
for i in range(len(energy)):
    residuals2_fastEL1.append([energy[i] - (pol2_params_fastEL1[0] + pol2_params_fastEL1[1]*ch_fastEL1[i] + pol2_params_fastEL1[2]*ch_fastEL1[i]**2)])
with open(str(residuals_dir / (name + '_pol2_fastEL1_residuals.txt')), 'w') as f:
    for item in residuals2_fastEL1:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol2_fastEL1_residuals.txt')))

x3 = ROOT.TCanvas("c3", "c3", 800, 600)
x3.SetGrid()
h3 = ROOT.TH2D("h3", "h3", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h3.Fill(energy[i], residuals2_fastEL1[i][0])
h3.Draw()
h3.GetXaxis().SetTitle("Energy (keV)")
h3.GetYaxis().SetTitle("Residuals (keV)")
h3.SetStats(0)
h3.SetMarkerStyle(20)
h3.SetMarkerSize(2)
x3.Update()
x3.SaveAs(str(residuals_dir / (name + '_pol2_fastEL1_residuals.root')))
x3.SaveAs(str(residuals_dir / (name + '_pol2_fastEL1_residuals.png')))
x3.Close()

residuals3_fastEL1 = []
x4 = ROOT.TCanvas("c4", "c4", 800, 600)
x4.SetGrid()
h4 = ROOT.TH2D("h4", "h4", 16350,0,16350,2000,0,2000)
for i in range(len(energy)):
    h4.Fill(ch_fastEL1[i], energy[i])
h4.Draw()
h4.GetXaxis().SetTitle("Channel (a.u.)")
h4.GetYaxis().SetTitle("Energy (keV)")
h4.Fit("pol3")
h4.GetFunction("pol3").SetLineColor(2)
h4.GetFunction("pol3").SetLineWidth(1)
h4.GetFunction("pol3").Draw("same")
h4.SetStats(0)
h4.SetMarkerStyle(20)
h4.SetMarkerSize(2)
h4.SetMarkerColor(1)
x4.Update()
x4.SaveAs(str(residuals_dir / (name + '_pol3_fastEL1_fit.root')))
x4.SaveAs(str(residuals_dir / (name + '_pol3_fastEL1_fit.png')))
x4.Close()

pol3_params_fastEL1 = []
for i in range(4):
    pol3_params_fastEL1.append(h4.GetFunction("pol3").GetParameter(i))
for i in range(len(energy)):
    residuals3_fastEL1.append([energy[i] - (pol3_params_fastEL1[0] + pol3_params_fastEL1[1]*ch_fastEL1[i] + pol3_params_fastEL1[2]*ch_fastEL1[i]**2 + pol3_params_fastEL1[3]*ch_fastEL1[i]**3)])
with open(str(residuals_dir / (name + '_pol3_fastEL1_residuals.txt')), 'w') as f:
    for item in residuals3_fastEL1:
        f.write("%s\n" % item)
print('Residuals saved to {}'.format(residuals_dir / (name + '_pol3_fastEL1_residuals.txt')))

x5 = ROOT.TCanvas("c5", "c5", 800, 600)
x5.SetGrid()
h5 = ROOT.TH2D("h5", "h5", 2000, 0, 2000, 50,-25,25)
for i in range(len(energy)):
    h5.Fill(energy[i], residuals3_fastEL1[i][0])
h5.Draw()
h5.GetXaxis().SetTitle("Energy (keV)")
h5.GetYaxis().SetTitle("Residuals (keV)")
h5.SetStats(0)
h5.SetMarkerStyle(20)
h5.SetMarkerSize(2)
x5.Update()
x5.SaveAs(str(residuals_dir / (name + '_pol3_fastEL1_residuals.root')))
x5.SaveAs(str(residuals_dir / (name + '_pol3_fastEL1_residuals.png')))
x5.Close()

# ______________________________________________________________________________ best fit
# Get chi-square values for each fit
chi2_1 = h.GetFunction("pol1").GetChisquare()
chi2_2 = h2.GetFunction("pol2").GetChisquare()
chi2_3 = h4.GetFunction("pol3").GetChisquare()

# Handle zero chi-square values
if chi2_1 == 0:
    residuals1_fastEL1 = [np.nan]
if chi2_2 == 0:
    residuals2_fastEL1 = [np.nan]
if chi2_3 == 0:
    residuals3_fastEL1 = [np.nan]

# Determine which polynomial fit has the smallest residuals
fits_fastEL1 = [(residuals1_fastEL1, "pol1"), (residuals2_fastEL1, "pol2"), (residuals3_fastEL1, "pol3")]
min_fit_fastEL1 = min(fits_fastEL1, key=lambda x: sum(abs(val) if isinstance(val, (int, float)) else sum(abs(inner_val) for inner_val in val) for val in x[0]))

# Assign the best fit function name
bestfit = min_fit_fastEL1[1]
if min_fit_fastEL1[1] == 'pol1':
    p0_fastEL1 = pol1_params_fastEL1[0]
    p1_fastEL1 = pol1_params_fastEL1[1]
    p2_fastEL1 = 0
    p3_fastEL1 = 0
elif min_fit_fastEL1[1] == 'pol2':
    p0_fastEL1 = pol2_params_fastEL1[0]
    p1_fastEL1 = pol2_params_fastEL1[1]
    p2_fastEL1 = pol2_params_fastEL1[2]
    p3_fastEL1 = 0
elif min_fit_fastEL1[1] == 'pol3':
    p0_fastEL1 = pol3_params_fastEL1[0]
    p1_fastEL1 = pol3_params_fastEL1[1]
    p2_fastEL1 = pol3_params_fastEL1[2]
    p3_fastEL1 = pol3_params_fastEL1[3]

# ______________________________________________________________________________ text file of fit parameters for the best fit
with open(str(residuals_dir / (name + '_bestfit_params.txt')), 'w') as f:
    f.write('Best fit: {}\n'.format(min_fit_fastEL0[1]))
    # write polynomial equation
    f.write('y = p0 + p1*x + p2*x^2 + p3*x^3\n')
    f.write('FastEL0: p3 = {}, p2 = {}, p1 = {}, p0 = {}\n'.format(p3_fastEL0, p2_fastEL0, p1_fastEL0, p0_fastEL0))
    f.write('SlowEL0: p3 = {}, p2 = {}, p1 = {}, p0 = {}\n'.format(p3_slowEL0, p2_slowEL0, p1_slowEL0, p0_slowEL0))

    f.write('FastEL1: p3 = {}, p2 = {}, p1 = {}, p0 = {}\n'.format(p3_fastEL1, p2_fastEL1, p1_fastEL1, p0_fastEL1))
    f.write('SlowEL1: p3 = {}, p2 = {}, p1 = {}, p0 = {}\n\n'.format(p3_slowEL1, p2_slowEL1, p1_slowEL1, p0_slowEL1))

    f.write('Both detectors E0 E1 - fast p2 {}, {}\n'.format(p2_fastEL0, p2_fastEL1))
    f.write('Both detectors E0 E1 - fast p1 {}, {}\n'.format(p1_fastEL0, p1_fastEL1))
    f.write('Both detectors E0 E1 - fast p0 {}, {}\n\n'.format(p0_fastEL0, p0_fastEL1))

    f.write('Both detectors E0 E1 - slow p2 {}, {}\n'.format(p2_slowEL0, p2_slowEL1))
    f.write('Both detectors E0 E1 - slow p1 {}, {}\n'.format(p1_slowEL0, p1_slowEL1))
    f.write('Both detectors E0 E1 - slow p0 {}, {}\n\n'.format(p0_slowEL0, p0_slowEL1))


print('Fit parameters saved to {}'.format(residuals_dir / (name + '_bestfit_params.txt')))


# ______________________________________________________________________________ calibration
# calibrate the data and update the open root file
calib = ROOT.TFile(file, "UPDATE")
t = calib.Get("LaBrData")

slowEL0calibpol1 = ROOT.TH1D("slowEL0calib_pol1", "slowEL0calib_pol1", 4000, 0, 4000)
fastEL0calibpol1 = ROOT.TH1D("fastEL0calib_pol1", "fastEL0calib_pol1", 4000, 0, 4000)
slowEL0calibpol2 = ROOT.TH1D("slowEL0calib_pol2", "slowEL0calib_pol2", 4000, 0, 4000)
fastEL0calibpol2 = ROOT.TH1D("fastEL0calib_pol2", "fastEL0calib_pol2", 4000, 0, 4000)
slowEL0calibpol3 = ROOT.TH1D("slowEL0calib_pol3", "slowEL0calib_pol3", 4000, 0, 4000)
fastEL0calibpol3 = ROOT.TH1D("fastEL0calib_pol3", "fastEL0calib_pol3", 4000, 0, 4000)

slowEL1calibpol1 = ROOT.TH1D("slowEL1calib_pol1", "slowEL1calib_pol1", 4000, 0, 4000)
fastEL1calibpol1 = ROOT.TH1D("fastEL1calib_pol1", "fastEL1calib_pol1", 4000, 0, 4000)
slowEL1calibpol2 = ROOT.TH1D("slowEL1calib_pol2", "slowEL1calib_pol2", 4000, 0, 4000)
fastEL1calibpol2 = ROOT.TH1D("fastEL1calib_pol2", "fastEL1calib_pol2", 4000, 0, 4000)
slowEL1calibpol3 = ROOT.TH1D("slowEL1calib_pol3", "slowEL1calib_pol3", 4000, 0, 4000)
fastEL1calibpol3 = ROOT.TH1D("fastEL1calib_pol3", "fastEL1calib_pol3", 4000, 0, 4000)

for i in range(t.GetEntries()):
    t.GetEntry(i)
    if t.slowEL0 > 16350 or t.slowEL0 > 16350:
        continue
    if scaleslowEL0 > 1:
        scale = int(scaleslowEL0)
        scale = scale/2
        ranges = [-scale, scale]
        slowEL0calibpol1.Fill((t.slowEL0*pol1_params_slowEL0[1] + pol1_params_slowEL0[0]) + r.Uniform(-scale,scale))
        slowEL0calibpol2.Fill((t.slowEL0**2*pol2_params_slowEL0[2] + t.slowEL0*pol2_params_slowEL0[1] + pol2_params_slowEL0[0]) + r.Uniform(-scale,scale))
        slowEL0calibpol3.Fill((t.slowEL0**3*pol3_params_slowEL0[3] + t.slowEL0**2*pol3_params_slowEL0[2] + t.slowEL0*pol3_params_slowEL0[1] + pol3_params_slowEL0[0]) + r.Uniform(-scale,scale))
    else:
        slowEL0calibpol1.Fill((t.slowEL0*pol1_params_slowEL0[1] + pol1_params_slowEL0[0]))
        slowEL0calibpol2.Fill((t.slowEL0**2*pol2_params_slowEL0[2] + t.slowEL0*pol2_params_slowEL0[1] + pol2_params_slowEL0[0]))
        slowEL0calibpol3.Fill((t.slowEL0**3*pol3_params_slowEL0[3] + t.slowEL0**2*pol3_params_slowEL0[2] + t.slowEL0*pol3_params_slowEL0[1] + pol3_params_slowEL0[0]))

    if scalefastEL0 > 1:
        scale = int(scalefastEL0)
        scale = scale/2
        ranges = [-scale, scale]
        fastEL0calibpol1.Fill((t.fastEL0*pol1_params_fastEL0[1] + pol1_params_fastEL0[0]) + r.Uniform(-scale,scale))
        fastEL0calibpol2.Fill((pol2_params_fastEL0[2]*t.fastEL0**2 + pol2_params_fastEL0[1]*t.fastEL0 + pol2_params_fastEL0[0]) + r.Uniform(-scale,scale))
        fastEL0calibpol3.Fill((pol3_params_fastEL0[3]*t.fastEL0**3 + pol3_params_fastEL0[2]*t.fastEL0**2 + pol3_params_fastEL0[1]*t.fastEL0 + pol3_params_fastEL0[0]) + r.Uniform(-scale,scale))
    else:
        fastEL0calibpol1.Fill((t.fastEL0*pol1_params_fastEL0[1] + pol1_params_fastEL0[0]))
        fastEL0calibpol2.Fill((pol2_params_fastEL0[2]*t.fastEL0**2 + pol2_params_fastEL0[1]*t.fastEL0 + pol2_params_fastEL0[0]))
        fastEL0calibpol3.Fill((pol3_params_fastEL0[3]*t.fastEL0**3 + pol3_params_fastEL0[2]*t.fastEL0**2 + pol3_params_fastEL0[1]*t.fastEL0 + pol3_params_fastEL0[0]))

    if scaleslowEL1 > 1:
        scale = int(scaleslowEL1)
        scale = scale/2
        ranges = [-scale, scale]
        slowEL1calibpol1.Fill((t.slowEL1*pol1_params_slowEL1[1] + pol1_params_slowEL1[0]) + r.Uniform(-scale,scale))
        slowEL1calibpol2.Fill((pol2_params_slowEL1[2]*t.slowEL1**2 + pol2_params_slowEL1[1]*t.slowEL1 + pol2_params_slowEL1[0]) + r.Uniform(-scale,scale))
        slowEL1calibpol3.Fill((pol3_params_slowEL1[3]*t.slowEL1**3 + pol3_params_slowEL1[2]*t.slowEL1**2 + pol3_params_slowEL1[1]*t.slowEL1 + pol3_params_slowEL1[0]) + r.Uniform(-scale,scale))
    else:
        slowEL1calibpol1.Fill((t.slowEL1*pol1_params_slowEL1[1] + pol1_params_slowEL1[0]))
        slowEL1calibpol2.Fill((pol2_params_slowEL1[2]*t.slowEL1**2 + pol2_params_slowEL1[1]*t.slowEL1 + pol2_params_slowEL1[0]))
        slowEL1calibpol3.Fill((pol3_params_slowEL1[3]*t.slowEL1**3 + pol3_params_slowEL1[2]*t.slowEL1**2 + pol3_params_slowEL1[1]*t.slowEL1 + pol3_params_slowEL1[0]))
    
    if scalefastEL1 > 1:
        scale = int(scalefastEL1)
        scale = scale/2
        ranges = [-scale, scale]
        fastEL1calibpol1.Fill((t.fastEL1*pol1_params_fastEL1[1] + pol1_params_fastEL1[0]) + r.Uniform(-scale,scale))
        fastEL1calibpol2.Fill((pol2_params_fastEL1[2]*t.fastEL1**2 + pol2_params_fastEL1[1]*t.fastEL1 + pol2_params_fastEL1[0]) + r.Uniform(-scale,scale))
        fastEL1calibpol3.Fill((pol3_params_fastEL1[3]*t.fastEL1**3 + pol3_params_fastEL1[2]*t.fastEL1**2 + pol3_params_fastEL1[1]*t.fastEL1 + pol3_params_fastEL1[0]) + r.Uniform(-scale,scale))
    else:
        fastEL1calibpol1.Fill((t.fastEL1*pol1_params_fastEL1[1] + pol1_params_fastEL1[0]))
        fastEL1calibpol2.Fill((pol2_params_fastEL1[2]*t.fastEL1**2 + pol2_params_fastEL1[1]*t.fastEL1 + pol2_params_fastEL1[0]))
        fastEL1calibpol3.Fill((pol3_params_fastEL1[3]*t.fastEL1**3 + pol3_params_fastEL1[2]*t.fastEL1**2 + pol3_params_fastEL1[1]*t.fastEL1 + pol3_params_fastEL1[0]))

slowEL0calibpol1.Write()
fastEL0calibpol1.Write()        
slowEL0calibpol2.Write()
fastEL0calibpol2.Write()
slowEL0calibpol3.Write()
fastEL0calibpol3.Write()

slowEL1calibpol1.Write()
fastEL1calibpol1.Write()
slowEL1calibpol2.Write()
fastEL1calibpol2.Write()

c10 = ROOT.TCanvas("c10", "c10", 800, 600)
c10.SetGrid()
slowEL0calibpol1.Draw()
slowEL0calibpol1.GetXaxis().SetTitle("Energy (keV)")
slowEL0calibpol1.GetYaxis().SetTitle("Counts (1/keV)")
slowEL0calibpol1.SetLineColor(1) # red
slowEL0calibpol1.SetLineWidth(1)
slowEL0calibpol1.SetStats(0)
c10.Update()
c10.SaveAs(str(calib_dir / (name + '_slowEL0calib_pol1.root')))
c10.SaveAs(str(calib_dir / (name + '_slowEL0calib_pol1.png')))

c11 = ROOT.TCanvas("c11", "c11", 800, 600)
c11.SetGrid()
fastEL0calibpol1.Draw()
fastEL0calibpol1.GetXaxis().SetTitle("Energy (keV)")
fastEL0calibpol1.GetYaxis().SetTitle("Counts (1/keV)")
fastEL0calibpol1.SetLineColor(1) # red
fastEL0calibpol1.SetLineWidth(1)
fastEL0calibpol1.SetStats(0)
c11.Update()
c11.SaveAs(str(calib_dir / (name + '_fastEL0calib_pol1.root')))
c11.SaveAs(str(calib_dir / (name + '_fastEL0calib_pol1.png')))

c12 = ROOT.TCanvas("c12", "c12", 800, 600)
c12.SetGrid()
slowEL0calibpol2.Draw()
slowEL0calibpol2.GetXaxis().SetTitle("Energy (keV)")
slowEL0calibpol2.GetYaxis().SetTitle("Counts (1/keV)")
slowEL0calibpol2.SetLineColor(1) # red
slowEL0calibpol2.SetLineWidth(1)
slowEL0calibpol2.SetStats(0)
c12.Update()
c12.SaveAs(str(calib_dir / (name + '_slowEL0calib_pol2.root')))
c12.SaveAs(str(calib_dir / (name + '_slowEL0calib_pol2.png')))

c15 = ROOT.TCanvas("c15", "c15", 800, 600)
c15.SetGrid()
fastEL0calibpol2.Draw()
fastEL0calibpol2.GetXaxis().SetTitle("Energy (keV)")
fastEL0calibpol2.GetYaxis().SetTitle("Counts (1/keV)")
fastEL0calibpol2.SetLineColor(1) # red
fastEL0calibpol2.SetLineWidth(1)
fastEL0calibpol2.SetStats(0)
c15.Update()
c15.SaveAs(str(calib_dir / (name + '_fastEL0calib_pol2.root')))
c15.SaveAs(str(calib_dir / (name + '_fastEL0calib_pol2.png')))

c13 = ROOT.TCanvas("c13", "c13", 800, 600)
c13.SetGrid()
slowEL0calibpol3.Draw()
slowEL0calibpol3.GetXaxis().SetTitle("Energy (keV)")
slowEL0calibpol3.GetYaxis().SetTitle("Counts (1/keV)")
slowEL0calibpol3.SetLineColor(1) # red
slowEL0calibpol3.SetLineWidth(1)
slowEL0calibpol3.SetStats(0)
c13.Update()
c13.SaveAs(str(calib_dir / (name + '_slowEL0calib_pol3.root')))
c13.SaveAs(str(calib_dir / (name + '_slowEL0calib_pol3.png')))

c14 = ROOT.TCanvas("c14", "c14", 800, 600)
c14.SetGrid()
fastEL0calibpol3.Draw()
fastEL0calibpol3.GetXaxis().SetTitle("Energy (keV)")
fastEL0calibpol3.GetYaxis().SetTitle("Counts (1/keV)")
fastEL0calibpol3.SetLineColor(1) # red
fastEL0calibpol3.SetLineWidth(1)
fastEL0calibpol3.SetStats(0)
c14.Update()
c14.SaveAs(str(calib_dir / (name + '_fastEL0calib_pol3.root')))
c14.SaveAs(str(calib_dir / (name + '_fastEL0calib_pol3.png')))

c16 = ROOT.TCanvas("c16", "c16", 800, 600)
c16.SetGrid()
slowEL1calibpol1.Draw()
slowEL1calibpol1.GetXaxis().SetTitle("Energy (keV)")
slowEL1calibpol1.GetYaxis().SetTitle("Counts (1/keV)")
slowEL1calibpol1.SetLineColor(1) # red
slowEL1calibpol1.SetLineWidth(1)
slowEL1calibpol1.SetStats(0)
c16.Update()
c16.SaveAs(str(calib_dir / (name + '_slowEL1calib_pol1.root')))
c16.SaveAs(str(calib_dir / (name + '_slowEL1calib_pol1.png')))

c17 = ROOT.TCanvas("c17", "c17", 800, 600)
c17.SetGrid()
fastEL1calibpol1.Draw()
fastEL1calibpol1.GetXaxis().SetTitle("Energy (keV)")
fastEL1calibpol1.GetYaxis().SetTitle("Counts (1/keV)")
fastEL1calibpol1.SetLineColor(1) # red
fastEL1calibpol1.SetLineWidth(1)
fastEL1calibpol1.SetStats(0)
c17.Update()
c17.SaveAs(str(calib_dir / (name + '_fastEL1calib_pol1.root')))
c17.SaveAs(str(calib_dir / (name + '_fastEL1calib_pol1.png')))

c18 = ROOT.TCanvas("c18", "c18", 800, 600)
c18.SetGrid()
slowEL1calibpol2.Draw()
slowEL1calibpol2.GetXaxis().SetTitle("Energy (keV)")
slowEL1calibpol2.GetYaxis().SetTitle("Counts (1/keV)")
slowEL1calibpol2.SetLineColor(1) # red
slowEL1calibpol2.SetLineWidth(1)
slowEL1calibpol2.SetStats(0)
c18.Update()
c18.SaveAs(str(calib_dir / (name + '_slowEL1calib_pol2.root')))
c18.SaveAs(str(calib_dir / (name + '_slowEL1calib_pol2.png')))

c19 = ROOT.TCanvas("c19", "c19", 800, 600)
c19.SetGrid()
fastEL1calibpol2.Draw()
fastEL1calibpol2.GetXaxis().SetTitle("Energy (keV)")
fastEL1calibpol2.GetYaxis().SetTitle("Counts (1/keV)")
fastEL1calibpol2.SetLineColor(1) # red
fastEL1calibpol2.SetLineWidth(1)
fastEL1calibpol2.SetStats(0)
c19.Update()
c19.SaveAs(str(calib_dir / (name + '_fastEL1calib_pol2.root')))
c19.SaveAs(str(calib_dir / (name + '_fastEL1calib_pol2.png')))

c20 = ROOT.TCanvas("c20", "c20", 800, 600)
c20.SetGrid()
slowEL1calibpol3.Draw()
slowEL1calibpol3.GetXaxis().SetTitle("Energy (keV)")
slowEL1calibpol3.GetYaxis().SetTitle("Counts (1/keV)")
slowEL1calibpol3.SetLineColor(1) # red
slowEL1calibpol3.SetLineWidth(1)
slowEL1calibpol3.SetStats(0)
c20.Update()
c20.SaveAs(str(calib_dir / (name + '_slowEL1calib_pol3.root')))
c20.SaveAs(str(calib_dir / (name + '_slowEL1calib_pol3.png')))
c21 = ROOT.TCanvas("c21", "c21", 800, 600)
c21.SetGrid()
fastEL1calibpol3.Draw()
fastEL1calibpol3.GetXaxis().SetTitle("Energy (keV)")
fastEL1calibpol3.GetYaxis().SetTitle("Counts (1/keV)")
fastEL1calibpol3.SetLineColor(1) # red
fastEL1calibpol3.SetLineWidth(1)
fastEL1calibpol3.SetStats(0)
c21.Update()
c21.SaveAs(str(calib_dir / (name + '_fastEL1calib_pol3.root')))
c21.SaveAs(str(calib_dir / (name + '_fastEL1calib_pol3.png')))

calib.Close()

print('Calibration complete.')
print(' ')
print(' ')
# Print the fit with the smallest residuals
print(' ')
print(' ')
print(f"The {min_fit_slowEL0[1]} fit has the smallest residuals - slowEL0.")
print(' ')
print(' ')
# Print the fit with the smallest residuals
print(' ')
print(' ')
print(f"The {min_fit_fastEL0[1]} fit has the smallest residuals - fastEL0.")
print(' ')
print(' ')
print(f"The {min_fit_slowEL1[1]} fit has the smallest residuals - slowEL1.")
print(' ')
print(' ')
# Print the fit with the smallest residuals
print(' ')
print(' ')
print(f"The {min_fit_fastEL1[1]} fit has the smallest residuals - fastEL1.")
print(' ')
print(' ')