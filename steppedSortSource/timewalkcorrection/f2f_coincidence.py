import pandas as pd
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ROOT
import os
from concurrent.futures import ThreadPoolExecutor
import imageio

plt.rcParams.update({'font.size': 18})
# Make axis labels bold
plt.rcParams["axes.labelweight"] = "bold"
# Label size
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)

#____________________________________________________________
# e.g. timesyncLaBrPOLARIS % ~/miniconda/envs/my_root_env/bin/python f2f_coincidence.py ~/Documents/PhD/exp/2024/timing_shaping_characterisation_june2024/CC_polaris_2inch/run7/ (path to parquet file)
#____________________________________________________________

# Function to find backward events where LaBr3 event precedes a POLARIS event
def find_backward_events_chunk(chunk_indices, data_df, time_window, energy_lower, energy_upper):
    backward_events = []
    for index_det5 in chunk_indices:
        for index_det0 in range(index_det5 - 1, -1, -1):
            time_difference = data_df.loc[index_det5, 'time'] - data_df.loc[index_det0, 'time']
            if time_difference > time_window:
                break
            #if data_df.loc[index_det0, 'detectorID'] == 0 and energy_lower <= data_df.loc[index_det0, 'energy'] <= energy_upper:
            else:
                backward_events.append((time_difference, data_df.loc[index_det0, 'energy'], data_df.loc[index_det5, 'energy']))
    return backward_events


# Function to find backward events where LaBr3 event follows a POLARIS event
def find_forward_events_chunk(chunk):
    forward_events = []
    time_window = 600e3  # 500 microseconds in nanoseconds

    for idy in chunk:
        polaristime = times[idy]
        polarisenergy = energies[idy]

        # Find the index in detector0_indices where the labr3time is just greater than the polaristime using binary search
        idx = np.searchsorted(times[detector0_indices], polaristime)
        
        while idx < len(detector0_indices) and times[detector0_indices[idx]] > polaristime:
            labr3time = times[detector0_indices[idx]]
            labr3energy = energies[detector0_indices[idx]]
            time_diff = labr3time - polaristime
            
            if 0 < time_diff <= time_window:
                forward_events.append((time_diff, labr3energy, polarisenergy))
                break
            idx += 1

    return forward_events
#____________________________________________________________________________

data_dir = sys.argv[1]
save_dir = data_dir + "/plots/22Na_f2f_coinc/"

# Check if save_dir exists, if not, create it
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Directory to save the PNG images
png_dir = os.path.join(save_dir, 'png')
if not os.path.exists(png_dir):
    os.makedirs(png_dir)


run_num = int(re.search('run(\d+)', data_dir).group(1))
df = pd.read_parquet(data_dir + 'run' + str(run_num) + '_mergedLaBrPOLARISdata.parquet')
df = df.sort_values(by='time', ascending=True)

# Define the energy range for 511 keV
energy_tolerance = 0.035 * 511  # 15% tolerance
energy_lower = 511 - energy_tolerance
energy_upper = 511 + energy_tolerance
data_df = df[(df['detectorID'].isin([0, 5])) ].copy()
#data_df = data_df[(data_df['energy'] >= energy_lower) & (data_df['energy'] <= energy_upper)]

# quick 511 calibration
df.loc[df['detectorID'] == 0, 'energy'] *= 0.9734
df.loc[df['detectorID'] == 5, 'energy'] *= 1.0081

# Sort by time to ensure correct order
df.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
df.reset_index(inplace=True, drop=True)

#print the time difference between detector 0 and detector 5 consecutive events that that are within +-3.5% of 511 keV
""" valid_indices = df[(df['detectorID'] == 5) & (df['detectorID'].shift(-1) == 0)].index
time_diffs = df.loc[valid_indices + 1, 'time'].values - df.loc[valid_indices, 'time'].values
energies = df.loc[valid_indices + 1, 'energy'].values
time_diffs = time_diffs[abs(energies - 511) <= energy_tolerance]
print("time_diffs", time_diffs)

# 1d histogram of time differences
canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
hist = ROOT.TH1D("hist", "Time difference between LaBr3 and POLARIS events", 100, 0, 1000000)
for time_diff in time_diffs:
    hist.Fill(time_diff)
hist.Draw()
hist.GetXaxis().SetTitle("Time difference (ns)")
hist.GetYaxis().SetTitle("Counts (1/10 #mus)")
canvas.Update()
canvas.SaveAs(f"{save_dir}time_difference_histogram.root") """

# Apply time correction in increments and plot
""" time_offsets = np.arange(0, 301e3, 10e3)  # from -10 microseconds to 300 microseconds in steps of 10 microseconds

# Create an empty list to store file names of saved plots
plot_filenames = []
png_filenames = []

for offset in time_offsets:
    df_corrected = df.copy()
    df_corrected.loc[df_corrected['detectorID'] == 5, 'time'] -= offset  # adjust time
    df_corrected.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
    df_corrected.reset_index(inplace=True, drop=True)

    # Identify valid events
    valid_indices = df_corrected[(df_corrected['detectorID'] == 0) & (df_corrected['detectorID'].shift(-1) == 5)].index # valid indices where LaBr3 event precedes a POLARIS event

    # Create boolean mask for time difference condition
    time_condition = abs(df_corrected.loc[valid_indices + 1, 'time'].values - df_corrected.loc[valid_indices, 'time'].values) <= 50e3

    # Filter valid_indices based on the boolean mask
    valid_indices = valid_indices[time_condition]
    # Extract energies for valid events
    labr3_energies = df_corrected.loc[valid_indices, 'energy']
    polaris_energies = df_corrected.loc[valid_indices + 1, 'energy']

    # Plot the 2D energy matrix
    canvas_2d = ROOT.TCanvas("canvas_2d", "canvas_2d", 800, 600)
    hist_2d = ROOT.TH2D("hist_2d", "Polaris offset: {} #mus".format(offset/1000), 130, 0, 1300, 130, 0, 1300)
    for labr3_energy, polaris_energy in zip(labr3_energies, polaris_energies):
        hist_2d.Fill(polaris_energy, labr3_energy)
    hist_2d.Draw("COLZ")
    hist_2d.SetMinimum(1)
    hist_2d.SetMaximum(2)
    hist_2d.GetXaxis().SetTitle("Polaris Energy (keV)")
    hist_2d.GetYaxis().SetTitle("LaBr_{3}:Ce Energy (keV)")
    canvas_2d.Update()

    # Save the plot as a PNG file
    plot_filename = f"{save_dir}energy_energy_matrix_offset_{int(offset/1000)}us.png"  # save with microseconds in filename
    canvas_2d.SaveAs(plot_filename)
    plot_filenames.append(plot_filename)

    png_filename = os.path.join(png_dir, f"energy_energy_matrix_offset_{int(offset/1000)}us.png")
    canvas_2d.SaveAs(png_filename)
    plot_filenames.append(png_filename) """


# calculate how many polaris events within 15 keV of 511 keV occur 5us, 10us, 15us, 20us, 25us, 30us, 35us, 40us, 45us, 50us ...500 us before a LaBr3 event within 50keV of 511 keV
#ACTUAL COINCIDENCE CODE:
time_offsets = np.arange(-101e3, 601e3, 100)  # from -500 microseconds to 500 microseconds in steps of 1 microsecond

# Create an empty list to store file names of saved plots
plot_filenames = []
png_filenames = []

df_corrected = df.copy()
df_corrected.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
df_corrected.reset_index(inplace=True, drop=True)

# print the ratio of detector 0 to detector 5 events within 22 keV of 511 keV in polaris and 64 keV of 511 keV in LaBr3
print("Ratio: ", len(df_corrected[(df_corrected['detectorID'] == 0) & (abs(df_corrected['energy'] - 511) <= 64)])/len(df_corrected[(df_corrected['detectorID'] == 5) & (abs(df_corrected['energy'] - 511) <= 22)]))
coincidence_rates = []
# Process each time offset
for offset in time_offsets:
    count = 0
    # Apply time offset correction
    df_corrected['time'] = df['time'].copy()
    df_corrected.loc[df_corrected['detectorID'] == 5, 'time'] -= offset  # adjust time

    # Filter events within the specified energy range
    polaris_energy_range = 22  # 6 sigma
    labr3_energy_range = 64  # 6 sigma
    df_corrected = df_corrected[(df_corrected['detectorID'] == 5) & (abs(df_corrected['energy'] - 511) <= polaris_energy_range) | (df_corrected['detectorID'] == 0) & (abs(df_corrected['energy'] - 511) <= labr3_energy_range)]
    df_corrected.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
    df_corrected.reset_index(inplace=True, drop=True)
    
    #calculate the time difference between detector 0 and detector 5 consecutive events that that are within 300 ns of each other. this means that an event might not be exactly above or below the other event but within 300 ns of each other
    if (offset >= 0):
        count = 0
        for index_det5 in range(len(df_corrected)):
            for index_det0 in range(index_det5 - 1, -1, -1):
                time_difference = (df_corrected.loc[index_det5, 'time'] - df_corrected.loc[index_det0, 'time'])
                if abs(time_difference) > 1e3:
                    break
                if df_corrected.loc[index_det0, 'detectorID'] == 0:
                    count += 1
                    coincidence_rates.append((offset, count))
                    break
    else:
        count = 0
        for index_det0 in range(len(df_corrected)):
            for index_det5 in range(index_det0 -1, -1, -1):
                time_difference = (df_corrected.loc[index_det0, 'time'] - df_corrected.loc[index_det5, 'time'])
                if abs(time_difference) > 1e3:
                    break
                if df_corrected.loc[index_det5, 'detectorID'] == 0:
                    count += 1
                    coincidence_rates.append((offset, count))
                    break

coincidence_rates = np.array(coincidence_rates)
# divide the number of coincidences by the time window to get the coincidence rate
coincidence_rates[:, 1] /= 1e3  # divide by 300 microseconds to get the coincidence rate in 1/10 microseconds

# Plot the number of coincidences vs offset
canvas4 = ROOT.TCanvas("canvas4", "canvas4", 800, 600) 
coinc_num  = ROOT.TH1D("coinc_num", "Number of coincidences for each offset", 700000, -100000, 600000) # from - 500 us to 500 us in steps of 1 us
for offset, rate in coincidence_rates:
    coinc_num.Fill(offset, rate)
coinc_num.Draw()
coinc_num.GetXaxis().SetTitle("Time Delay (ns)")
coinc_num.GetYaxis().SetTitle("Coincidence rate (1/10 #mus)")
canvas4.Update()
canvas4.SaveAs(f"{save_dir}coincidence_rate_histogram.root")



df_processed = df.copy()
df_processed.loc[df_processed['detectorID'] == 5, 'time'] -= 145000  # adjust time
df_processed.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
df_processed.reset_index(inplace=True, drop=True)

# calculate the time difference between detector 0 and detector 5 consecutive events that that are within 500 ns of each other
coincidence_rates = []
same_energy = []
index_det5 = 0
while index_det5 < len(df_processed):
    for index_det0 in range(index_det5 - 1, -1, -1):
        time_difference = (df_processed.loc[index_det5, 'time'] - df_processed.loc[index_det0, 'time'])
        energy_polaris = df_processed.loc[index_det5, 'energy']
        energy_labr = df_processed.loc[index_det0, 'energy']
        if ((abs(time_difference) > 1)  | (abs(time_difference) == 0)):
            break
        if df_processed.loc[index_det0, 'detectorID'] == 0:  # if the event is a LaBr3 event
            coincidence_rates.append((time_difference, energy_polaris, energy_labr))
            # Drop the rows to avoid double counting
            df_processed.drop([index_det5, index_det0], inplace=True)
            # Reset the index to maintain consistency
            df_processed.reset_index(drop=True, inplace=True)
            index_det5 -= 1  # Decrement index_det5 because we dropped a row before it
            if (energy_labr == energy_polaris):
                same_energy.append((energy_labr, time_difference))
                #print("energy", energy_labr, "time", time_difference)
            break  # Exit the inner loop after finding a coincidence
    index_det5 += 1  # Move to the next index

coincidence_rates = np.array(coincidence_rates)
# Plot the 2D energy matrix
canvas_2d = ROOT.TCanvas("canvas_2d", "canvas_2d", 800, 600)
hist_2d = ROOT.TH2D("hist_2d", "Time window: 500 ns", 130, 0, 1300, 130, 0, 1300)
for energy_polaris, energy_labr in coincidence_rates[:, 1:]:
    hist_2d.Fill(energy_polaris, energy_labr)
hist_2d.Draw("COLZ")
hist_2d.SetMinimum(1)
hist_2d.GetXaxis().SetTitle("Polaris Energy (keV)")
hist_2d.GetYaxis().SetTitle("LaBr_{3}:Ce Energy (keV)")
canvas_2d.Update()
canvas_2d.SaveAs(f"{save_dir}energy_energy_matrix__delay150us_500nswindow.root")

# plot the time difference between events with the same energy
same_energy = np.array(same_energy)
canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
hist = ROOT.TH2D("hist", "Time difference between events with the same energy", 130, 0, 1300, 10,0,10)
for energy, time in same_energy:
    hist.Fill(energy, time)
hist.Draw("COLZ")
hist.GetXaxis().SetTitle("Energy (keV)")
hist.GetYaxis().SetTitle("Time difference (ns)")
canvas.Update()
canvas.SaveAs(f"{save_dir}sameenergy_time_matrix_delay150us_500nswindow.root")




""" df.loc[df['detectorID'] == 5, 'time']   # 125 microseconds
data_df.sort_values(by='time', ascending=True, ignore_index=True, inplace=True)
data_df.reset_index(inplace=True, drop=True)

print("Ratio of detector 0 to detector 5 events: ", len(data_df[data_df['detectorID'] == 0])/len(data_df[data_df['detectorID'] == 5]))

#subset_size = int(len(data_df) * 0.8)
#data_df = data_df.iloc[:subset_size]

times = data_df['time'].values
energies = data_df['energy'].values
detectorIDs = data_df['detectorID'].values

detector5_indices = np.where(detectorIDs == 5)[0]
detector0_indices = np.where(detectorIDs == 0)[0]

# Split the detector5_indices into chunks for parallel processing
num_chunks = 10  # Number of threads
chunks = np.array_split(detector5_indices, num_chunks)
nbins = 30
xmin = 0
xmax = 300000 # in nanoseconds



# List to store the filenames of the PNG images
png_filenames = []

for coincidence_window in range(5000, 600001, 5000):
    nbins = int(coincidence_window / 50)
    xmin = 0
    xmax = int(coincidence_window)

    with ThreadPoolExecutor(max_workers=num_chunks) as executor:
        futures = [executor.submit(find_backward_events_chunk, chunk, data_df, coincidence_window, energy_lower, energy_upper) for chunk in chunks]
        results = [future.result() for future in futures]

    backward_events = [event for result in results for event in result]

    backward_events_df = pd.DataFrame(backward_events, columns=['timediff', 'labr3energy', 'polarisenergy'])

    time_diffs_backward = [event[0] for event in backward_events]
    labr3energy_backward = [event[1] for event in backward_events]
    polarisenergy_backward = [event[2] for event in backward_events]

    canvas_2d = ROOT.TCanvas("canvas_2d", "canvas_2d", 800, 600)
    hist_2d = ROOT.TH2D("hist_2d", "Coincidence window: {} #mus".format(coincidence_window/1000), 130, 0, 1300, 130, 0, 1300)
    for labr3energy, polarisenergy in zip(labr3energy_backward, polarisenergy_backward):
        hist_2d.Fill(polarisenergy, labr3energy)
    hist_2d.Draw("COLZ")
    hist_2d.GetXaxis().SetTitle("Polaris energy (keV)")
    hist_2d.GetYaxis().SetTitle("LaBr_{3}:Ce energy (keV)")
    canvas_2d.Update()

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    hist = ROOT.TH1D("hist", "LaBr3 event precedes a POLARIS event", nbins, xmin, xmax)
    for timediff in time_diffs_backward:
        hist.Fill(timediff)
    hist.Draw()
    hist.GetXaxis().SetTitle("Time difference (ns)")
    hist.GetYaxis().SetTitle("Counts (1/10 #mus)")
    canvas.Update()
    canvas.SaveAs(f"{save_dir}coinc_time_labrbefore_histogram_{coincidence_window}.root") 

    backward_events_df['timediff'] = backward_events_df['timediff'].fillna(1000000000)

    canvas_3d = ROOT.TCanvas("canvas_3d", "canvas_3d", 800, 600)
    hist_3d = ROOT.TH3D("hist_3d", "LaBr3 event precedes a POLARIS event", 160, 0, 1600, 160, 0, 1600, 100, 0, 1000000)
    for labr3energy, polarisenergy, time_diff in zip(labr3energy_backward, polarisenergy_backward, time_diffs_backward):
        hist_3d.Fill(polarisenergy, labr3energy, time_diff)
    hist_3d.Draw("BOX")
    hist_3d.GetXaxis().SetTitle("POLARIS Energy")
    hist_3d.GetYaxis().SetTitle("LaBr3 Energy")
    hist_3d.GetZaxis().SetTitle("Time Difference (ns)")
    canvas_3d.Update()
    canvas_3d.SaveAs(f"{save_dir}coinc_time_labrbefore_3d_histogram_{coincidence_window}.root")


    png_filename = os.path.join(png_dir, f"coinc_time_labrbefore_2d_histogram_{coincidence_window}.png")
    canvas_2d.SaveAs(png_filename)
    png_filenames.append(png_filename)
 """


#_______________________________________________________________________________________

