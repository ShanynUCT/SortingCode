import pandas as pd
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ROOT
import os
from concurrent.futures import ThreadPoolExecutor

plt.rcParams.update({'font.size': 18})
# Make axis labels bold
plt.rcParams["axes.labelweight"] = "bold"
# Label size
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)

#____________________________________________________________
# e.g. timesyncLaBrPOLARIS % ~/miniconda/envs/my_root_env/bin/python f2f_coincidence.py ~/Documents/PhD/exp/2024/timing_shaping_characterisation_june2024/CC_polaris_2inch/run11/ (path to parquet file)
#____________________________________________________________

data_dir = sys.argv[1]
save_dir = data_dir + "/plots/22Na_f2f_coinc/"

# Check if save_dir exists, if not, create it
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

run_num = int(re.search('run(\d+)', data_dir).group(1))
df = pd.read_parquet(data_dir + 'run' + str(run_num) + '_mergedLaBrPOLARISdata.parquet')

# Sort the DataFrame by the 'time' column in ascending order
df = df.sort_values(by='time', ascending=True)

# Filter the data to retain only detectorID 0 and 5 and energies close to 511 keV
energy_tolerance = 0.15 * 511  # 15% tolerance
#df_filtered = df[(df['detectorID'].isin([0, 5])) & ((np.abs(df['energy'] - 511) <= energy_tolerance) | (np.abs(df['energy'] + 511) <= energy_tolerance))]
df_filtered = df[(df['detectorID'].isin([0, 5])) ]
df_filtered.reset_index(inplace=True, drop=True)

# Use only the first 35% of the DataFrame for performance reasons
subset_size = int(len(df_filtered) * 0.8)
df_filtered = df_filtered.iloc[:subset_size]

# Extract relevant columns and convert to numpy arrays for speed
times = df_filtered['time'].values
energies = df_filtered['energy'].values
detectorIDs = df_filtered['detectorID'].values

# Get indices for detectorID 5 and 0
detector5_indices = np.where(detectorIDs == 5)[0]
detector0_indices = np.where(detectorIDs == 0)[0]

# Function to find backward events where LaBr3 event precedes a POLARIS event
def find_backward_events_chunk(chunk):
    backward_events = []
    time_window = 600e3  # 500 microseconds in nanoseconds

    for idy in chunk:
        polaristime = times[idy]
        polarisenergy = energies[idy]

        # Find the index in detector0_indices where the labr3time is just less than the polaristime using binary search
        idx = np.searchsorted(times[detector0_indices], polaristime) - 1
        
        while idx >= 0 and times[detector0_indices[idx]] < polaristime:
            labr3time = times[detector0_indices[idx]]
            labr3energy = energies[detector0_indices[idx]]
            time_diff = polaristime - labr3time
            
            if 0 < time_diff <= time_window:
                backward_events.append((time_diff, labr3energy, polarisenergy))
                break
            idx -= 1

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

# Split the detector5_indices into chunks for parallel processing
num_chunks = 10  # Number of threads
chunks = np.array_split(detector5_indices, num_chunks)
nbins = 600
xmin = 0
xmax = 600e3  # in nanoseconds


""" with ThreadPoolExecutor(max_workers=num_chunks) as executor:
    backwards = [executor.submit(find_backward_events_chunk, chunk) for chunk in chunks]
    results = [backward.result() for backward in backwards]

backward_events = [event for result in results for event in result]

backward_events_df = pd.DataFrame(backward_events, columns=['timediff', 'labr3energy', 'polarisenergy'])
time_diffs_backward = [event[0] for event in backward_events]
labr3energy_backward = [event[1] for event in backward_events]
polarisenergy_backward = [event[2] for event in backward_events]


plt.figure(figsize=(10, 6))
sns.histplot(data=backward_events_df, x='timediff', bins=np.arange(0, 300000 + 1000, 1000), kde=False)
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts (1/\u00B5s)')
plt.title('LaBr3 event precedes a POLARIS event')
plt.grid(True)
plt.show() 

# Save the histogram with ROOT
canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
hist = ROOT.TH1D("hist", "LaBr3 event precedes a POLARIS event", nbins, xmin, xmax)
for timediff in time_diffs_backward:
    hist.Fill(timediff)
hist.Draw()
hist.GetXaxis().SetTitle("Time Difference (ns)")
hist.GetYaxis().SetTitle("Counts (1/#mu s)")
canvas.Update()
canvas.SaveAs(save_dir + "coinc_time_labrbefore_histogram.root") 

# 2D Histogram for energy vs time difference
canvas_2d = ROOT.TCanvas("canvas_2d", "canvas_2d", 800, 600)
hist_2d = ROOT.TH2D("hist_2d", "LaBr3 event precedes a POLARIS event", 1600, 0, 1600, 1600, 0, 1600)
for labr3energy, polarisenergy in zip(labr3energy_backward, polarisenergy_backward):
    hist_2d.Fill(polarisenergy, labr3energy)
hist_2d.Draw("COLZ")
hist_2d.GetXaxis().SetTitle("POLARIS Energy")
hist_2d.GetYaxis().SetTitle("LaBr3 Energy")
canvas_2d.Update()
canvas_2d.SaveAs(save_dir + "coinc_time_labrbefore_2d_histogram.root")  """
#____________________________________________________________________________
 
# Use ThreadPoolExecutor to process chunks in parallel
""" with ThreadPoolExecutor(max_workers=num_chunks) as executor:
    forwards = [executor.submit(find_forward_events_chunk, chunk) for chunk in chunks]
    results = [forward.result() for forward in forwards]

# Combine results from all threads
forward_events = [event for result in results for event in result]

# Convert to DataFrame for analysis
forward_events_df = pd.DataFrame(forward_events, columns=['timediff', 'labr3energy', 'polarisenergy'])
time_diffs_forward= [event[0] for event in forward_events]
labr3energy_forward= [event[1] for event in forward_events]
polarisenergy_forward= [event[2] for event in forward_events]

# Plot the histogram of the time resolution using Seaborn with 1 microsecond bins (1000 nanoseconds)
plt.figure(figsize=(10, 6))
sns.histplot(data=forward_events_df, x='timediff', bins=np.arange(0, 300000 + 1000, 1000), kde=False)
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts (1/\u00B5s)')
plt.title('LaBr3 event follows a POLARIS event')
plt.grid(True)
plt.show()

canvas1 = ROOT.TCanvas("canvas1", "canvas1", 800, 600)
hist1 = ROOT.TH1D("hist1", "LaBr3 event follows a POLARIS event", nbins, xmin, xmax)
for timediff in time_diffs_forward:
    hist1.Fill(timediff)
hist1.Draw()
hist1.GetXaxis().SetTitle("Time Difference (ns)")
hist1.GetYaxis().SetTitle("Counts (1/#mu s)")
canvas1.Update()
canvas1.SaveAs(save_dir + "coinc_time_labrafter_histogram.root") 

# 2D Histogram for energy vs time difference
canvas_2d = ROOT.TCanvas("canvas_2d", "canvas_2d", 800, 600)
hist_2d = ROOT.TH2D("hist_2d", "LaBr3 event precedes a POLARIS event", 1600, 0, 1600, 1600, 0, 1600)
for labr3energy, polarisenergy in zip(labr3energy_forward, polarisenergy_forward):
    hist_2d.Fill(polarisenergy, labr3energy)
hist_2d.Draw("COLZ")
hist_2d.GetXaxis().SetTitle("POLARIS Energy")
hist_2d.GetYaxis().SetTitle("LaBr3 Energy")
canvas_2d.Update()
canvas_2d.SaveAs(save_dir + "coinc_time_labrafter_2d_histogram.root")  """


#   __________________________________________________
energy_lower = 511 - energy_tolerance
energy_upper = 511 + energy_tolerance

def find_energy_matched_events_chunk(chunk_indices, times, energies, detector0_indices, time_window=500e3):
    forward_events = []
    for idy in chunk_indices:
        polaristime = times[idy]
        polarisenergy = energies[idy]
        if not (energy_lower <= polarisenergy <= energy_upper):
            continue

        # Find the index in detector0_indices where the labr3time is just less than the polaristime using binary search
        idx = np.searchsorted(times[detector0_indices], polaristime) - 1

        while idx >= 0 and times[detector0_indices[idx]] < polaristime:
            labr3time = times[detector0_indices[idx]]
            labr3energy = energies[detector0_indices[idx]]
            time_diff = polaristime - labr3time

            if 0 < time_diff <= time_window and (energy_lower <= labr3energy <= energy_upper):
                forward_events.append((time_diff, labr3energy, polarisenergy))
                break

            idx -= 1

    return forward_events

# Number of chunks for multithreading
chunk_size = len(detector5_indices) // num_chunks
chunks = [detector5_indices[i:i + chunk_size] for i in range(0, len(detector5_indices), chunk_size)]

with ThreadPoolExecutor(max_workers=num_chunks) as executor:
    futures = [executor.submit(find_energy_matched_events_chunk, chunk, times, energies, detector0_indices) for chunk in chunks]
    results = [future.result() for future in futures]

forward_events = [event for result in results for event in result]

forward_events_df = pd.DataFrame(forward_events, columns=['timediff', 'labr3energy', 'polarisenergy'])
time_diffs = [event[0] for event in forward_events]
labr3energy = [event[1] for event in forward_events]
polarisenergy = [event[2] for event in forward_events]

# Define the bin edges for the histogram
x_edges = np.linspace(0, 1300, 1300)
y_edges = np.linspace(0, 1300, 1300)

# Create a 2D histogram without weights (only counts)
hist, x_edges, y_edges = np.histogram2d(labr3energy, polarisenergy, bins=[x_edges, y_edges])

# Plot the histogram
plt.figure(figsize=(8, 6))
plt.imshow(hist.T, origin='lower', cmap='viridis', aspect='auto', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]])
plt.colorbar(label='Counts')
plt.xlabel('LaBr3 Energy')
plt.ylabel('POLARIS Energy')
plt.title('2D Histogram of LaBr3 Energy vs POLARIS Energy')
plt.grid(True)
plt.show()

# Save the 2D histogram with ROOT
canvas_2d = ROOT.TCanvas("canvas_2d", "canvas_2d", 800, 600)
hist_2d = ROOT.TH2D("hist_2d", "LaBr3 Energy vs POLARIS Energy", 1300, 0, 1300, 1300, 0, 1300)
for labr3_e, polaris_e in zip(labr3energy, polarisenergy):
    hist_2d.Fill(labr3_e, polaris_e)
hist_2d.Draw("COLZ")
hist_2d.GetXaxis().SetTitle("LaBr3 Energy")
hist_2d.GetYaxis().SetTitle("POLARIS Energy")
canvas_2d.Update()
canvas_2d.SaveAs(save_dir + "coinc_energy_energy_2d_histogram.root")

# Plot the 1D histogram of time differences
plt.figure(figsize=(10, 6))
sns.histplot(time_diffs, bins=np.arange(0, 300000 + 1000, 1000), kde=False)
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts (1/μs)')
plt.title('Time Difference Histogram')
plt.grid(True)
plt.show()

# Save the 1D histogram with ROOT
canvas_1d = ROOT.TCanvas("canvas_1d", "canvas_1d", 800, 600)
hist_1d = ROOT.TH1D("hist_1d", "Time Difference Histogram", 300, 0, 300000)
for timediff in time_diffs:
    hist_1d.Fill(timediff)
hist_1d.Draw()
hist_1d.GetXaxis().SetTitle("Time Difference (ns)")
hist_1d.GetYaxis().SetTitle("Counts (1/μs)")
canvas_1d.Update()
canvas_1d.SaveAs(save_dir + "time_difference_histogram.root")