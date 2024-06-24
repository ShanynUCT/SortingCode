import pandas as pd
import sys
import re
import fastparquet
import numpy as np
import matplotlib.pyplot as plt
import ROOT as ROOT
import seaborn as sns


plt.rcParams.update({'font.size': 18})
# make axis labels bold
plt.rcParams["axes.labelweight"] = "bold"
# label size
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)

#____________________________________________________________
# eg. timesyncLaBrPOLARIS % ~/miniconda/envs/my_root_env/bin/python f2f_coincidence.py ~/Documents/PhD/exp/2024/timing_shaping_characterisation_june2024/CC_polaris_2inch/run11/ (path to parquet file)
#____________________________________________________________


data_dir = sys.argv[1]
run_num = int(re.search('run(\d+)', data_dir).group(1))
df = pd.read_parquet(data_dir + 'run' + str(run_num) + '_mergedLaBrPOLARISdata.parquet')

energy_tolerance = 0.15 * 511  # keV

# Define a function to find forward events starting at each detectorID == 5 index
def find_forward_events(df, time_limit, detector5_indices, detector0_indices):
    forward_events = []

    for idx in detector0_indices:
        labr3time = df.loc[idx, 'time']
        labr3energy = df.loc[idx, 'energy']
        for idy in detector5_indices:
            if idy > idx:
                polaristime = df.loc[idy, 'time']
                polarisenergy = df.loc[idx, 'energy']
                time_diff = polaristime - labr3time
                if time_diff < time_limit:
                    forward_events.append((time_diff, labr3energy, polarisenergy))
                else:
                    break

    return forward_events

#_____

# Filter the data to retain only detectorID 0 and 5
df_filtered = df[df['detectorID'].isin([0, 5])]
#df_filtered = df_filtered[(df_filtered['energy'] >= (511 - energy_tolerance)) & (df_filtered['energy'] <= (511 + energy_tolerance))]

df_filtered = df_filtered.sort_values(by='time', ascending=True)
df_filtered.reset_index(inplace=True, drop=True)
# Find forward events
time_limit = 500e3 # 600 us
# List of indices where detectorID == 5 and 0
detector5_indices = df_filtered.index[df_filtered['detectorID'] == 5].tolist()
detector0_indices = df_filtered.index[df_filtered['detectorID'] == 0].tolist()

forward_events = find_forward_events(df_filtered, time_limit, detector5_indices, detector0_indices)

# Convert to DataFrame for analysis
forward_events_df = pd.DataFrame(forward_events, columns=['timediff', 'labr3energy', 'polarisenergy'])
mean_time_difference = np.mean(forward_events['timediff'])
median_time_difference = np.median(forward_events['timediff'])
time_resolution = np.std(forward_events['timediff'])
print("\nTime Resolution std dev:", time_resolution)
print("\nTime Resolution mean:", mean_time_difference)
print("\nTime Resolution median:", median_time_difference)

# Plot the histogram of the time resolution using Seaborn
plt.figure(figsize=(10, 6))
sns.histplot(data=forward_events_df, x='timediff', bins=1000, kde=False)
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts')
plt.title('Time Resolution Histogram')
plt.grid(True)
plt.show()

# Plot the time differences as a scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(np.arange(len(forward_events)), forward_events['timediff'], marker='.', color='b')
plt.axhline(np.mean(forward_events['timediff']), color='r', linestyle='--', label='Mean')
plt.axhline(np.median(forward_events['timediff']), color='g', linestyle='--', label='Median')
plt.xlabel('Event Index')
plt.ylabel('Time Difference (ns)')
plt.title('Time Differences between detectorID == 5 and detectorID == 0 Events')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Extract labr3energy and polarisenergy from forward_events
labr3energy = [event[1] for event in forward_events]
polarisenergy = [event[2] for event in forward_events]

# Plot 2D histogram
plt.figure(figsize=(8, 6))
plt.hist2d(labr3energy, polarisenergy, bins=30, cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel('Labr3 Energy')
plt.ylabel('Polaris Energy')
plt.title('2D Histogram of Labr3 Energy vs Polaris Energy')
plt.grid(True)
plt.show()


#____

# Adjust time values for detectorID == 5 to account for the offset
df_filtered.loc[df['detectorID'] == 5, 'time'] -= mean_time_difference
df_filtered = df.sort_values(by='time', ascending=True)
df_filtered.reset_index(inplace=True, drop=True)
# List of indices where detectorID == 5 and 0
detector5_indices = df_filtered.index[df['detectorID'] == 5].tolist()
detector0_indices = df_filtered.index[df['detectorID'] == 0].tolist()
#time_limit = 100e3 # 100us

# Recalculate time differences with adjusted times
forward_events_with_correction = find_forward_events(df_filtered, time_limit, detector5_indices, detector0_indices)


forward_events_df_with_correction = pd.DataFrame(forward_events_with_correction, columns=['timediff'])
time_resolution_with_correction = np.median(forward_events_with_correction)
print("Time Resolution with Correction median:", time_resolution_with_correction)


plt.figure(figsize=(10, 6))
sns.histplot(data=forward_events_df_with_correction, x='timediff', bins=1000, kde=False)
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts')
plt.title('Time Resolution Histogram with Correction')
plt.grid(True)
plt.show()


#_____


def find_coincident_events(df, time_window_ns):
    coincident_events = []
    for index, row in df.iterrows():
        # Find all events within the time window
        time_diff = (df['time'] - row['time']).abs()
        coincident = df[(time_diff <= time_window_ns) & (df['detectorID'] != row['detectorID'])]
        
        for _, coin_event in coincident.iterrows():
            coincident_events.append((row, coin_event))
    
    return coincident_events

# Define the time window in nanoseconds
time_window_ns = 300  # Adjust as necessary

# Find coincident events
coincident_events = find_coincident_events(df_filtered, time_window_ns)

# Calculate time differences
time_differences = [abs(e1['time'] - e2['time']) for e1, e2 in coincident_events]

# Calculate time resolution (e.g., FWHM or standard deviation)
time_resolution = np.std(time_differences)
print("Time Resolution (std dev):", time_resolution)

# Plot the stepped histogram of the time resolution
plt.figure(figsize=(10, 6))
plt.hist(time_differences, bins=100, histtype='step', linewidth=1.5)
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts')
plt.title('Time Resolution Histogram')
plt.grid(True)
plt.show()


# Extract energies
energies1 = [e1['energy'] for e1, _ in filtered_coincident_events]
energies2 = [e2['energy'] for _, e2 in filtered_coincident_events]

# Plot the energy-energy matrix
plt.hist2d(energies1, energies2, bins=100, cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel('Energy of Detector 0 (keV)')
plt.ylabel('Energy of Detector 5 (keV)')
plt.title('Energy-Energy Matrix')
plt.show()