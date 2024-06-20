import sys
import uproot
import numpy as np
import pandas as pd

inputFileName = sys.argv[1]

file = uproot.open(inputFileName)
tree = file["LaBrData"]

timeFast0 = tree["timeFL0"].array()
timeFast0 = timeFast0[timeFast0 != 0]
timeFast0 = timeFast0/10
timeFast0 = pd.DataFrame(timeFast0, columns = ['time'])
# sort ascending by time
#timeFast0 = timeFast0.drop_duplicates(subset = ['time'])

timeFast0 = timeFast0.sort_values(by = ['time'], ascending = True, ignore_index = True)
DeltaFast0 = np.diff(timeFast0['time'])

print("time duration of measurement fast: " + str((timeFast0['time'].iloc[-1] - timeFast0['time'].iloc[0])*1e-9/60))

timeSlow0 = tree["timeSL0"].array()
timeSlow0 = timeSlow0/10
energySlow0 = tree["slowECalibL0"].array()
data = {'time': timeSlow0, 'energy': energySlow0}
timeSlow0 = pd.DataFrame(data)
timeSlow0 = timeSlow0.sort_values(by = ['time'], ascending = True, ignore_index = True)
# select rows where energy is not 0 and time is not 0
timeSlow0 = timeSlow0[ (timeSlow0['energy'] != 0) & (timeSlow0['time'] != 0)]
# keep only unique values
timeSlow0 = timeSlow0.drop_duplicates(subset = ['time'])
DeltaSlow0 = np.diff(timeSlow0['time'])

print("time duration of measurement slow: " + str((timeSlow0['time'].iloc[-1] - timeSlow0['time'].iloc[0])*1e-9/60))

data_dir = sys.argv[1].replace(sys.argv[1].split('/')[-1], '')
data_filename = 'mod51.txt'

mod_event_data = pd.read_csv(data_dir + data_filename, sep = '	', header = None)  
mod_event_data.columns = ['scatters', 'x', 'y', 'z', 'energy', 'time']   
mod_event_data = mod_event_data[ (mod_event_data['scatters'] == 1)]
mod_event_data_df = pd.DataFrame(mod_event_data)
mod_event_data_df['time'] = mod_event_data_df['time'] * 10
mod_event_data_df = mod_event_data_df[mod_event_data_df['time'] != 0]
# keep only the first 8.54 minutes of data ( find the first time value and add 8.54 minutes to it)
mod_event_data_df = mod_event_data_df[mod_event_data_df['time'] < (mod_event_data_df['time'].iloc[0] + 8.54*60*1e9)]
# sort ascending by time
mod_event_data_df = mod_event_data_df.sort_values(by = ['time'], ascending = True, ignore_index = True)
DeltaMod = np.diff(mod_event_data_df['time'])

print("time duration of measurement mod: " + str((mod_event_data_df['time'].iloc[-1] - mod_event_data_df['time'].iloc[0])*1e-9/60))

print("Mean time difference between consecutive events POLARIS:" + str(np.mean(DeltaMod)))
print("Mean time difference between consecutive events FAST:" + str(np.mean(DeltaFast0)))
print("Mean time difference between consecutive events SLOW:" + str(np.mean(DeltaSlow0)))

# print the number of events in each histogram
print("Number of events in mod: " + str(len(DeltaMod)))
print("Number of events in fast: " + str(len(DeltaFast0)))
print("Number of events in slow: " + str(len(DeltaSlow0)))

# # Plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches

#plot all three histograms on the same axes
plt.figure(1)
# step histogram
plt.hist(DeltaMod, bins = 4000, range = (0, 40000),  color = 'blue', label = 'Polaris', histtype = 'step', linewidth = 2)
plt.hist(DeltaFast0, bins = 4000, range = (0, 40000), color = 'red', label = 'LaBr$_3$(Ce) Fast Signal', histtype = 'step', linewidth = 2)
plt.hist(DeltaSlow0, bins = 4000, range = (0, 40000), color = 'green', label = 'LaBr$_3$(Ce) Slow Signal', histtype = 'step', linewidth = 2)
# do the labels in latex font
plt.rc('font', family = 'serif')
plt.xlabel('Time difference between consecutive events (ns)', fontsize = 20)
plt.ylabel('Counts (1/10 ns)', fontsize = 20)
plt.legend(loc = 'upper right', fontsize = 20)
plt.tick_params(axis = 'both', which = 'major', labelsize = 20)
# increase number of xticks
plt.locator_params(axis = 'x', nbins = 20)
plt.xlim(0, 30000)
plt.yscale('log')
plt.show()

