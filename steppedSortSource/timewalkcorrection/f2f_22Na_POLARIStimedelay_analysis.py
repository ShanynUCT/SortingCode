################################################################################
#                               to run:                                        #
# e.g. timesyncLaBrPOLARIS % ~/miniconda/envs/my_root_env/bin/python f2f_22Na_POLARIStimedelay_analysis.py ~/Documents/PhD/exp/2024/timing_shaping_characterisation_june2024/CC_polaris_2inch/run7/ (path to parquet file)
__author__ = "Shanyn Hart"
__date__ = "2022-08-24"
__version__ = "1.0"
                                            #
################################################################################
#----
# PYTHON IMPORT STATEMENTS
#----
import ROOT as root
import numpy as np
import pandas as pd
import sys, os
import re
import time
from concurrent.futures import ThreadPoolExecutor

# ----
# ********************* DEFINE FUNCTIONS ******************************
# ----
def data_run_duration(data_df):
    # calculate the run duration
    run_duration = data_df['time'].max() - data_df.loc[data_df['time'] != 0, 'time'].min()
    # convert the run duration to seconds
    run_duration = run_duration * 1e-9 / 60
    print('Run duration: ' + str(run_duration) + ' min\n')
    return run_duration

def find_coincidence_events_chunk(chunk_indices, data_df, time_window, energy_lower, energy_upper):
    coincidence_events = []
    for index_det5 in chunk_indices:
        for index_det0 in range(index_det5 - 1, -1, -1):
            time_difference = data_df.loc[index_det5, 'time'] - data_df.loc[index_det0, 'time']
            if time_difference > time_window:
                break
            if data_df.loc[index_det0, 'detectorID'] == 0 and energy_lower <= data_df.loc[index_det0, 'energy'] <= energy_upper:
                coincidence_events.append((time_difference, data_df.loc[index_det0, 'energy'], data_df.loc[index_det5, 'energy']))
    return coincidence_events

def coincidence(data_df, save_dir, run_num):
    coincidence_window = 1000e3  # 1000 us = 1 ms

    # Define the energy range for 511 keV
    energy_tolerance = 0.15 * 511  # 15% tolerance
    energy_lower = 511 - energy_tolerance
    energy_upper = 511 + energy_tolerance

    # Ensure 'energy' column exists
    if 'energy' not in data_df.columns:
        raise ValueError("The 'energy' column is missing from the DataFrame")

    data_df = data_df.loc[(data_df['energy'] >= energy_lower) & (data_df['energy'] <= energy_upper)]
    data_df.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
    data_df.reset_index(inplace=True, drop=True)

    print('Dataframe after energy range cut:\n', data_df)

    indices_det5 = data_df.loc[data_df['detectorID'] == 5].index
    indices_det0 = data_df.loc[data_df['detectorID'] == 0].index

    print('Number of events in POLARIS: ' + str(len(indices_det5)))
    print('Number of events in LaBr3 det 0: ' + str(len(indices_det0)))
    print('Event ratio labr3/polaris:', str(len(indices_det0) / len(indices_det5)))

    num_chunks = 10
    chunk_size = len(indices_det5) // num_chunks
    chunks = [indices_det5[i:i + chunk_size] for i in range(0, len(indices_det5), chunk_size)]

    with ThreadPoolExecutor(max_workers=num_chunks) as executor:
        futures = [executor.submit(find_coincidence_events_chunk, chunk, data_df, coincidence_window, energy_lower, energy_upper) for chunk in chunks]
        results = [future.result() for future in futures]

    coincidence_events = [event for result in results for event in result]

    coincidence_df = pd.DataFrame(coincidence_events, columns=['time_difference', 'energy_det0', 'energy_det5'])

    return coincidence_df, coincidence_window

def plot_energy_histograms(data_df, save_dir, run_num):
    canvas = root.TCanvas("canvas", "canvas", 800, 600)

    hist_det0 = root.TH1D("hist_det0", "LaBr3:Ce det 0 Energy Spectrum", 1700, 0, 1700)
    for energy in data_df.loc[data_df['detectorID'] == 0, 'energy']:
        hist_det0.Fill(energy)
    hist_det0.Draw()
    hist_det0.GetXaxis().SetTitle("Energy (keV)")
    hist_det0.GetYaxis().SetTitle("Counts (1/keV)")
    hist_det5 = root.TH1D("hist_det5", "Polaris Energy Spectrum", 1700, 0, 1700)
    for energy in data_df.loc[data_df['detectorID'] == 5, 'energy']:
        hist_det5.Fill(energy)
    hist_det5.Draw("same")
    hist_det5.GetXaxis().SetTitle("Energy (keV)")
    hist_det5.GetYaxis().SetTitle("Counts (1/keV)")
    L = root.TLegend(0.6, 0.6, 0.9, 0.9)
    L.AddEntry(hist_det0, "LaBr_{3}:Ce L0", "l")
    L.AddEntry(hist_det5, "Polaris", "l")
    L.Draw()
    canvas.Update()
    canvas.SaveAs(save_dir + '22Na_energy_spectrum_Polaris_L0_run' + str(run_num) + '.root')

def plot_time_difference_histogram(coincidence_df, save_dir, run_num, coincidence_window):
    canvas = root.TCanvas("canvas", "canvas", 800, 600)

    hist_time_diff = root.TH1D("hist_time_diff", "Time Difference Histogram", 300, 0, coincidence_window)
    for time_diff in coincidence_df['time_difference']:
        hist_time_diff.Fill(time_diff)
    hist_time_diff.Draw()
    hist_time_diff.GetXaxis().SetTitle("Time Difference (ns)")
    hist_time_diff.GetYaxis().SetTitle("Counts")
    canvas.SaveAs(save_dir + '22Na_time_difference_histogram_run' + str(run_num) + '.root')

def main():
    start_time = time.time()

    # Read in file
    data_dir = sys.argv[1]
    run_num = str(re.search('run(\d+)', data_dir).group(1))
    data_df = pd.read_parquet(data_dir + 'run' + str(run_num) + '_mergedLaBrPOLARISdata.parquet')
    run_num = int(run_num)
    data_df = pd.DataFrame(data_df)
    data_df.columns = ['detectorID', 'energy', 'time', 'x', 'y', 'z']
    data_df.drop(['x', 'y', 'z'], axis=1, inplace=True)
    data_df.sort_values(by='time', inplace=True, ignore_index=True, ascending=True)
    data_df.reset_index(inplace=True, drop=True)

    save_dir = data_dir + '/plots/' + '22Na_f2f_coinc/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        print('Creating directory: {}'.format(save_dir))

    # Call functions
    # data_run_duration(data_df)
    coincidence_df, coincidence_window = coincidence(data_df, save_dir, run_num)

    plot_energy_histograms(data_df, save_dir, run_num)
    plot_time_difference_histogram(coincidence_df, save_dir, run_num, coincidence_window)

    end_time = time.time()
    print('Time to read in file: ' + str((end_time - start_time) / 60) + ' minutes')

if __name__ == '__main__':
    main()
