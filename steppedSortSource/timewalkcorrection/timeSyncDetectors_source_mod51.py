################################################################################
#                               to run:                                        #
# python3 timeSyncDetectors.py /directory/to/rootfile/  (R**_rawData.root)   
# eg. timesyncLaBrPOLARIS % ~/miniconda/envs/my_root_env/bin/python timeSyncDetectors_source_mod51.py ~/Documents/PhD/exp/2024/timing_shaping_characterisation_june2024/CC_polaris_2inch/run11/ R11_aligned.root(obtained from running sort1labr.C and then sort2labr.C)
__author__ = "Shanyn Hart"
__date__ = "2023-05-26"
__version__ = "2.0"

# adapted from the other file with the same name in /timesyncLaBrPOLARIS folder 

######################## CODE UPDATES: ########################
# 2022-08-22: Created script to run all detectors for all runs
# 2022-08-22: Added in all relevant transformations for each run
# 2022-08-22: Added in all relevant detectors for each run
# 2023-03-30: Changed LaBr3 TTree branch names to correspond to new RXX.root file
#             which has been updated to only use insync events with RF for background reduction.
# 2023-05-26: All LaBr3:Ce detectors share a tree. Code updated to reflect this.
# 2023-05-26: Code improved to run with functions                                
################################################################################

#------------------------------------------------------------------
# python import statements
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
import csv as csv
import io as io
import sys, os
from math import sin, cos, pi, log, floor
import cProfile
import re
import pylab
import pyroot as pr
import decimal
import time as time
import fastparquet
import pyarrow
#os.environ["OPENBLAS_NUM_THREADS"] = "200"
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()}, language_level=3)
run_dir = os.getcwd()
LaBrPOLARIS_utils =  run_dir + "/LaBrPOLARIS_utils.pyx"
import LaBrPOLARIS_utils as utils #type: ignore
from ydata_profiling import ProfileReport

# run this script in multiple threads
os.environ["NUMEXPR_MAX_THREADS"] = "16"

pd.options.mode.chained_assignment = None  # default='warn'

#_______________________________________________________________________________________________
# ________________  RELEVANT TRANSFORMATIONS FOR D1 FOR EACH RUN 28/02/2023 ____________________ 
#_______________________________________________________________________________________________
#NOT DONE YET
transformation_run1 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run2 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run3 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run4 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run5 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run6 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run7 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run8 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]
transformation_run9 = [90.0, 0.0, 0.0, -11.0, 192.0, 0.0]
transformation_run10 = [90.0, 0.0, 0.0, -11.0, 192.0, 0.0]
transformation_run11 = [90.0, 0.0, 0.0, -11.0, 192.0, 0.0]
transformation_run12 = [0.0, 0.0, 0.0, -11.0, 0.0, -192.0]
transformation_run13 = [0.0, 0.0, 0.0, -11.0, 0.0, -192.0]
transformation_run14 = [0.0, 0.0, 0.0, -11.0, 0.0, -192.0]
transformation_run15 = [0.0, 0.0, 0.0, -11.0, 0.0, -192.0]
transformation_run16 = [45.0, 0.0, 0.0, -11.0, 135.8, -135.8]
transformation_run17 = [45.0, 0.0, 0.0, -11.0, 135.8, -135.8]
transformation_run18 = [45.0, 0.0, 0.0, -11.0, 213.5, -213.5]
transformation_run19 = [45.0, 0.0, 0.0, -11.0, 213.5, -213.5]
transformation_run20 = [90.0, 0.0, 0.0, -11.0, 362.0, 0.0]

#make a dictionary of the transformations (automates the process for each run)
transformation_dict = {1: transformation_run1, 2: transformation_run2, 3: transformation_run3, 4: transformation_run4, 5: transformation_run5, 6: transformation_run6, 7: transformation_run7, 8: transformation_run8, 9: transformation_run9, 10: transformation_run10, 11: transformation_run11, 12: transformation_run12, 13: transformation_run13, 14: transformation_run14, 15: transformation_run15, 16: transformation_run16, 17: transformation_run17, 18: transformation_run18, 19: transformation_run19, 20: transformation_run20}


#_______________________________________________________________________________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________________________________________________________________________
def read_polaris_data(data_dir):
    print('\033[91m' +'________________________________________') 
    print('\033[91m' +'_____________ READ IN POLARIS DATA _____________')
    print('\033[91m' +'________________________________________', '\n')

    data_filename = 'mod51.txt'
    run_num = int(re.search('run(\d+)', data_dir).group(1))

    start_time = time.time()

    print('-------------->> Run number: ', run_num)

    mod_event_data = pd.read_csv(data_dir + data_filename, sep='	', header=None)
    mod_event_data.columns = ['scatters', 'x', 'y', 'z', 'energy', 'time']
    mod_event_data = mod_event_data[(mod_event_data['scatters'] == 122) | (mod_event_data['scatters'] == 1)]
    mod_event_data_df = pd.DataFrame(mod_event_data)

    evts_all = len(mod_event_data.index)
    time_total = (mod_event_data.iloc[-1]['time']*10 - mod_event_data.iloc[0]['time']*10) * 1e-9
    evt_rate_all = evts_all / time_total

    print('\033[91m' +'{:30} {:<10d}'.format('\n nEvents POLARIS: ', evts_all),
          '{:30} {:<10.1f}'.format('\n Total time of Run POLARIS (minutes): ', time_total/60),
          '{:30} {:<10.1f}'.format('\n Event rate POLARIS (events/sec): ', evt_rate_all), '\n\n')

    if data_filename == 'mod51.txt':
        det_index = 0
    mod_event_data_df['detector'] = det_index
    mod_event_data_df = mod_event_data_df[['detector', 'scatters', 'energy', 'x', 'y', 'z', 'time']]
    mod_event_data_df['scatters'] = mod_event_data_df['scatters'].astype(int)
    mod_event_data_df = mod_event_data_df.reset_index(drop=True)
    mod_event_data_df['time'] = mod_event_data_df['time'].apply(lambda x: int(x) if not np.isnan(x) else x)

    print('\033[91m' +'Dataframe with POLARIS data created. \n', mod_event_data_df)

    return run_num, mod_event_data_df, mod_event_data

def apply_coordinate_transformations(run_num,mod_event_data_df):
    print('\033[92m' +'________________________________________')
    print('\033[92m' +'Applying transformations to convert data into isocentric coordinate system...')
    print('\033[92m' +'________________________________________', '\n')
    transformation_matrices = {}  # Creating transformation matrices
    transformation_matrices[0] = utils.get_transformation_matrix_array(transformation_dict[run_num])
    mod_event_data_df.drop(['detector'], inplace=True, axis=1)
    mod_event_data_df = mod_event_data_df[['scatters', 'energy', 'time', 'x', 'y', 'z']]
    mod_event_data_df['energy'] = mod_event_data_df['energy'] / 1e3  # convert to keV
    mod_event_data_df['x'] = mod_event_data_df['x'] / 1e3  # convert to mm
    mod_event_data_df['y'] = mod_event_data_df['y'] / 1e3  # convert to mm
    mod_event_data_df['z'] = mod_event_data_df['z'] / 1e3  # convert to mm
    mod_event_data_df['time'] = mod_event_data_df['time']*10  # convert to seconds

    print('\033[92m' +'Dataframe with POLARIS data created and coordinates transformed. \n', mod_event_data_df)

    return mod_event_data_df

def find_POLARIS_sync_pulses(mod_event_data, mod_event_data_df, data_dir):
    print('\033[93m' +'________________________________________') 
    print ('\033[93m' +'_____________ FIND SYNC PULSES IN POLARIS DATAFRAME _____________ ')
    print('\033[93m' +'________________________________________', '\n')
    mod_sync_pulses = mod_event_data.loc[mod_event_data['scatters'] == 122]
    mod_sync_pulses.drop(['z', 'energy', 'time', 'detector'], inplace=True, axis=1)
    mod_sync_pulses.columns = ['sync_flag', 'sync_index', 'sync_timestamp']
    mod_sync_pulses.reset_index(drop=True, inplace=True)
    mod_sync_pulses['sync_timestamp'] = mod_sync_pulses['sync_timestamp']*10
    mod_sync_pulses['sync_energy'] = mod_event_data['energy'].iloc[mod_sync_pulses['sync_index'] + 1]
    mod_sync_pulses_df = pd.DataFrame(mod_sync_pulses)
    mod_sync_pulses_df.reset_index(drop=True, inplace=True)
    print('\033[93m' +'____________  Sync Pulses POLARIS before selecting first valid index ____________ \n', mod_sync_pulses_df)
    #mod_sync_pulses_df.to_parquet(data_dir + 'syncPulsesPOLARIS.parquet')

    firstValidSyncPOLARIS = mod_sync_pulses_df.loc[mod_sync_pulses_df['sync_energy'].first_valid_index()]['sync_index']
    firstValidSyncPOLARIS_timestamp = mod_sync_pulses_df.loc[firstValidSyncPOLARIS - 1]['sync_timestamp']
    mod_sync_pulses_df = mod_sync_pulses_df.loc[mod_sync_pulses_df['sync_timestamp'] >= firstValidSyncPOLARIS_timestamp]
    mod_sync_pulses_df.reset_index(drop=True, inplace=True)
    #print('\033[93m' +'____________  Sync Pulses POLARIS  ____________ \n', mod_sync_pulses_df)

    sync_time_diff_POLARIS = (mod_sync_pulses_df['sync_timestamp']).diff()
    #print('\033[93m' +'____________  Time Diff Between POLARIS Sync Pulses  ____________ \n', sync_time_diff_POLARIS)
    mod_sync_pulses_df['sync_time_diff_POLARIS'] = sync_time_diff_POLARIS
    sync_time_diff_POLARIS = sync_time_diff_POLARIS[1:]
    print('\033[93m' +'First Valid Timestamp: ', firstValidSyncPOLARIS_timestamp,
          '\n First Valid Index: ', firstValidSyncPOLARIS,
          ' \n Average POLARIS sync difference: {} s +- {} s'.format(np.mean(mod_sync_pulses_df['sync_time_diff_POLARIS']),
                                                                       np.sqrt(
                                                                           np.var(mod_sync_pulses_df['sync_time_diff_POLARIS']))),
          '\n\n')

    mod_event_data_df = mod_event_data_df.loc[mod_event_data_df['time'] >= firstValidSyncPOLARIS_timestamp]
    mod_event_data_df.reset_index(drop=True, inplace=True)

    print('\033[93m' +'____________ POLARIS df  ____________ \n', mod_event_data_df)
    return mod_sync_pulses_df,firstValidSyncPOLARIS,mod_event_data_df

def read_LaBr_data(data_dir):
    print('\033[94m' + '________________________________________') 
    print('\033[94m' + '_____________ READ IN LaBr3:Ce DATA _____________')
    print('\033[94m' + '________________________________________', '\n')
    

    LaBr_file = sys.argv[2]
    LaBr_data = uproot.open(data_dir + LaBr_file)
    LaBr_tree = LaBr_data["dir_aligned/AlignedData"]
    detectorID = "detectorID"  
    globalTime = "globalTime"
    timeF = "timeF"
    energyF = "energyF"
    timeS = "timeS"
    energyS = "energyS"

    detectorID = LaBr_tree[detectorID].array()
    globalTime = LaBr_tree[globalTime].array()
    timeF = LaBr_tree[timeF].array()
    energyF = LaBr_tree[energyF].array()
    timeS = LaBr_tree[timeS].array()
    energyS = LaBr_tree[energyS].array()
    

    LaBr_data_df = pd.DataFrame({'detectorID': detectorID, 'globalTime': globalTime, 'timeF': timeF, 'energyF': energyF, 'timeS': timeS, 'energyS': energyS})
    LaBr_data_df.reset_index(drop=True, inplace=True)

    # DROP detector 4 (RF)
    LaBr_data_df = LaBr_data_df.loc[LaBr_data_df['detectorID'] != 4]
    #drop negative times
    LaBr_data_df = LaBr_data_df.loc[LaBr_data_df['timeF'] > 0]
    LaBr_data_df = LaBr_data_df.sort_values(by=['globalTime'], ascending=True)
    LaBr_data_df.reset_index(drop=True, inplace=True)
    # take the time difference between the first non zero globalTime and the last non zero globalTime
    first_time = LaBr_data_df.loc[LaBr_data_df['timeF'] != 0.0].iloc[0]['timeF']
    last_time = LaBr_data_df.loc[LaBr_data_df['timeF'] != 0.0].iloc[-1]['timeF']

    time_duration_labr = (last_time - first_time)*1e-9/60 
    # if the time difference is > 50 minutes, then the bool tensOfNanoseconds is True and we divide the times by 10
    if (time_duration_labr >= 50.0):
        LaBr_data_df['globalTime'] = LaBr_data_df['globalTime']/10
        LaBr_data_df['timeF'] = LaBr_data_df['timeF']/10
        LaBr_data_df['timeS'] = LaBr_data_df['timeS']/10
    elif (time_duration_labr >= 500.0):
        LaBr_data_df['globalTime'] = LaBr_data_df['globalTime']/100
        LaBr_data_df['timeF'] = LaBr_data_df['timeF']/100
        LaBr_data_df['timeS'] = LaBr_data_df['timeS']/100
    else:
        pass # do nothing


    print('\nTime duration LaBr3 data: ', time_duration_labr, 'min\n')

    min_globalTime = LaBr_data_df.loc[LaBr_data_df['globalTime'] != 0.0].min()['globalTime']
    max_globalTime = LaBr_data_df['globalTime'].max()

    print('\n Max - Min duration gloal time LaBr3 ', (max_globalTime - min_globalTime)*1e-9/60, 'min\n')

    print('\033[94m' +'____________ LaBr3:Ce df  ____________ \n', LaBr_data_df)


    return LaBr_data_df

def find_LaBr3_sync_pulses(firstValidSyncPOLARIS,LaBr_data_df,mod_sync_pulses_df):
    print('\033[95m' + '__________________________________________________________________')
    print ('\033[95m' + '_____________ FIND SYNC PULSES IN LaBr3:Ce DATAFRAME _____________ ')
    print('\033[95m' + '__________________________________________________________________', '\n')

    sync_pulses_LaBr_df = pd.DataFrame(LaBr_data_df)
    sync_pulses_LaBr_df = sync_pulses_LaBr_df.loc[sync_pulses_LaBr_df['detectorID'] == 5] # POLARIS is detector 5
    sync_pulses_LaBr_df = sync_pulses_LaBr_df.sort_values(by=['timeF'], ascending=True)
    sync_pulses_LaBr_df.reset_index(drop=True, inplace=True)

    print("\n sync pulses df: " ,sync_pulses_LaBr_df,  )
    first_sync_pulse_time = sync_pulses_LaBr_df['timeF'].iloc[0]  # Get the time of the first sync pulse

    # Calculate the time difference for every sync pulse relative to the first
    sync_pulses_LaBr_df['time_difference'] = (sync_pulses_LaBr_df['timeF'] - first_sync_pulse_time)

    # Divide the time difference by 200000006 to get the ratio
    sync_pulses_LaBr_df['time_ratio'] = sync_pulses_LaBr_df['time_difference'] / 2000000060 # 2000000060 and not 200000006 because POLARIS sends pulses in tens of nanoseconds (1e8)
    #sync_pulses_LaBr_df = sync_pulses_LaBr_df.drop_duplicates(subset = ['timeF'], keep = 'first')

    # group 0.0 to <1.0 together, 1.0 to <2.0 together, etc
    sync_pulses_LaBr_df['time_ratio'] = sync_pulses_LaBr_df['time_ratio'].apply(lambda x: floor(x))

    print('\033[95m' +'____________  Time Ratio of LaBr3:Ce Sync Pulses ____________ \n', sync_pulses_LaBr_df['time_ratio'])

    sync_pulses_LaBr_df = sync_pulses_LaBr_df.groupby('time_ratio').first().reset_index()
    sync_pulses_LaBr_df.reset_index(drop=True, inplace=True)

    print('\033[95m' +'____________  Sync Pulses LaBr3:Ce before selecting first valid index ____________ \n', sync_pulses_LaBr_df)


    # plt.figure()
    # plt.plot(sync_pulses_LaBr_df['timeF'], sync_pulses_LaBr_df['energyF'], 'o', markersize=0.5)
    # plt.xlabel('timeF')
    # plt.ylabel('energyF')
    # plt.title('LaBr3:Ce Energy vs Time')
    # plt.text(0.5, 0.5, 'Number of entries: ' + str(len(sync_pulses_LaBr_df)), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
    # plt.show()

    firstValidSyncPOLARIS = int(firstValidSyncPOLARIS)
    sync_pulses_LaBr_df = sync_pulses_LaBr_df.iloc[(firstValidSyncPOLARIS-1):]# for instance index 3 is actually 2 etc
    excess_LaBr_sync_pulses = len(sync_pulses_LaBr_df)-len(mod_sync_pulses_df)  
    sync_pulses_LaBr_df = sync_pulses_LaBr_df.drop(['time_ratio', 'time_difference'], axis=1)

    print('\033[95m' +'____________ Time Duration of LaBr3 Sync Pulses____________\n', (sync_pulses_LaBr_df['timeF'].iloc[-1] - sync_pulses_LaBr_df['timeF'].iloc[0])*1e-9/60, "minutes")
    print('\033[95m' +'____________ Discrepancy Sync Pulses (LaBr3:Ce - POLARIS) ____________\n', excess_LaBr_sync_pulses)

    # remove excess sync pulses from the end of the dataframe because the POLARIS continued to send pulses after DAQ had stopped and LaBr3 was still recording
    if (excess_LaBr_sync_pulses > 0):
        sync_pulses_LaBr_df = sync_pulses_LaBr_df.iloc[:-(excess_LaBr_sync_pulses)]
    else :
        mod_sync_pulses_df = mod_sync_pulses_df.iloc[:-(abs(excess_LaBr_sync_pulses))]
    syncTimeDiff = sync_pulses_LaBr_df['timeF'].diff()
    sync_pulses_LaBr_df['syncTimeDiff'] = syncTimeDiff
    sync_pulses_LaBr_df['sync_flag'] = 122
    firstValidSyncLaBr = sync_pulses_LaBr_df['timeF'].iloc[0]
    
    print('\033[95m' +'____________ LaBr3:Ce sync pulses  ____________ \n', sync_pulses_LaBr_df)
    # sync_pulses_LaBr_df has columns : 'detectorID', 'globalTime', 'timeF', 'energyF', 'timeS', 'energyS', 'syncTimeDiff', 'sync_flag'

    sync_pulses_LaBr_df = sync_pulses_LaBr_df[['detectorID', 'globalTime', 'timeF', 'energyF', 'timeS', 'energyS', 'sync_flag', 'syncTimeDiff']]
    sync_pulses_LaBr_df = sync_pulses_LaBr_df.sort_values(by=['globalTime', 'timeF'], ascending=True)
    sync_pulses_LaBr_df.reset_index(drop=True, inplace=True)

    return sync_pulses_LaBr_df,firstValidSyncLaBr, mod_sync_pulses_df
    
def time_walk_correction(data_dir,mod_event_data_df,LaBr_data_df,sync_pulses_LaBr_df,firstValidSyncLaBr,mod_sync_pulses_df):   
    print('\033[96m' + '________________________________________') 
    print('\033[96m' + '_________ TIME WALK CORRECTION _________')
    print('\033[96m' + '________________________________________', '\n')  

    # LaBr_data_df has columns: 'detectorID', 'globalTime', 'timeF', 'energyF', 'timeS', 'energyS'

    # remove all detectorID == 5 (POLARIS) from LaBr3:Ce dataframe
    LaBr_data_df = LaBr_data_df.loc[LaBr_data_df['detectorID'] != 5]
    LaBr_data_df.reset_index(drop=True, inplace=True)

    LaBr_data_df.insert(0, 'sync_flag', 0)
    LaBr_data_df.insert(1, 'syncTimeDiff', 0) 
    LaBr_data_df = LaBr_data_df[['detectorID', 'globalTime', 'timeF', 'energyF', 'timeS', 'energyS', 'sync_flag', 'syncTimeDiff']]

    # sort by globalTime
    LaBr_data_df = LaBr_data_df.sort_values(by=['globalTime', 'timeF'], ascending=True)
    LaBr_data_df.reset_index(drop = True, inplace = True)
    LaBr_data_df['timeDiffLaBr'] = LaBr_data_df.loc[(LaBr_data_df['sync_flag'] == 0),'timeF'].diff()

    # INSERT the sync pulses in the LaBr3:Ce data (POLARIS is detector 5, which now also refers to sync pulses in general)
    sync_pulses_LaBr_df['timeDiffLaBr'] = 0
    LaBr_data_df = pd.concat([LaBr_data_df, sync_pulses_LaBr_df], ignore_index = True)
    LaBr_data_df = LaBr_data_df.drop(['energyF'], axis = 1)
    LaBr_data_df = LaBr_data_df.drop(['timeS'], axis = 1)
    LaBr_data_df['timeDiffLaBr'] = LaBr_data_df['timeDiffLaBr'].fillna(0)    

    # sort by ALL TIME
    LaBr_data_df.sort_values(by = ['globalTime', 'timeF'], inplace = True)
    LaBr_data_df['sync_flag'] = LaBr_data_df['sync_flag'].astype(int)
    LaBr_data_df = LaBr_data_df.loc[(LaBr_data_df['timeF'] >= sync_pulses_LaBr_df['timeF'].iloc[0]) & (LaBr_data_df['timeF'] <= sync_pulses_LaBr_df['timeF'].iloc[-1])]
    LaBr_data_df.reset_index(drop = True, inplace = True)

    sync_pulses_LaBr_df = sync_pulses_LaBr_df.drop(['detectorID', 'globalTime', 'energyF', 'timeS', 'energyS', 'syncTimeDiff', 'sync_flag'], axis = 1)
    sync_pulses_LaBr_df.reset_index(drop = True, inplace = True)    

    sync_pulses_LaBr_df = pd.concat([sync_pulses_LaBr_df['timeF'], mod_sync_pulses_df['sync_timestamp']], axis = 1)    
    sync_pulses_LaBr_df.columns = ['LaBr3_time', 'POLARIS_time']
    sync_pulses_LaBr_df['index'] = LaBr_data_df[LaBr_data_df['sync_flag'] == 122].index 
    sync_pulses_LaBr_df.set_index('index', inplace=True)

    print('\033[96m' +'____________ LaBr3:Ce sync pulses  ____________ \n', sync_pulses_LaBr_df)

    LaBr_data_df.loc[(LaBr_data_df['sync_flag'] == 0),'syncTimeDiff'] = np.nan
    LaBr_data_df['syncTimeDiff'] = LaBr_data_df['syncTimeDiff'].fillna(method = 'bfill')
    LaBr_data_df.loc[(LaBr_data_df['sync_flag'] == 0),'timeF'] =  (LaBr_data_df.loc[(LaBr_data_df['sync_flag'] == 0),'timeF']/LaBr_data_df.loc[(LaBr_data_df['sync_flag'] == 0),'syncTimeDiff'])*2000000060 # divide the timeF by the syncTimeDiff and multiply by 2000000060 to get the timeF corrected for time walk
    
    index_list = sync_pulses_LaBr_df.index.tolist()
    len_index_list = len(index_list)
    LaBr_data_df.loc[LaBr_data_df['sync_flag'] == 122, 'timeF'] = sync_pulses_LaBr_df['POLARIS_time'].values

    for i in range(len_index_list):
        if i < len_index_list - 1:
            LaBr_data_df.loc[index_list[i]:index_list[i+1]-1, 'timeF'] = LaBr_data_df.loc[index_list[i], 'timeF'] + LaBr_data_df.loc[index_list[i]:index_list[i+1]-1, 'timeDiffLaBr'].cumsum() 
        else:
            LaBr_data_df.loc[index_list[i]:, 'timeF'] = LaBr_data_df.loc[index_list[i], 'timeF'] + LaBr_data_df.loc[index_list[i]:, 'timeDiffLaBr'].cumsum()
    
    LaBr_data_df = LaBr_data_df.loc[(LaBr_data_df['sync_flag'] != 122) & (LaBr_data_df['detectorID'] != 5)]
    LaBr_data_df.drop(['sync_flag', 'syncTimeDiff', 'timeDiffLaBr', 'globalTime'], axis = 1, inplace = True)
    LaBr_data_df['x'] = np.nan
    LaBr_data_df['y'] = np.nan
    LaBr_data_df['z'] = np.nan
    LaBr_data_df.rename(columns = {'energyS': 'energy', 'timeF': 'time'}, inplace = True)
    LaBr_data_df.sort_values(by = ['time'], inplace = True)
    LaBr_data_df.reset_index(drop=True, inplace=True)

    first_time = LaBr_data_df.loc[LaBr_data_df['time'].notnull()].iloc[0]['time']
    last_time = LaBr_data_df.loc[LaBr_data_df['time'].notnull()].iloc[-1]['time']
    time_duration_labr = (last_time - first_time)*1e-9/60
    LaBr_data_df = LaBr_data_df[['detectorID', 'energy', 'time', 'x', 'y', 'z']]
    print('\033[96m' +'____________ LaBr3:Ce df time walk corrected ____________ \n', LaBr_data_df)
    print('\033[96m' +'____________ LaBr3:Ce time duration ____________ \n', time_duration_labr, 'minutes')
    
    #LaBr_data_df has columns: 'detectorID', 'energy', 'time', 'x', 'y', 'z'

    return LaBr_data_df

def merge_two_detector_dataframes(data_dir,mod_event_data_df,LaBr_data_df,run_num):
    print('\033[89m' + '________________________________________') 
    print('\033[89m' + '_____________ MERGE POLARIS AND LaBr3:Ce DATAFRAMES _____________')
    print('\033[89m' + '________________________________________', '\n')

    mod_event_data_df = mod_event_data_df.loc[mod_event_data_df['scatters'] != 122]
    mod_event_data_df.reset_index(drop = True, inplace = True)
    mod_event_data_df.drop(['scatters'], axis = 1, inplace = True)
    mod_event_data_df['detectorID'] = 5
    mod_event_data_df =  mod_event_data_df[['detectorID', 'energy', 'time', 'x', 'y', 'z']]
    mod_event_data_df.sort_values(by = ['time'], ascending=True, inplace = True) 
    mod_event_data_df.reset_index(drop = True, inplace = True)

    df_final = pd.concat([LaBr_data_df, mod_event_data_df], ignore_index = True, axis=0) # the dataframes are concatenated vertically
    df_final.sort_values(by = ['time'], ascending=True, inplace = True)
    df_final.reset_index(drop = True, inplace = True)
    
    print('\033[89m' +'____________ df_final  ____________ \n', df_final)
    print('\033[89m' +'____________ df_final time duration ____________ \n', (df_final['time'].iloc[-1] - df_final['time'].iloc[0])*1e-9/60, 'minutes')

    # save the final dataframe to a parquet file
    df_final.to_parquet(data_dir + 'run' + str(run_num) + '_mergedLaBrPOLARISdata.parquet')

    
    # ______________________________

    # profile = ProfileReport(df_final, title="Profiling Report")
    # import subprocess
    # profile.to_file("profiling_report.html")
    # html_file = 'profiling_report.html'
    # subprocess.call(['cp', html_file, data_dir])

    return df_final


def time_diff_events(data_dir, LaBr_data_df,mod_event_data_df):
    polaris_time_dff = mod_event_data_df['time'].diff()
    polaris_time_dff = pd.DataFrame(polaris_time_dff/1000000) # microseconds
    LaBr0 = LaBr_data_df.loc[LaBr_data_df['detectorID'] == 0]
    slow_time_dff = LaBr0['timeS'].diff()
    fast_time_dff = LaBr0['timeF'].diff()
    slow_time_dff = pd.DataFrame(slow_time_dff/1000000) # microseconds
    fast_time_dff = pd.DataFrame(fast_time_dff/1000000) # microseconds
    polaris_time_dff = polaris_time_dff.iloc[1:]
    slow_time_dff = slow_time_dff.iloc[1:]
    fast_time_dff = fast_time_dff.iloc[1:]
    print('\033[96m' +'____________  Time Diff Between Consecutive Events  ____________ \n', 'POLARIS: \n', polaris_time_dff, '\n LaBr3:Ce Slow Signal: \n', slow_time_dff, '\n LaBr3:Ce Fast Signal: \n', fast_time_dff)

    fig = plt.figure()
    plt.plot(polaris_time_dff, label = 'POLARIS')
    plt.xlabel('Time difference ($\mu$s)', fontsize = 25)
    plt.ylabel('Counts', fontsize = 25)
    plt.title('Time Difference Between Consecutive Events')
    plt.legend(loc = 'upper right', fontsize = 30)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 20)
    plt.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    plt.show()

    fig2 = plt.figure()
    plt.plot(slow_time_dff, label = 'LaBr$_{3}$ Slow Signal')
    plt.xlabel('Time difference ($\mu$s)', fontsize = 25)
    plt.ylabel('Counts', fontsize = 25)
    plt.title('Time Difference Between Consecutive Events')
    plt.legend(loc = 'upper right', fontsize = 30)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 20)
    plt.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    plt.show()

    fig3 = plt.figure()
    plt.plot(fast_time_dff, label = 'LaBr$_{3}$ Fast Signal')
    plt.xlabel('Time difference ($\mu$s)', fontsize = 25)
    plt.ylabel('Counts', fontsize = 25)
    plt.title('Time Difference Between Consecutive Events')
    plt.legend(loc = 'upper right', fontsize = 30)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 20)
    plt.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    plt.show()

def doubles(df_final):

    df_final = df_final.loc[(df_final['detectorID'] == 5) | (df_final['detectorID'] == 0)]
    df_final.reset_index(drop=True, inplace=True)

    # Initialize variables to track consecutive pairs
    pairs = []
    current_pair = []

    for index, row in df_final.iterrows():
        if row['detectorID'] == 5:
            current_pair = [row]
        elif row['detectorID'] == 0 and current_pair:
            current_pair.append(row)
            pairs.append(current_pair)
            current_pair = []

    # Create histograms for each pair and plot
    for pair in pairs:
        energy_values = pair[0]['energy'], pair[1]['energy']
        plt.hist(energy_values, bins=10, alpha=0.5, label=f'Detector {pair[0]["detectorID"]} to {pair[1]["detectorID"]}')

    plt.xlabel('Energy (keV) of Detector 0 and 5')
    plt.legend()
    plt.title('Histogram of Energy for Detector Pairs')
    plt.show()

    # Create a histogram of the energy detector 0 vs energy detector 5
    plt.hist2d(df_final.loc[df_final['detectorID'] == 5]['energy'], df_final.loc[df_final['detectorID'] == 0]['energy'], bins=1000, range=[[0, 1000], [0, 1000]])
    plt.xlabel('POLARIS Energy (keV)')
    plt.ylabel('LaBr3:Ce Energy (keV)')
    cbar = plt.colorbar(label='Frequency')
    cbar.cmap.set_under('white')
    plt.title('LaBr3:Ce Energy vs POLARIS Energy (2D Histogram)')
    plt.show()




    
# ______________________________________________________________________________________________________________________
def main():
    start = time.time()
    data_dir = sys.argv[1]
    run_num, mod_event_data_df,mod_event_data = read_polaris_data(data_dir)
    LaBr_data_df = read_LaBr_data(data_dir) 

    mod_event_data_df = apply_coordinate_transformations(run_num, mod_event_data_df)
    mod_sync_pulses_df,firstValidSyncPOLARIS,mod_event_data_df = find_POLARIS_sync_pulses(mod_event_data,mod_event_data_df, data_dir)
    sync_pulses_LaBr_df,firstValidSyncLaBr, mod_sync_pulses_df = find_LaBr3_sync_pulses(firstValidSyncPOLARIS,LaBr_data_df,mod_sync_pulses_df)
    time_diff_events(data_dir, LaBr_data_df,mod_event_data_df)
    LaBr_data_df = time_walk_correction(data_dir, mod_event_data_df, LaBr_data_df, sync_pulses_LaBr_df,firstValidSyncLaBr,mod_sync_pulses_df)
    df_final = merge_two_detector_dataframes(data_dir,mod_event_data_df,LaBr_data_df,run_num)
    #doubles(df_final)

    end = time.time()

    print('\033[92m' + '________________________________________')
    print('\033[92m' + 'Run time: ', (end - start)/60, 'minutes')

if __name__ == '__main__':
    main()

