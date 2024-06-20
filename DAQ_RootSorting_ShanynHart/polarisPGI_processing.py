#!/usr/bin/python
"""
Description:
  Python script that reads in raw data from Polaris (typically named: AllEventsCombined.txt)
    and can: perform a coordinate transformation on the data
             does Big Energy sorting on the data
             do Compton line filtering on the data
             remove non-physical events (double scatters only)
             randomly shuffle events (doubles & triples only)
             look for coincidence events (multiple-stage) - not yet functioning
             creates csv data files
             output lots of interesting plots

  To Run: python polarisPGI_processing.py &> output-log.txt &
    From main program, call function: filter_data(data_dir, data_filename, transforms)
              data_dir            - Required : Path to combined detector data .txt file (Str)
              data_filename       - Required : Name of combined detector data .txt file to read data from (Str)
              transforms          - Required : Array of transformation arrays, used for coordinate transformations (Float[][])
                                               format: (rot[0], rot[1], rot[2], pos[0], pos[1], pos[2]) -> units are degrees and mm

    Data (and plots) will be saved to the given data directory and plots are
      saved to the following directories (plots_raw, plots_transformed, plots_final_doubles)

    Additional settings are handled through Global Variable under SETTINGS / PARAMETERS (just before main)
      - PRODUCE_RAW_PLOTS                      = True (create raw plots) / False (do not create plots)
      - PRODUCE_TRANSFORMED_PLOTS              = True (create transformed plots) / False (do not create plots)
      - PRODUCE_FINAL_DOUBLE_PLOTS             = True (create final double energy plots) / False (do not create plots)
      - APPLY_BIG_ENERGY_SORTING               = apply big energy sorting (automatically arrange events so that bigger energy deposited comes first, default is False)
      - REMOVE_UNPHYSICAL_EVENTS               = remove unphysical events, i.e. Compton scatter angle == nan (checks both ordering of two scatters and flips if original ordering in unphysical, default is False)
      - APPLY_COMPTON_LINE_FILTERING           = apply Compton line filtering (checks if E1 and theta1 are consistent with Compton formula, default is False)
      - COMPTON_LINE_RANGE                     = sets accepted energy range (min / max values (percentage) for Compton line filtering)
      - COMPTON_LINE_ENERGIES                  = sets expected gamma energies to use for filtering (can take multiple arguments, see list of options below)
      - SHUFFLE_SCATTER_EVENTS                 = shuffle the 2x & 3x scatter events (to overcome timing errors like runs 6, 7, 8 for 180315 data, default is False)
      - APPLY_COINCIDENCE_GROUPING_MULTISTAGE  = groups events from different modules based on energy and time stamp (NOT FUNCTIONING, default is False)

    Possible gamma energies for Compton Line Filtering
      - Oxygen Prompts: 2.742, 5.240, 6.129, 6.916, 7.116
      - Carbon: 4.444
      - Nitrogen: 1.635, 2.313, 5.269, 5.298
      - Boron: 0.718
      - Cobalt-60: 1.173, 1.332
      - Cesium-137: 0.6617
      - Sodium-22 : 0.5110, 1.274

NOTES:
- written using python v2.7.15 / updated to also work in python 3.6.2
- processPolarisData_utils.pyx required to run this code
- some aspects of the code only function for double-scatters

Code borrowed heavily from Dennis Mackin <dsmackin@mdanderson.org> & Nicholas Hyslop <HYSNIC007@myuct.ac.za>
"""
__author__ = "Steve Peterson <steve.peterson@uct.ac.za>"
__date__ = "July 29, 2019"
__version__ = "$Revision: 1.0.0$"

#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------

import os
import sys
import timeit

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import cython utils library from polaris_data_processor
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
# sys.path.append('/Users/nicholas/Google Drive/Work/code/polaris/polaris_data_processor/')       # Nicholas latop
# sys.path.append('/home/nicholas/polarisPEPT/polaris/polaris_data_processor/')                 # Steve server
# sys.path.append('/Volumes/BATMAN/PGI/code/polaris/polaris_data_processor/')                 # Steve laptop
import processPolarisData_utils as utils

#import pyximport; pyximport.install()
#import polarisPGI_processing_utils as utils


#------------------------------------------------------------------
# FUNCTION/CLASS DEFINITIONS
#------------------------------------------------------------------

# Read in polaris data from a text file and convert into pandas array, including transformations
def read_polaris_data(data_dir, data_filename, transformations):
    '''
    Takes in an AllEventsCombined.txt file and reads the information into a
    python pandas array. It then applies a coordinate transformation according
    to transformations parameter, and returns the DataFrame.

    @params:
        data_dir            - Required : Path to combined detector data .txt file (Str)
        data_filename       - Required : Name of combined detector data .txt file to read data from (Str)
        transforms          - Required : Array of transformation arrays, used for coordinate transformations (Float[][])
                                         format: (rot[0], rot[1], rot[2], pos[0], pos[1], pos[2]) -> units are degrees and mm

    @returns:
        allevents_data      - pandas array of transformed AllEventsCombined.txt data
    '''

    print('  Reading in data from', data_dir + data_filename + '.', '\n  This may take several minutes...')
    allevents_data = pd.read_csv(data_dir + data_filename, sep = '	', header = None)                        # read in text file as a csv with tab separation
    print('  Done.')
    allevents_data.columns = ['scatters', 'detector', 'energy', 'x', 'y', 'z', 'time']       # separate into columns and name them

    print()
    print('{:30} {:10} {:10} {:10}'.format('', 'All', 'D0', 'D1'))
    evts_all = len(allevents_data.index);  evts_d0 = len(allevents_data[allevents_data.detector == 0]);  evts_d1 = len(allevents_data[allevents_data.detector == 1]);
    print('{:30} {:<10d} {:<10d} {:<10d}'.format('Total number of events:', evts_all, evts_d0, evts_d1))
    time_total = (allevents_data.iloc[-1]['time'] - allevents_data.iloc[0]['time']) * 1e-8
    print('{:30} {:<10.1f}'.format('Total time (seconds):', time_total))
    evt_rate_all = evts_all / time_total;  evt_rate_D0 = evts_d0 / time_total;  evt_rate_D1 = evts_d1 / time_total
    print('{:30} {:<10.1f} {:<10.1f} {:<10.1f}'.format('Event rate (events/sec):', evt_rate_all, evt_rate_D0, evt_rate_D1))

    # Plotting basic detector plots using raw data
    if PRODUCE_RAW_PLOTS:
        plot_raw_data(data_dir, data_filename, allevents_data)

    print('\nApplying transformations to convert data into isocentric coordinate system...')
    print('-----------------------------------------------------------------------------')
    transformation_matrices = {}

    # Creating transformation matrices
    for i, arr in enumerate(transformations):
        transformation_matrices[i] = utils.get_transformation_matrix_array(arr)

    #  applying coordinate transformation to data
    allevents_data = utils.apply_transformation(allevents_data, transformation_matrices)

    #  convert detector, scatters, time columns from float (default) into ints
    allevents_data['scatters'] = allevents_data['scatters'].astype(int)
    allevents_data['detector'] = allevents_data['detector'].astype(int)
    allevents_data['time'] = allevents_data['time'].astype(int)

    print('\n  Done.')

    # print(allevents_data.iloc[:3])

    return allevents_data


# Plotting basic detector plots using raw data
def add_tags_to_filename(data_filename, ending):

    tagged_filename = data_filename.replace(".txt", "").replace(".dat", "").replace(".csv", "") + "-scatters"
    if APPLY_BIG_ENERGY_SORTING:                tagged_filename = tagged_filename + "-BigE"
    if REMOVE_UNPHYSICAL_EVENTS:                tagged_filename = tagged_filename + "-UnPhy"
    if APPLY_COMPTON_LINE_FILTERING:            tagged_filename = tagged_filename + "-ComL"
    if SHUFFLE_SCATTER_EVENTS:                  tagged_filename = tagged_filename + "-Shuf"
    if APPLY_COINCIDENCE_GROUPING_MULTISTAGE:   tagged_filename = tagged_filename + "-MScoin"
    return tagged_filename + ending


# Plotting basic detector plots using raw data
def plot_raw_data(data_dir, data_filename, cc_dataframe):

    #  checking if plots directory exists, if not, create
    plots_dir = data_dir + add_tags_to_filename(data_filename, '-plots_raw/')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print ('Creating directory: {}'.format(plots_dir))


    ### PLOTTING BASIC DETECTOR DATA ###
    print ('\n--- Making Basic Detector Plots (using raw data) ---')

    #  plotting raw data (1D profiles in x, y, z, energy + 2D profiles in xy, yz, xz)
    utils.make_basic_detector_plots(cc_dataframe, 'raw_all', plots_dir)
    #  alternate plots for raw data of 2D xy profiles
    d0_xPos = cc_dataframe[cc_dataframe.detector == 0].x.to_numpy()
    d0_yPos = cc_dataframe[cc_dataframe.detector == 0].y.to_numpy()
    d1_xPos = cc_dataframe[cc_dataframe.detector == 1].x.to_numpy()
    d1_yPos = cc_dataframe[cc_dataframe.detector == 1].y.to_numpy()
    plt.cla(); plt.clf()
    plt.plot(d0_yPos, d0_xPos, color='green', marker='o', linestyle='None', markersize=0.5)
    plt.title('XY Position of all scatters / D0 / Raw Data, total: %s'%(len(d0_yPos)))
    plt.xlabel('Y-Position (mm)'); plt.ylabel('X-Position (mm)')
    utils.save_to_file("raw_all_D0_y_x_alternate", plots_dir)
    plt.cla(); plt.clf()
    plt.plot(d1_yPos, d1_xPos, color='green', marker='o', linestyle='None', markersize=0.5)
    plt.title('XY Position of all scatters / D1 / Raw Data, total: %s'%(len(d1_yPos)))
    plt.xlabel('Y-Position (mm)'); plt.ylabel('X-Position (mm)')
    utils.save_to_file("raw_all_D1_y_x_alternate", plots_dir)
    #  alternate plots for energy deposition of raw data
    alternate_max_energy = (np.amax(COMPTON_LINE_ENERGIES) * 1000) + 100   # convert max Compton Line energy to keV and add 100 keV
    print ('  Creating 1D Plots - using alternate maximum energy of {} keV . . .'.format(alternate_max_energy))
    plt.cla(); plt.clf()
    utils.create_plot_range(cc_dataframe.energy, 'Energy Deposited (keV)', 300, 0, alternate_max_energy)
    plt.title('Total energy deposted by all scatters, total: %s'%(len(cc_dataframe.energy)))
    utils.save_to_file("raw_all_Energy_Counts_alternate", plots_dir)
    plt.cla(); plt.clf()
    d0_eng = cc_dataframe[cc_dataframe.detector == 0].energy.to_numpy()
    utils.create_plot_range(d0_eng, 'Energy Deposited (keV)', 300, 0, alternate_max_energy)
    plt.title('Total energy deposted by all scatters / D0, total: %s'%(len(d0_eng)))
    utils.save_to_file("raw_all_D0_Energy_Counts_alternate", plots_dir)
    d1_eng = cc_dataframe[cc_dataframe.detector == 1].energy.to_numpy()
    plt.cla(); plt.clf()
    utils.create_plot_range(d1_eng, 'Energy Deposited (keV)', 300, 0, alternate_max_energy)
    plt.title('Total energy deposted by all scatters / D1, total: %s'%(len(d1_eng)))
    utils.save_to_file("raw_all_D1_Energy_Counts_alternate", plots_dir)


# Plotting basic detector plots using transformed data
def plot_transformed_data(data_dir, data_filename, cc_dataframe):

    #  checking if plots directory exists, if not, create
    plots_dir = data_dir + add_tags_to_filename(data_filename, '-plots_transformed/')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print ('  Creating directory: {}'.format(plots_dir))

    ### PLOTTING TRANSFORMED DETECTOR DATA ###
    print ('\n--- Making Basic Detector Plots (using transformed data) ---')

    #  plotting transformed data (1D profiles in x, y, z, energy + 2D profiles in xy, yz, xz)
    utils.make_basic_detector_plots(cc_dataframe, 'transformed_all', plots_dir)
    #  plotting based on number of pixel interactions (1, 2, 3, 4+)
    utils.make_basic_detector_plots(cc_dataframe[cc_dataframe.scatters == 1], 'transformed_1px', plots_dir)
    utils.make_basic_detector_plots(cc_dataframe[cc_dataframe.scatters == 2], 'transformed_2px', plots_dir)
    utils.make_basic_detector_plots(cc_dataframe[cc_dataframe.scatters == 3], 'transformed_3px', plots_dir)
    utils.make_basic_detector_plots(cc_dataframe[cc_dataframe.scatters > 3], 'transformed_4px+', plots_dir)
    #  plotting 1D profiles (x, y, z) across all detectors
    utils.plot_1D(cc_dataframe.x.to_numpy(), 200, "x", "Counts", "transformed_all", plots_dir)
    utils.plot_1D(cc_dataframe.y.to_numpy(), 200, "y", "Counts", "transformed_all", plots_dir)
    utils.plot_1D(cc_dataframe.z.to_numpy(), 200, "z", "Counts", "transformed_all", plots_dir)


# Plotting detector plots using final data
def plot_final_doubles_data(data_dir, data_filename, scatters_2x_out, scatters_2x_d0_out, scatters_2x_d1_out):

    #  checking if plots directory exists, if not, create
    plots_dir = data_dir + add_tags_to_filename(data_filename, '-plots_final_doubles/')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print ('  Creating directory: {}'.format(plots_dir))

    ### PLOTTING FINAL DETECTOR DATA ###
    print ('\n--- Making Basic Detector Plots (using final doubles data) ---')

    #  producing required data
    MeCsq = 0.5109989461  # electron mass in energy units (MeV)
    #   - scatters_2x_out
    ds_1st_eng = scatters_2x_out[:, 0];  ds_2nd_eng = scatters_2x_out[:, 4];  ds_tot_eng = ds_1st_eng + ds_2nd_eng
    ds_all_eng = np.concatenate((ds_1st_eng, ds_2nd_eng), axis = 0)
    ds_theta1 = np.arccos( 1 + MeCsq * ( 1.0/(ds_tot_eng) - 1.0/(ds_tot_eng - ds_1st_eng) ) )
    #   - scatters_2x_d0_out
    ds_d0_1st_eng = scatters_2x_d0_out[:, 0];  ds_d0_2nd_eng = scatters_2x_d0_out[:, 4];  ds_d0_tot_eng = ds_d0_1st_eng + ds_d0_2nd_eng
    ds_d0_all_eng = np.concatenate((ds_d0_1st_eng, ds_d0_2nd_eng), axis = 0)
    ds_d0_theta1 = np.arccos( 1 + MeCsq * ( 1.0/(ds_d0_tot_eng) - 1.0/(ds_d0_tot_eng - ds_d0_1st_eng) ) )
    #   - scatters_2x_d1_out
    ds_d1_1st_eng = scatters_2x_d1_out[:, 0];  ds_d1_2nd_eng = scatters_2x_d1_out[:, 4];  ds_d1_tot_eng = ds_d1_1st_eng + ds_d1_2nd_eng
    ds_d1_all_eng = np.concatenate((ds_1st_eng, ds_d1_2nd_eng), axis = 0)
    ds_d1_theta1 = np.arccos( 1 + MeCsq * ( 1.0/(ds_d1_tot_eng) - 1.0/(ds_d1_tot_eng - ds_d1_1st_eng) ) )

    alternate_max_energy = np.amax(COMPTON_LINE_ENERGIES) + 0.1   # use max Compton Line energy and add 100 keV
    print ('  Creating 1D Plots - using alternate maximum energy of {} MeV . . .'.format(alternate_max_energy))

    #  plots for energy deposition of final data
    #   - scatters_2x_out
    plt.cla(); plt.clf()
    utils.create_plot_range(ds_tot_eng, 'Energy Deposited (MeV)', 300, 0, alternate_max_energy)
    plt.title('Energy deposited by all double scatters (summed total), total: %s'%(len(ds_tot_eng)))
    utils.save_to_file("final_doubles_all_Energy_Sum_Counts_alternate", plots_dir)
    plt.cla(); plt.clf()
    utils.create_plot_range(ds_all_eng, 'Energy Deposited (MeV)', 300, 0, alternate_max_energy)
    plt.title('Energy deposited by all double scatters (individual scatters), total: %s'%(len(ds_all_eng)))
    utils.save_to_file("final_doubles_all_Energy_Each_Counts_alternate", plots_dir)
    #   - scatters_2x_d0_out
    plt.cla(); plt.clf()
    utils.create_plot_range(ds_d0_tot_eng, 'Energy Deposited (MeV)', 300, 0, alternate_max_energy)
    plt.title('Energy deposited by D0 double scatters (summed total), total: %s'%(len(ds_d0_tot_eng)))
    utils.save_to_file("final_doubles_D0_Energy_Sum_Counts_alternate", plots_dir)
    plt.cla(); plt.clf()
    utils.create_plot_range(ds_d0_all_eng, 'Energy Deposited (MeV)', 300, 0, alternate_max_energy)
    plt.title('Energy deposited by D0 double scatters (individual scatters), total: %s'%(len(ds_d0_all_eng)))
    utils.save_to_file("final_doubles_D0_Energy_Each_Counts_alternate", plots_dir)
    #   - scatters_2x_d1_out
    plt.cla(); plt.clf()
    utils.create_plot_range(ds_d1_tot_eng, 'Energy Deposited (MeV)', 300, 0, alternate_max_energy)
    plt.title('Energy deposited by D1 double scatters (summed total), total: %s'%(len(ds_d1_tot_eng)))
    utils.save_to_file("final_doubles_D1_Energy_Sum_Counts_alternate", plots_dir)
    plt.cla(); plt.clf()
    utils.create_plot_range(ds_d1_all_eng, 'Energy Deposited (MeV)', 300, 0, alternate_max_energy)
    plt.title('Energy deposited by D1 double scatters (individual scatters), total: %s'%(len(ds_d1_all_eng)))
    utils.save_to_file("final_doubles_D1_Energy_Each_Counts_alternate", plots_dir)

    #  plotting E1 vs Theta1
    print ('  Creating energy 1 vs theta 1 plots . . .')
    alternate_max_energy_cl_plots = np.amax(COMPTON_LINE_ENERGIES)   # use maxCompton Line energy
    print ('    - using alternate maximum energy of {} MeV . . .'.format(alternate_max_energy_cl_plots))
    #   - scatters_2x_out
    plt.cla(); plt.clf()
    plt.plot(ds_theta1, ds_1st_eng, color='blue', marker='o', linestyle='None', markersize=0.5)
    plt.title('First Energy Deposited vs Calculated First Scatter Angle / All, total: %s'%(len(ds_theta1)))
    plt.xlabel('Calculated First Scatter Angle (rad)');  plt.ylabel('Energy Deposited in First Scatter (MeV)')
    plt.axis([0, 3, 0, alternate_max_energy_cl_plots])
    utils.save_to_file("final_doubles_all_Energy1_Theta1", plots_dir)
    #   - scatters_2x_d0_out
    plt.cla(); plt.clf()
    plt.plot(ds_d0_theta1, ds_d0_1st_eng, color='blue', marker='o', linestyle='None', markersize=0.5)
    plt.title('First Energy Deposited vs Calculated First Scatter Angle / D0, total: %s'%(len(ds_d0_theta1)))
    plt.xlabel('Calculated First Scatter Angle (rad)');  plt.ylabel('Energy Deposited in First Scatter (MeV)')
    plt.axis([0, 3, 0, alternate_max_energy_cl_plots])
    utils.save_to_file("final_doubles_D0_Energy1_Theta1", plots_dir)
    #   - scatters_2x_d1_out
    plt.cla(); plt.clf()
    plt.plot(ds_d1_theta1, ds_d1_1st_eng, color='blue', marker='o', linestyle='None', markersize=0.5)
    plt.title('First Energy Deposited vs Calculated First Scatter Angle / D1, total: %s'%(len(ds_d1_theta1)))
    plt.xlabel('Calculated First Scatter Angle (rad)');  plt.ylabel('Energy Deposited in First Scatter (MeV)')
    plt.axis([0, 3, 0, alternate_max_energy_cl_plots])
    utils.save_to_file("final_doubles_D1_Energy1_Theta1", plots_dir)


# Process the all_events data and extract coincidences
def filter_data(data_dir, data_filename, transforms):
    '''
    Takes in an AllEventsCombined.txt file, reads the information into a pandas
    array and then filters data using the flagged parameters

    @params:
        data_dir            - Required : Path to combined detector data .txt file (Str)
        data_filename       - Required : Name of combined detector data .txt file to read data from (Str)
        transforms          - Required : Array of transformation arrays, used for coordinate transformations (Float[][])
                                         format: (rot[0], rot[1], rot[2], pos[0], pos[1], pos[2]) -> units are degrees and mm

    @returns:
        scatters_2x_out     - numpy array of the processed double scatter data recorded by both polaris detectors
        scatters_2x_d0_out  - numpy array of the processed double scatter data recorded by polaris detector 0
        scatters_2x_d1_out  - numpy array of the processed double scatter data recorded by polaris detector 1
        scatters_3x_out     - numpy array of the processed triple scatter data recorded by both polaris detectors
    '''

    #  measuring time required to run code - start time
    start = timeit.default_timer()

    allevents_data = read_polaris_data(data_dir, data_filename, transforms)       # read in the raw polaris data
    tot_raw_events = len(allevents_data.index)          # get the total number of events before coincidence processing

    #  output time details
    #   - convert time into seconds (t given in units of 10 ns clock cycles), so tS = t * 10E-9
    time_in_sec = allevents_data['time'] * (10E-9)
    time_start = np.amin(time_in_sec)
    time_end = np.amax(time_in_sec)
    #   - print to screen
    print ('   - Start time (s): {} | End time (s): {} | Run duration (s): {}'.format(time_start, time_end, (time_end - time_start)))
    print ('    - run activity (events/s): {}'.format(len(allevents_data) / (time_end - time_start)))

    # Counters for run summary
    num_1px, num_2px, num_3px, num_4px = [float('nan')] * 4                # count events based on pixel number
    num_1px_d0, num_1px_d1 = [float('nan')] * 2                            # count single pixel events based on detector
    num_ss, num_ss_et, num_ss_fi = [float('nan')] * 3                      # count events for single scatters
    num_ds, num_ds_d0, num_ds_d1, num_ts = [float('nan')] * 4              # count total number of double/triple scatters
    num_ds_et, num_ds_d0_et, num_ds_d1_et, num_ts_et = [float('nan')] * 4  # count number of energy threshold events
    num_ds_pe, num_ds_d0_pe, num_ds_d1_pe = [float('nan')] * 3             # count number of physical events
    num_ds_cl, num_ds_d0_cl, num_ds_d1_cl = [float('nan')] * 3             # count number of Compton line filtered events
    num_ds_ms, num_ds_d0_ms, num_ds_d1_ms = [float('nan')] * 3             # count number of multistage coincidence events
    num_ds_fi, num_ds_d0_fi, num_ds_d1_fi, num_ts_fi = [float('nan')] * 4  # final number of double/triple scatters

    # count number of pixel events
    num_1px = len(allevents_data[allevents_data.scatters == 1]);  num_2px = len(allevents_data[allevents_data.scatters == 2])
    num_1px_d0 = len(allevents_data[(allevents_data.scatters == 1) & (allevents_data.detector == 0)]);  num_1px_d1 = len(allevents_data[(allevents_data.scatters == 1) & (allevents_data.detector == 1)]);
    num_3px = len(allevents_data[allevents_data.scatters == 3]);  num_4px = len(allevents_data[allevents_data.scatters  > 3])

    print('\nFiltering data')
    print('-----------------------')

    # Plotting basic detector plots using transformed data
    if PRODUCE_TRANSFORMED_PLOTS:
        plot_transformed_data(data_dir, data_filename, allevents_data)

    print ('\n--- Grouping Compton Scatter Data ---')

    if APPLY_BIG_ENERGY_SORTING:
        print ('  Sorting scatters to place largest energy deposition first . . .')

    #  returns interaction data (E, X, Y, Z) for the specified number of scatters and detector number
    #   format -> def get_interaction_data(df, energy_sort = True, num_scatters = 2, det_num = None)
    scatters_1x = utils.get_interaction_data(allevents_data, APPLY_BIG_ENERGY_SORTING, 1)
    scatters_2x = utils.get_interaction_data(allevents_data, APPLY_BIG_ENERGY_SORTING, 2)
    scatters_2x_d0 = utils.get_interaction_data(allevents_data, APPLY_BIG_ENERGY_SORTING, 2, 0)
    scatters_2x_d1 = utils.get_interaction_data(allevents_data, APPLY_BIG_ENERGY_SORTING, 2, 1)
    scatters_3x = utils.get_interaction_data(allevents_data, APPLY_BIG_ENERGY_SORTING, 3)

    ### FORMATTING DATA FOR OUTPUT ###

    print ('\n--- Formatting Compton Scatter Data for Output ---')
    print ('  Re-arranging data array for CSV output / re-scaling energy from keV to MeV . . .')

    #  re-arranges interaction data into CSV format (entire event on one line)
    #   format -> def format_data_for_output(scatters, num_scatters, det_num = None):
    scatters_1x_out = scatters_1x  # no rearrangment required
    scatters_2x_out = utils.format_data_for_output(scatters_2x, 2)
    scatters_2x_d0_out = utils.format_data_for_output(scatters_2x_d0, 2, 0)
    scatters_2x_d1_out = utils.format_data_for_output(scatters_2x_d1, 2, 1)
    scatters_3x_out = utils.format_data_for_output(scatters_3x, 3)

    # count initial number of double/triple scatters
    num_ss = len(scatters_1x)
    num_ds = len(scatters_2x_out);  num_ds_d0 = len(scatters_2x_d0_out)
    num_ds_d1 = len(scatters_2x_d1_out);  num_ts = len(scatters_3x_out)

    ### FILTERING DETECTOR DATA ###

    print ('\n--- Filtering Detector Data ---')

    #  removing unphysical events (Compton scatter angle == nan)
    #   - checks both ordering of two scatters and flips if original ordering in unphysical
    #   - modified code from Matt Leigh's Filter.py & CPUFunctions.cu
    if REMOVE_UNPHYSICAL_EVENTS:
        print ('  Removing unphysical events (double scatters) . . .')

        print ('   - Filtering {} | total events checked: {}'.format('scatters_2x_out', len(scatters_2x_out)))
        scatters_2x_out = utils.filtering_unphysical_double_scatters(scatters_2x_out)
        print ('   - Filtering {} | total events checked: {}'.format('scatters_2x_d0_out', len(scatters_2x_d0_out)))
        scatters_2x_d0_out = utils.filtering_unphysical_double_scatters(scatters_2x_d0_out)
        print ('   - Filtering {} | total events checked: {}'.format('scatters_2x_d1_out', len(scatters_2x_d1_out)))
        scatters_2x_d1_out = utils.filtering_unphysical_double_scatters(scatters_2x_d1_out)

        # count number of physical double scatters
        num_ds_pe = len(scatters_2x_out);  num_ds_d0_pe = len(scatters_2x_d0_out);  num_ds_d1_pe = len(scatters_2x_d1_out)


    #  applying Compton line filtering (checks if E1 and theta1 are consistent with Compton formula)
    #   - use COMPTON_LINE_RANGE (accepted energy range) & COMPTON_LINE_ENERGIES (expected gamma energies) global variables
    #   - modified code from Matt Leigh's Filter.py
    if APPLY_COMPTON_LINE_FILTERING:
        print ('\n  Applying Compton line filtering / Energies (MeV): {} / Range of values: {}'.format(COMPTON_LINE_ENERGIES, COMPTON_LINE_RANGE))

        print ('   - Filtering {} | total events checked: {}'.format('scatters_2x_out', len(scatters_2x_out)))
        scatters_2x_out = utils.compton_line_filtering(scatters_2x_out, COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES)
        print ('   - Filtering {} | total events checked: {}'.format('scatters_2x_d0_out', len(scatters_2x_d0_out)))
        scatters_2x_d0_out = utils.compton_line_filtering(scatters_2x_d0_out, COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES)
        print ('   - Filtering {} | total events checked: {}'.format('scatters_2x_d1_out', len(scatters_2x_d1_out)))
        scatters_2x_d1_out = utils.compton_line_filtering(scatters_2x_d1_out, COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES)

        # count number of Compton line filtered double scatters
        num_ds_cl = len(scatters_2x_out);  num_ds_d0_cl = len(scatters_2x_d0_out);  num_ds_d1_cl = len(scatters_2x_d1_out)


    #  shuffle the 2x & 3x scatter events
    #   - to overcome timing errors (like runs 6, 7, 8 for 180315 data)
    if SHUFFLE_SCATTER_EVENTS:

        #  - np.random.shuffle: Modify a sequence in-place by shuffling its contents.
        np.random.shuffle(scatters_2x_out)
        np.random.shuffle(scatters_3x_out)


    ### GROUPING COINCIDENCE EVENTS ###

    print ('\n--- Grouping Coincidence Events ---')

    #  groups events from different modules based on energy and time stamp
    #   - modified code from Paul Maggi's coincCheckMod function
    if APPLY_COINCIDENCE_GROUPING_MULTISTAGE:
        print ('  Grouping multi-stage coincidence events . . .')

        print ('    ######  CURRENTLY NOT FUNCTIONING  ######')

        # count number of multi-stage coincidence double scatters
        # num_ds_ms = len(scatters_2x_out);  num_ds_d0_ms = len(scatters_2x_d0_out);  num_ds_d1_ms = len(scatters_2x_d1_out)

    # Plotting detector plots using final doubles data
    if PRODUCE_FINAL_DOUBLE_PLOTS:
        plot_final_doubles_data(data_dir, data_filename, scatters_2x_out, scatters_2x_d0_out, scatters_2x_d1_out)

    ### SAVING DATA TO CSV ###

    print ('\n--- Saving Polaris Data to CSV ---')

    #  creating output file names
    #   - add tags based on filter parameters
    output_file = data_dir + add_tags_to_filename(data_filename, '.csv')
    output_file_2x = output_file.replace(".csv", "_2x.csv")
    output_file_2x_d0 = output_file.replace(".csv", "_2x_d0.csv")
    output_file_2x_d1 = output_file.replace(".csv", "_2x_d1.csv")
    output_file_3x = output_file.replace(".csv", "_3x.csv")
    output_file_2x_3x = output_file.replace(".csv", "_2x+3x.csv")

    #  outputting data to CSV
    #   format: eng1, x1, y1, z1, eng2, x2, y2, z2, (eng3, x3, y3, z3) <- if triple scatter
    #   units: time in us and pos in mm
    print ('  Saving {} events from data array ({}) to file (CSV format): {}'.format(len(scatters_2x_out), 'scatters_2x_out', output_file_2x))
    np.savetxt(output_file_2x, scatters_2x_out, delimiter=',', fmt='%.5f')
    print ('  Saving {} events from data array ({}) to file (CSV format): {}'.format(len(scatters_2x_d0_out), 'scatters_2x_d0_out', output_file_2x_d0))
    np.savetxt(output_file_2x_d0, scatters_2x_d0_out, delimiter=',', fmt='%.5f')
    print ('  Saving {} events from data array ({}) to file (CSV format): {}'.format(len(scatters_2x_d1_out), 'scatters_2x_d1_out', output_file_2x_d1))
    np.savetxt(output_file_2x_d1, scatters_2x_d1_out, delimiter=',', fmt='%.5f')
    print ('  Saving {} events from data array ({}) to file (CSV format): {}'.format(len(scatters_3x_out), 'scatters_3x_out', output_file_3x))
    np.savetxt(output_file_3x, scatters_3x_out, delimiter=',', fmt='%.5f')

    #  saving combined 2x and 3x data (using file.open in order to append data to file)
    print ('  Saving {} events from data array ({}) to file (CSV format): {}'.format(len(scatters_2x_out), 'scatters_2x_out', output_file_2x_3x))
    ofile = open(output_file_2x_3x, "w")
    np.savetxt(ofile, scatters_2x_out, delimiter=',', fmt='%.5f')
    print ('   - saving {} additional events from data array ({}) to file (CSV format): {}'.format(len(scatters_3x_out), 'scatters_3x_out', output_file_2x_3x))
    np.savetxt(ofile, scatters_3x_out, delimiter=',', fmt='%.5f')
    ofile.close()

    # count final number of double/triple scatters
    num_ss_fi = len(scatters_1x_out)
    num_ds_fi = len(scatters_2x_out);  num_ds_d0_fi = len(scatters_2x_d0_out)
    num_ds_d1_fi = len(scatters_2x_d1_out);  num_ts_fi = len(scatters_3x_out)

    #  measuring time required to run code - stop time
    stop = timeit.default_timer()

    #  save run summary to file
    summary_file = output_file.replace(".txt", "").replace(".dat", "").replace(".csv", "") + "-summary.txt"
    f = open(summary_file, "w+")
    f.write ('--- Completed processing / time required {} s'.format(stop - start))
    f.write ('\n-------  SUMMARY  -------')
    f.write ('\n-- Total events: {:.0f} ({} interactions), in file: {}'.format((num_1px + (num_2px/2) + (num_3px/3) + (num_4px/4)), len(allevents_data.index), data_filename))
    f.write ('\n--- Start time (s): {:.1f} | End time (s): {:.1f} | Run duration (s): {:.1f} ({:.3f} hrs) | Run activity (events/s): {:.1f}'.format(time_start, time_end, (time_end - time_start), (time_end - time_start)/3600, len(allevents_data.index) / (time_end - time_start)))
    f.write ('\n-- Number of 1x pixel interactions: {} (D0: {}, D1: {}), 2x pixel interactions: {} ({:.1f} double scatter events), 3x pixel interactions: {} ({:.1f} triple scatter events), 4x+ pixel interactions: {}'.format(num_1px, num_1px_d0, num_1px_d1, num_2px, num_2px/2, num_3px, num_3px/3, num_4px))
    f.write ('\n--- Detector transformations: D0 -> {} / D1 -> {}'.format(transforms[0], transforms[1]))
    f.write ('\n-- Run Options -> Big Energy Sort: {} / Apply Energy Threshold: {} / Remove Unphysical Events: {} / Compton Line Filtering: {} /  Shuffle Events: {} / Multistage Coincidence: {}'.format(APPLY_BIG_ENERGY_SORTING, APPLY_ENERGY_THRESHOLD, REMOVE_UNPHYSICAL_EVENTS, APPLY_COMPTON_LINE_FILTERING, SHUFFLE_SCATTER_EVENTS, APPLY_COINCIDENCE_GROUPING_MULTISTAGE))
    if APPLY_COMPTON_LINE_FILTERING: f.write ('\n--- Compton line filter settings -> Range [min, max]: {} / Energies (MeV): {}'.format(COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES))
    if APPLY_COINCIDENCE_GROUPING_MULTISTAGE: f.write ('\n--- Multistage coincidence settings -> Range [min, max]: {} / Energies (MeV): {}'.format(COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES))
    f.write ('\n-- Single Scatters Events       -> Initial: {} | Above Eng Thres: {} | Final: {}'.format(num_ss, num_ss_et, num_ss_fi))
    f.write ('\n-- Double Scatters Events       -> Initial: {} | Above Eng Thres: {} | Physically valid: {} | CL filtered: {} | MS coincidence: {} | Final: {} ({:.1f}%)'.format(num_ds, num_ds_et, num_ds_pe, num_ds_cl, num_ds_ms, num_ds_fi, ((float(num_ds_fi)/num_ds) * 100.0)))
    f.write ('\n-- Double Scatters Events in D0 -> Initial: {} | Above Eng Thres: {} | Physically valid: {} | CL filtered: {} | MS coincidence: {} | Final: {} ({:.1f}%)'.format(num_ds_d0, num_ds_d0_et, num_ds_d0_pe, num_ds_d0_cl, num_ds_d0_ms, num_ds_d0_fi, ((float(num_ds_d0_fi)/num_ds_d0) * 100.0)))
    f.write ('\n-- Double Scatters Events in D1 -> Initial: {} | Above Eng Thres: {} | Physically valid: {} | CL filtered: {} | MS coincidence: {} | Final: {} ({:.1f}%)'.format(num_ds_d1, num_ds_d1_et, num_ds_d1_pe, num_ds_d1_cl, num_ds_d1_ms, num_ds_d1_fi, ((float(num_ds_d1_fi)/num_ds_d1) * 100.0)))
    f.write ('\n-- Triple Scatters Events       -> Initial: {} | Above Eng Thres: {} | Final: {}'.format(num_ts, num_ts_et, num_ts_fi))
    f.write ('\n--------------------------')
    f.close()

    #  print out run summary to screen
    print ('\n--- Completed processing / time required {} s'.format(stop - start))
    print ('\n-------  SUMMARY  -------')
    print ('-- Total events: {:.0f} ({} interactions), in file: {}'.format((num_1px + (num_2px/2) + (num_3px/3) + (num_4px/4)), len(allevents_data.index), data_filename))
    print ('--- Start time (s): {:.1f} | End time (s): {:.1f} | Run duration (s): {:.1f} ({:.3f} hrs) | Run activity (events/s): {:.1f}'.format(time_start, time_end, (time_end - time_start), (time_end - time_start)/3600, len(allevents_data.index) / (time_end - time_start)))
    print ('-- Number of 1x pixel interactions: {} (D0: {}, D1: {}), 2x pixel interactions: {} ({:.1f} double scatter events), 3x pixel interactions: {} ({:.1f} triple scatter events), 4x+ pixel interactions: {}'.format(num_1px, num_1px_d0, num_1px_d1, num_2px, num_2px/2, num_3px, num_3px/3, num_4px))
    print ('--- Detector transformations: D0 -> {} / D1 -> {}'.format(transforms[0], transforms[1]))
    print ('-- Run Options -> Big Energy Sort: {} / Apply Energy Threshold: {} / Remove Unphysical Events: {} / Compton Line Filtering: {} / Shuffle Events: {} / Multistage Coincidence: {}'.format(APPLY_BIG_ENERGY_SORTING, APPLY_ENERGY_THRESHOLD, REMOVE_UNPHYSICAL_EVENTS, APPLY_COMPTON_LINE_FILTERING, SHUFFLE_SCATTER_EVENTS, APPLY_COINCIDENCE_GROUPING_MULTISTAGE))
    if APPLY_COMPTON_LINE_FILTERING: print ('--- Compton line filter settings -> Range [min, max]: {} / Energies (MeV): {}'.format(COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES))
    if APPLY_COINCIDENCE_GROUPING_MULTISTAGE: print ('--- Multistage coincidence settings -> Range [min, max]: {} / Energies (MeV): {}'.format(COMPTON_LINE_RANGE, COMPTON_LINE_ENERGIES))
    print ('-- Single Scatters Events       -> Initial: {} | Above Eng Thres: {} | Final: {}'.format(num_ss, num_ss_et, num_ss_fi))
    print ('-- Double Scatters Events       -> Initial: {} | Above Eng Thres: {} | Physically valid: {} | CL filtered: {} | MS coincidence: {} | Final: {} ({:.1f}%)'.format(num_ds, num_ds_et, num_ds_pe, num_ds_cl, num_ds_ms, num_ds_fi, ((float(num_ds_fi)/num_ds) * 100.0)))
    print ('-- Double Scatters Events in D0 -> Initial: {} | Above Eng Thres: {} | Physically valid: {} | CL filtered: {} | MS coincidence: {} | Final: {} ({:.1f}%)'.format(num_ds_d0, num_ds_d0_et, num_ds_d0_pe, num_ds_d0_cl, num_ds_d0_ms, num_ds_d0_fi, ((float(num_ds_d0_fi)/num_ds_d0) * 100.0)))
    print ('-- Double Scatters Events in D1 -> Initial: {} | Above Eng Thres: {} | Physically valid: {} | CL filtered: {} | MS coincidence: {} | Final: {} ({:.1f}%)'.format(num_ds_d1, num_ds_d1_et, num_ds_d1_pe, num_ds_d1_cl, num_ds_d1_ms, num_ds_d1_fi, ((float(num_ds_d1_fi)/num_ds_d1) * 100.0)))
    print ('-- Triple Scatters Events       -> Initial: {} | Above Eng Thres: {} | Final: {}'.format(num_ts, num_ts_et, num_ts_fi))
    print ('--------------------------')


# Read in the all_events data, perform coordinate transformation and save in AllEventsCombined format
def transform_raw_data(data_dir, data_filename, transforms):
    '''
    Takes in an AllEventsCombined.txt file, reads the information into a pandas
    array, transforms the coordinates, and then saves as txt file (with same format as AllEventsCombined.txt)
    NOTE: The ordering the 2x and 3x scatters may change (due to sorting in the utils.apply_transformation function)

    @params:
        data_dir            - Required : Path to combined detector data .txt file (Str)
        data_filename       - Required : Name of combined detector data .txt file to read data from (Str)
        transforms          - Required : Array of transformation arrays, used for coordinate transformations (Float[][])
                                         format: (rot[0], rot[1], rot[2], pos[0], pos[1], pos[2]) -> units are degrees and mm

    @returns:
        scatters_2x_out     - numpy array of the processed double scatter data recorded by both polaris detectors
        scatters_2x_d0_out  - numpy array of the processed double scatter data recorded by polaris detector 0
        scatters_2x_d1_out  - numpy array of the processed double scatter data recorded by polaris detector 1
        scatters_3x_out     - numpy array of the processed triple scatter data recorded by both polaris detectors
    '''

    #  measuring time required to run code - start time
    start = timeit.default_timer()
    data_filename_transformed = data_filename.replace(".txt", "").replace(".dat", "").replace(".csv", "") + "_transformed.txt"
    print('\n  Transforming interactions position in ' + data_filename + ' and creating ' + data_filename_transformed)

    allevents_data = read_polaris_data(data_dir, data_filename, transforms)       # read in the raw polaris data

    # save transformed pandas array allevents_data as AllEventsCombined_transformed.txt
    #   - allevents_data format is: ['scatters', 'detector', 'energy', 'x', 'y', 'z', 'time']
    #   - Each line in the .txt file has the following format:
    #       Number of interactions
    #       Detector number
    #       Energy, 1 decimal place, (keV)
    #       x position (mm), 2 decimal places
    #       y position (mm), 2 decimal places
    #       z position (mm), 2 decimal places
    #       Time stamp (10*nanoseconds)
    #   Where each value is separated by a tab (\t)

    print('\n  Writing transformed data to ' + data_dir + data_filename_transformed)
    print ('   - Number of interactions: {}'.format(len(allevents_data.index)))

    # Write the final results to file
    with open(data_dir + data_filename_transformed, 'w') as write_file:

        # iterate through each row in pandas array
        for index, row in allevents_data.iterrows():

            write_line = "{:.0f}\t{:.0f}\t{:.1f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.0f}".format(
                row["scatters"],
                row["detector"],
                row["energy"],
                row["x"],
                row["y"],
                row["z"],
                row["time"]
            )
            # print(write_line)
            write_file.write(write_line + "\n")

    print("  Write complete.")

    #  measuring time required to run code - stop time
    stop = timeit.default_timer()

    #  print out run summary to screen
    print ('\n--- Completed processing / time required {} s'.format(stop - start))
    print ('--------------------------')



#------------------------------------------------------------------
# SETTINGS / PARAMETERS
#------------------------------------------------------------------

### Plotting
PRODUCE_RAW_PLOTS = True             # produces 1D & 2D position and energy plots from raw data, produces plots for individual detectors
PRODUCE_TRANSFORMED_PLOTS = False     # produces 1D & 2D position and energy plots from coordinate-transformed data, also produces plots for each scatter number
PRODUCE_FINAL_DOUBLE_PLOTS = True    # produces 1D energy plots of filtered double events, also produced Compton Line plots (E1 vs Theta1)

### Data filtering
APPLY_BIG_ENERGY_SORTING = False
REMOVE_UNPHYSICAL_EVENTS = False
APPLY_ENERGY_THRESHOLD = False       # default value is True
MINIMUM_ENERGY_THRESHOLD = 0.05      # any events below this energy (in MeV) will be removed / default is 0.05 MeV
APPLY_COMPTON_LINE_FILTERING = False
COMPTON_LINE_RANGE = np.array( [0.960, 1.040] )     # min / max values (percentage) for Compton line filtering
COMPTON_LINE_ENERGIES = np.array( [0.5110, 1.274] )  # multiple values possible, i.e. np.array( [0.5110, 0.6617, 1.173, 1.274, 1.332] )
SHUFFLE_SCATTER_EVENTS = False

### Coincidence grouping
APPLY_COINCIDENCE_GROUPING_MULTISTAGE = False

#------------------------------------------------------------------
# MAIN PROGRAM
#------------------------------------------------------------------

def main():
    print ('\n--- Running polarisPGI_processing.py ---')

    """
    #### 180315 Data ####
    # ----------------- #
    #  NOTE: Data shuffling required for runs 6, 7, 8 (problem with timing between detectors)
    d0_transformation = [180.0, 90.0, 0.0, -92.0, 0.0, 0.0]  # all runs
    d1_transformation = [90.0, 0.0, 90.0, 0.0, 94.0, 0.0]    # all runs
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run1-co60_at_15_-05_11mm-oem2_x_-92mm-oem3_y_94mm-3hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run2-co60_at_05_-05_11mm-oem2_x_-92mm-oem3_y_94mm-3hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run3-co60_at_-05_-05_11mm-oem2_x_-92mm-oem3_y_94mm-3hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run4-co60_at_-15_-05_11mm-oem2_x_-92mm-oem3_y_94mm-3hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run5-cs137_at_-15_-05_11mm-oem2_x_-92mm-oem3_y_94mm-2hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run6-cs137_at_-05_-05_11mm-oem2_x_-92mm-oem3_y_94mm-2hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run7-cs137_at_05_-05_11mm-oem2_x_-92mm-oem3_y_94mm-2hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-prestream/180315/J_180315-run8-cs137_at_15_-05_11mm-oem2_x_-92mm-oem3_y_94mm-2hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    """
    """
    #### 190516 Data ####
    # ----------------- #
    # NOTE: D0 and D1 were swapped from the usual configuration
    d0_transformation = [90.0, 0.0, 90.0, 0.0, 94.0, 0.0]    # all runs
    d1_transformation = [180.0, 90.0, 0.0, -93.0, 0.0, 0.0]  # all runs
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190516/190516-run1-co60_at_-3_6_11mm-oem2_y_93mm-oem3_x_-94mm-6hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190516/190516-run3-na22_at_-3_6_11mm-oem2_y_93mm-oem3_x_-94mm-6hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #transform_raw_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190516/190516-run6-co60_at_-3_6_11mm-oem2_y_93mm-oem3_x_-94mm-900s/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #transform_raw_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190516/190516-run2-cs137_at_-3_6_11mm-oem2_y_93mm-oem3_x_-94mm-6hrs/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    """
    """
    #### 190711 Data ####
    # ----------------- #
    d0_transformation = [180.0, 90.0, 0.0, -45.5, 0.0, 0.0]  # run 1
    d1_transformation = [0.0, 270.0, 0.0, 45.5, 0.0, 0.0]    # run 1
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190711/190711-run1-na22_oem2_y_+48_oem3_y_-48mm-z_trav/tina_source_at_centre_data/', 'AllEventsCombined_0mm.txt', [d0_transformation, d1_transformation])
    """
    """
    #### 190805 Data ####
    # ----------------- #
    d0_transformation = [180.0, 90.0, 0.0, -48, 0.0, 0.0]  # run 1
    d1_transformation = [0.0, 270.0, 0.0, 48, 0.0, 0.0]    # run 1
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190805/191004-run1-neutron_beam_on_d1-noSubMmPixel-15min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190805/190805-run2-na22_at_-1_0_11mm-oem2_x_-48mm-oem3_x_48mm-NOsubmm-10min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190805/190805-run3-na22_at_-1_0_11mm-oem2_x_-48mm-oem3_x_48mm-submm-10min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/190805/190805-run4-na22_at_-1_0_11mm-oem2_x_-48mm-oem3_x_48mm-NOsubmm-10min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    """
    """
    #### 191004 Data ####
    # ----------------- #
    #d0_transformation = [180.0, 90.0, 0.0, +475.0, -338.0, -11.0]  # runs 1-3
    #d1_transformation = [180.0, 90.0, 0.0, +475.0,  +12.0, -11.0]  # runs 1-3
    d0_transformation = [0.0, 270.0, 0.0, +437.0, -338.0, -11.0]   # runs 4-5
    d1_transformation = [0.0, 270.0, 0.0, +437.0,  +12.0, -11.0]   # runs 4-5
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/191004/191004-run1-neutron_beam_on_d1-noSubMmPixel-15min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/191004/191004-run5-neutron_beam_on_d1-withSubMmPixel-15min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    """
    """
    #### 200724 Data ####
    # ----------------- #
    d0_transformation = [180.0, 90.0, 0.0, -42.5, 0.0, 0.0]  # all runs
    d1_transformation = [0.0, 270.0, 0.0, 42.5, 0.0, 0.0]    # all runs
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/200724/200724-run1-na22_at_0_2.5_11mm-oem2_x_-40.5mm-oem3_x_40.5mm-submm-5min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    #filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/200724/200724-run3-na22_at_0_2.5_11mm-oem2_x_-40.5mm-oem3_x_40.5mm-nosubmm-5min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    filter_data('/Volumes/Batman/PGI/data/uct/polaris-streaming/200724/200724-run4-na22_at_0_2.5_11mm-oem2_x_-40.5mm-oem3_x_40.5mm-nosubmm-source_insert-30sec/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])
    """

    #### 191203 Data ####
    # ----------------- #
    d0_transformation = [180.0, 90.0, 0.0, -41.23, 0.0, 0.0]  # all runs
    d1_transformation = [0.0, 270.0, 0.0, 41.23, 0.0, 0.0]    # all runs
    transform_raw_data('/home/steve/PGI/data/uct/polaris-streaming/191203/191203-run19-ga68_at_0_0_11mm-oem2_x_-41mm-oem3_x_41mm-20min/', 'AllEventsCombined.txt', [d0_transformation, d1_transformation])


if __name__ == "__main__":

    main()
