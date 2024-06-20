#!/usr/bin/python

# 210721 - Added to force use of python3 verion of Cython
#import include
# cython: language_level=3

"""
Description:
Contains functions that will be used with processPolarisData.py for Compton camera event processing.

NOTES:
- written using python v2.7.15 / updated to also work in python 3.6.2

Code borrowed heavily from Dennis Mackin <dsmackin@mdanderson.org>
"""
__author__ = "Steve Peterson <steve.peterson@uct.ac.za>"
__date__ = "November 01, 2018"
__version__ = "$Revision: 4.0.0$"

#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------
import sys, os
import numpy
cimport numpy as np
from math import sin, cos, pi, log, floor
import pandas
import cProfile
import re
import matplotlib
import matplotlib.pyplot as plt
import pylab


#------------------------------------------------------------------
# CONSTANTS
#------------------------------------------------------------------
MeCsq = 0.5109989461  # electron mass in energy units (MeV)


#------------------------------------------------------------------
# FUNCTION/CLASS DEFINITIONS
#------------------------------------------------------------------

#  takes 3D rotation and 3D translation and produces transformation matrix, input is an array
def get_transformation_matrix_array(TM):

    # break up array into individual elements -> units are degrees and mm
    rot_x, rot_y, rot_z = TM[0], TM[1], TM[2]
    trans_x, trans_y, trans_z = TM[3], TM[4], TM[5]

    # convert rotations into radians
    rot_x *= pi/180.0
    rot_y *= pi/180.0
    rot_z *= pi/180.0

    Rx = numpy.matrix([
            [1, 0, 0, 0],
            [0, cos(rot_x), sin(rot_x), 0],
            [0, -sin(rot_x), cos(rot_x), 0],
            [0, 0, 0, 1]
            ]
    )

    Ry = numpy.matrix([
            [cos(rot_y), 0, -sin(rot_y), 0],
            [0, 1, 0, 0],
            [sin(rot_y), 0, cos(rot_y), 0],
            [0, 0, 0, 1]
            ]
    )

    Rz = numpy.matrix([
            [cos(rot_z), sin(rot_z), 0, 0],
            [-sin(rot_z), cos(rot_z), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
            ]
    )

    Trans = numpy.matrix(
        [
            [1, 0, 0, -trans_x],
            [0, 1, 0, -trans_y],
            [0, 0, 1, -trans_z],
            [0, 0, 0, 1]
        ]
    )

    return Rz*Ry*Rx*Trans


#  takes 3D rotation and 3D translation and produces transformation matrix, input is 6 values
def get_transformation_matrix(rot_x, rot_y, rot_z, trans_x, trans_y, trans_z):
    rot_x *= pi/180.0
    rot_y *= pi/180.0
    rot_z *= pi/180.0

    Rx = numpy.matrix([
            [1, 0, 0, 0],
            [0, cos(rot_x), sin(rot_x), 0],
            [0, -sin(rot_x), cos(rot_x), 0],
            [0, 0, 0, 1]
            ]
    )

    Ry = numpy.matrix([
            [cos(rot_y), 0, -sin(rot_y), 0],
            [0, 1, 0, 0],
            [sin(rot_y), 0, cos(rot_y), 0],
            [0, 0, 0, 1]
            ]
    )

    Rz = numpy.matrix([
            [cos(rot_z), sin(rot_z), 0, 0],
            [-sin(rot_z), cos(rot_z), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
            ]
    )

    Trans = numpy.matrix(
        [
            [1, 0, 0, -trans_x],
            [0, 1, 0, -trans_y],
            [0, 0, 1, -trans_z],
            [0, 0, 0, 1]
        ]
    )

    return Rz*Ry*Rx*Trans


#  transforms xyz coordinates
def apply_transformation(df, transformation_matrices):

    #  takes df.detector values and turns into a set of possible values, i.e. [0, 1] if two detectors
    detector_numbers = set(df.detector)
    #  error checking - if the number of detector_numbers is greater than the number of transform_matrices, code will break
    assert(len(detector_numbers) <= len(transformation_matrices))

    #  splits df data into individual groups for each detector
    #   - to_numpy() -> Convert the frame to its Numpy-array representation
    #   - also re-ordering events based detector number [destroying ordering of events] <- re-sort by timestamp at the end
    detector_slices = [df[df['detector'] == n].to_numpy() for n in detector_numbers]
    #  "cdef" re-assigns the variable xyz is a cython array in order to gain the speed-up benefit
    cdef np.ndarray[np.float_t, ndim=2] xyz

    #  loop through detector_slices [which have been rearranged by detector number]
    for i, n in enumerate(detector_numbers):
        xyz = detector_slices[i][:, 3:6]
        #xyz = detector_slices[i][:, 4:7]

        #  the POLARIS detectors have left handed Coordinate system (so we first flip y-axis)
        xyz[:,1] *= -1  # flip y-axis
        print ('   - detector {} - converted to right handed coordinates (flip y-axis)'.format(i))
        xyz = xyz.T
        xyz_four_length = numpy.append(xyz, [numpy.ones(xyz.shape[1])], axis=0)

        #  the detector to the reconstruction transformation is defined as the inverse transformation
        M = numpy.linalg.inv(transformation_matrices[n])
        xyz_transformed = numpy.dot(M, xyz_four_length)
        detector_slices[i][:, 3:6] = xyz_transformed[:-1, :].T
        #detector_slices[i][:, 4:7] = xyz_transformed[:-1, :].T
        print ('   - detector {} - successful coordinate transformation!'.format(i))

    #  recombines transformed coordinates back into full dataFrame
    new_df = pandas.DataFrame(numpy.concatenate(detector_slices), columns=df.columns)

    #  re-sort data by time stamp [important for coincidence grouping]
    new_df = new_df.sort_values(by = ['time', 'y'], ascending = [True, True])
    #  not sure if this is important
    new_df = new_df.reset_index(drop = True)

    return new_df


#  returns interaction data (E, X, Y, Z) for the specified number of scatters and detector number
def get_interaction_data(df, energy_sort = False, num_scatters = 2, det_num = None):

    print ('   - Grouping data for detector {} with scatter number {}'.format(det_num, num_scatters))

    #  stores the data with the relevant number of scatters & detector number
    if det_num is not None:
        chosen = df.loc[numpy.isclose(df.scatters, num_scatters) & numpy.isclose(df.detector, det_num)]
    else:
        chosen = df[numpy.isclose(df.scatters, num_scatters)]

    #  sorts the scatter of each event by energy deposited, highest first
    if energy_sort:
        chosen = chosen.sort_values(by = ['time', 'energy'], ascending = [True, False])

    #  re-formatting scatter data
    #   - to_numpy() -> Convert the frame to its Numpy-array representation
    scatters = chosen.to_numpy()
    #  only store columns 2, 3, 4, 5 (i.e. eng, x, y, z)
    scatters = scatters[:, [2, 3, 4, 5]]
    #  FIXED THIS -> only store columns 1, 4, 5, 6 (i.e. eng, x, y, z) <- columns reordered alphabetically by label
    #scatters = scatters[:, [1, 4, 5, 6]]

    #  checks that the final number of event is evenly divisible by the number of scatters <- disabled warning
    if len(scatters) % num_scatters != 0:
        # Change Exception to Warning
        #raise Exception("Wrong number of scatters (%d) for scatters (%d), mod = (%d)" \
        #                % (len(scatters), num_scatters, len(scatters) % num_scatters))
        print ("!!! WARNING: Wrong number of scatters (%d) for scatters (%d), mod = (%d)" % (len(scatters), num_scatters, len(scatters) % num_scatters))

    #  print counts to screen
    print ('    - found {} events out of {}'.format(len(scatters), len(df)))

    return scatters


#  rescales energy column from keV to MeV
def scale_energy(scatters):
    '''
    Args:
        scatters in 2D array with rows energy (keV), x, y, z

    Returns:
        scatters in 2D array with rows energy (MeV), x, y, z
    '''
    if numpy.average(scatters[:,0]) > 10:
        scatters[:,0] = scatters[:,0] * 0.001

    return scatters


#  re-arranges scatter data from consecutive lines to all Compton scatter data on a single line for CSV output / also re-scales energy
def format_data_for_output(scatters, num_scatters, det_num = None):

    print ('   - Formatting {} events from detector {} and number of scatters: {}'.format(len(scatters), det_num, num_scatters))

    #  check if any data in scatters array
    if scatters.shape[0] == 0:
        return scatters

    #  re-ordering interactions [CURRENTLY EMPTY FUNCTION]
    #scatters = put_interactions_in_order(scatters, num_scatters)
    #  converting energy values from keV into MeV
    scatters = scale_energy(scatters)

    #  re-arranges double scatters into E1, X1, Y1, Z1, E2, X2, Y2, Z2 format
    if num_scatters == 2:
        #  moves next event (i + 1) to the same line as current event
        temp = numpy.concatenate((scatters[0:-1], scatters[1:]), axis = 1)
        #  stores data from every other line
        out = temp[0::2]

    #  re-arranges triple scatters into E1, X1, Y1, Z1, E2, X2, Y2, Z2, E3, X3, Y3, Z3 format
    if num_scatters == 3:
        #  moves next two events (i + 1 & i + 2) to the same line as current event
        temp = numpy.concatenate((scatters[0:-2], scatters[1:-1], scatters[2:]), axis = 1)
        #  stores data from every third line
        out = temp[0::3]

    return out


#  plots 2D histogram
def plot_2D(np.ndarray[np.float_t, ndim=1] x, x_label, np.ndarray[np.float_t, ndim=1] y, y_label, title, output_folder=".", is_log=False):

    print ('  Creating 2D Plot - {} v {} with title: {} . . .'.format(x_label, y_label, title))
    plt.clf()
    ax = plt.subplot(111)
    ax.set_title(title)

    nbins = int(max(max(x) - min(x), max(y) - min(y)))
    if nbins < 1: nbins = 1

    H, xedges, yedges = numpy.histogram2d(y, x, bins=(nbins, nbins))
    if(is_log): H = numpy.log(H)
    im = plt.imshow(H, interpolation='none', origin='lower', extent=[min(yedges), max(yedges), min(xedges), max(xedges)])

    plt.title(title + ', total = {}'.format(len(x)), fontsize = 18)
    ax.xaxis.set_tick_params(labeltop = 'off')
    ax.xaxis.set_tick_params(labelbottom = 'on')

    ax.grid(True,linestyle='-',which='both', color='0.50')

    plt.xlabel(x_label, fontsize = 18)
    plt.ylabel(y_label, fontsize = 18)
    ax.xaxis.set_label_position('bottom')
    ax.tick_params(axis='both', bottom = 'on', top = 'off')

    plot_name = "%s/%s_%s_%s.png" % (os.path.abspath(output_folder), title.replace(":", "_").replace(" ", "_"), x_label, y_label)
    plt.savefig(plot_name)

    return plt


#  manually bins data
def bin_data(np.ndarray[np.float_t, ndim=1] x, int num_bins):

    cdef float range_min = min(x)
    cdef float range_max = max(x) + 1.0E-2 # make bin edge bigger than largest x.
    cdef float range_length = range_max - range_min
    print ('    - plot details: range_max = {}, range_min = {}, range_length = {}'.format(range_max, range_min, range_length))
    if (range_length * range_length < 1.0E-6):
        print ("ERROR: Range length is 0 (min=%.3e, max=%.3e, range=%.3e) . .. " % (range_min, range_max, range_length))

    cdef float step_factor = float(num_bins - 1)/range_length

    cdef np.ndarray[np.float_t, ndim=1] xshifted = x - range_min
    cdef np.ndarray[np.float_t, ndim=1] x_bin_number = xshifted * step_factor

    bins = numpy.linspace(0, num_bins - 1, num_bins)
    hist_dict = dict(zip(bins, numpy.zeros(len(bins))))
    hist_bin_values = bins

    def index_dict(int key):
        hist_dict[key] += 1

    map(index_dict, x_bin_number)

    cdef np.ndarray[np.float_t, ndim=1] x_vals = bins * range_length/num_bins + range_min
    cdef np.ndarray[np.float_t, ndim=1] y_vals = numpy.array([float(hist_dict[key]) for key in bins])

    return x_vals, y_vals


#  plots 1D profiles, uses manually binned data
def plot_1D(x, num_bins, x_label, y_label, title, output_folder = ".", isLog = False):
    '''
    The builtin python histogram function choke if the number of values gets too large.
    This function has a built in work around. It handles the binning and counting itself
    and then plots an x, y line plot.
    '''

    print ('  Creating 1D Plot - {} v {} with title: {} . . .'.format(x_label, y_label, title))
    plt.clf()
    plt.gcf().set_size_inches(8, 8)
    ax = plt.subplot(111)

    try:
        x = x.to_numpy()
    except AttributeError:
        True

    # not using this (doesn't seem to work in python3)
    """
    x_vals, y_vals = bin_data(x, num_bins)
    if (isLog): plt.hist(x_vals, num_bins, weights = y_vals, log = True)
    #else: plt.hist(x_vals, num_bins, weights = y_vals)
    """
    # using this instead
    if (isLog): plt.hist(x, num_bins, log = True)
    else: plt.hist(x, num_bins)

    ax.set_title(title + ', total = {}'.format(len(x)), fontsize = 18)
    plt.xlabel(x_label, fontsize = 18)
    plt.ylabel(y_label, fontsize = 18)

    ax.grid(True, linestyle='-', which='both', color='0.750')

    if (isLog): plot_name = "%s/%s_%s_%s_log.png" % (os.path.abspath(output_folder), title.replace(":", "_").replace(" ", "_"), x_label, y_label)
    else: plot_name = "%s/%s_%s_%s.png" % (os.path.abspath(output_folder), title.replace(":", "_").replace(" ", "_"), x_label, y_label)
    plt.savefig(plot_name)

    return plt


#  function to create series of plots from basic detector data
def make_basic_detector_plots(df, moniker, output_folder):

    # check for data in array (if array is full, will return True)
    if not df.empty:

        #  creates list of unique detector indices
        #    pandas.unique - returns unique values of the Series object and are returned in order of appearance
        detectors = pandas.unique(df.detector)
        print ('  set of detectors: {}'.format(detectors))

        #  loop through each detector
        for detector in detectors:
            #  pull out data for each detector
            df_detector = df[df.detector == detector]
            print ('    - detector {} / label: {}'.format(detector, moniker))
            #  plot 1D profiles for each detector
            plot_1D(df_detector.x.to_numpy(), 200, "x", "Counts", "%s_D%d" % (moniker, detector), output_folder)
            plot_1D(df_detector.y.to_numpy(), 200, "y", "Counts", "%s_D%d" % (moniker, detector), output_folder)
            plot_1D(df_detector.z.to_numpy(), 200, "z", "Counts", "%s_D%d" % (moniker, detector), output_folder)
            #  plot 2D profiles for each detector
            plot_2D(df_detector.y.to_numpy(), "y", df_detector.x.to_numpy(), "x", "%s_D%d" % (moniker, detector), output_folder, True)
            plot_2D(df_detector.y.to_numpy(), "y", df_detector.z.to_numpy(), "z", "%s_D%d" % (moniker, detector), output_folder, True)
            plot_2D(df_detector.z.to_numpy(), "z", df_detector.x.to_numpy(), "x", "%s_D%d" % (moniker, detector), output_folder, True)

        #  plot energy for all detectors
        plot_1D(df.energy.to_numpy(), 200, "Energy", "Counts", "%s" % moniker, output_folder, False)
        #  plot energy for all detectors (log)
        plot_1D(df.energy.to_numpy(), 200, "Energy", "Counts", "%s" % moniker, output_folder, True)


#  function used by filtering_unphysical_double_scatters() to find unphysical events
#   - calculate the inner product of Compton equation to find theta1 from E0 and E1
def physical_energy_ordering_double(E1, E2):

    #  check if energy ordering of double scatter event produces physical scatter angle
    E0 = E1 + E2

    if numpy.abs(1 + MeCsq * ( 1.0/(E0) - 1.0/(E0 - E1) ) ) < 1:
        return 1
    else:
        return 0


#  checks for unphysical events (theta1 = nan) from double scatter event data
#   - input format (scatters): eng1, x1, y1, z1, eng2, x2, y2, z2
#   - returns data in same format
def filtering_unphysical_double_scatters(scatters):

    scatters_physical = []
    count_both, count_one, count_flip, count_none = [0] * 4

    #  loop through list of events
    for index in range( len(scatters) ):

        #  check original energy ordering
        order1 = physical_energy_ordering_double(scatters[index][0], scatters[index][4])

        #  check flipped energy ordering
        order2 = physical_energy_ordering_double(scatters[index][4], scatters[index][0])

        #  storing appropriate events into final output array: scatters_physical
        if (order1 == 1 & order2 == 1):
            scatters_physical.append(scatters[index])
            count_both += 1
        elif (order1 == 1):
            scatters_physical.append(scatters[index])
            count_one += 1
        elif (order2 == 1):
            SE = scatters[index]
            flipped_scatter = numpy.array((SE[4], SE[5], SE[6], SE[7], SE[0], SE[1], SE[2], SE[3]))
            scatters_physical.append(flipped_scatter)
            count_one += 1; count_flip += 1
        else:
            count_none += 1

    #  print counts to screen
    print ('    - results of event ordering -> number of physical events returned: {} | both work: {} | only one order works: {} | order flipped: {} | neither work: {}'.format(len(scatters_physical), count_both, count_one, count_flip, count_none))

    #  return physical events
    #    vstack takes list of numpy arrays and converts into single numpy array
    return numpy.vstack(scatters_physical)



#  function used by compton_line_filtering() to calculate the energy of the first scatter
def calculate_expected_first_scatter_energy(E0 , theta):

    #  re-arrangement for Compton scatter equation to solve for E1
    alpha = E0 / MeCsq
    beta = alpha * ( 1 - numpy.cos(theta) )
    return E0 * beta / ( 1 + beta )


#  filtering events using Compton Line Filtering (based on expected gamma energies, input variable: CL_energies)
#   - input format (scatters): eng1, x1, y1, z1, eng2, x2, y2, z2
#   - returns data in same format
def compton_line_filtering(scatters, CL_range, CL_energies):

    scatters_filtered = []
    count_compton = [0] * len(CL_energies)

    #  loop through list of events
    for index in range( len(scatters) ):

        #  use the first and second energies to calculate E0 and theta1
        E1 = scatters[index][0]; E2 = scatters[index][4]
        E0 = E1 + E2
        theta1 = numpy.arccos( 1 + MeCsq * ( 1.0/(E0) - 1.0/(E0-E1) ) );

        #  looping through expected gamma energies
        for i, e in enumerate(CL_energies):
            #  calculating range of Compton values for given total energy and first scatter angle
            minE1 = CL_range[0] * calculate_expected_first_scatter_energy(e , theta1)
            maxE1 = CL_range[1] * calculate_expected_first_scatter_energy(e , theta1)

            #  check if energy falls within expected range
            if E1 > minE1 and E1 < maxE1:

                #  add filtered data to output array: scatters_filtered
                scatters_filtered.append(scatters[index])
                count_compton[i] += 1

    #  total up number of Compton filtered events
    total_compton = numpy.sum(count_compton)

    #  print results to screen
    print ('    - results of Compton filtering -> number of filtered events returned: {} -> counts: {} by energy {}, respectively'.format(total_compton, count_compton, CL_energies))

    #  return filtered events
    #    vstack takes list of numpy arrays and converts into single numpy array
    return numpy.vstack(scatters_filtered)



def create_plot_range(plot_data, x_label, nb_bins, rng_min, rng_max):
    plt.hist(plot_data, bins=nb_bins, range=[rng_min,rng_max], alpha = 1.0)
    plt.xlabel(x_label);
    plt.ylabel('Number of Gammas');


def create_plot_range_log(plot_data, x_label, nb_bins, rng_min, rng_max):
    plt.figure(figsize=(8, 4)) # (width, height in inches)
    plt.hist(plot_data, bins=nb_bins, range=[rng_min,rng_max], alpha = 1.0, log=True)
    plt.xlabel(x_label);
    plt.ylabel('Number of Gammas');


def create_plot_bins(plot_data, x_label, nb_bins):
    plt.hist(plot_data, bins=nb_bins)
    plt.xlabel(x_label);
    plt.ylabel('Number of Gammas');


# source: https://tomspur.blogspot.co.za/2015/08/publication-ready-figures-with.html
def save_to_file(filename, dir, fig=None):
    """Save to @filename with a custom set of file formats.

    By default, this function takes to most recent figure,
    but a @fig can also be passed to this function as an argument.
    """
    formats = [
                #"pdf",
                #"eps",
                "png",
                #"svg",
                #"pgf",   # ERROR: LatexError: LaTeX returned an error, probably missing font or error in preamble:
              ]
    if fig is None:
        for form in formats:
            plt.savefig("%s/%s.%s"%(dir, filename, form), dpi=300)
    else:
        for form in formats:
            fig.savefig("%s/%s.%s"%(dir, filename, form), dpi=300)

#  code from paul maggi [modified only in formatting]
def coincCheckMod(timeStamps, ene, modList, xList, yList, zList, cutOffClk, peakMin, peakMax):
    '''
    deltaT, mod1, mod2, doubs = coincCheckMod(timeStamps, ene, modList, xList, yList, zList, cutOffClk, peakMin, peakMax)

    this function should only be passed events that are single pixel events, e.g. ene[npx==1]. For each interaction, it looks up to cutOffClk beyond that point
	for a combination of events that is within the energy windows specified by peakMin and peakMax.

    input:
        timeStamps, ene, modList, xList, yList, zList:
            these are the t, edep, mdl, x, y, and z array from the output of readInC, after only selecting one pixel events
        cutOffClk:
            the number of 10 ns clock cycles beyond a point to look for coincidence events.
				example:
					cutOffSec = 1E-6 %this is a 1 us time window
					cutOffClk = cutOffSec/(10E-9) %this gives a value of 100
    output:
        deltaT:
			time difference (in clk cycles) between the two accepted coincidence points
		mod1, mod2:
			module numbers of the first and second interactions of a coincidence pair, respectively
		doubs:
			list of doubles in the standard format (edep1, x1, y1, z1, edep2, x2, y2, z2)

    '''
    # note: maxPts is the max number of coincidence Singles->doubles this code looks for.
    maxPts = 500000

    deltaT = numpy.zeros((maxPts,))
    mod1 = numpy.zeros((maxPts,))
    mod2 = numpy.zeros((maxPts,))
    yList1 = numpy.zeros((maxPts,))
    yList2 = numpy.zeros((maxPts,))
    ed1 = numpy.zeros((maxPts,))
    ed2 = numpy.zeros((maxPts,))
    xList1 = numpy.zeros((maxPts,))
    xList2 = numpy.zeros((maxPts,))
    zList1 = numpy.zeros((maxPts,))
    zList2 = numpy.zeros((maxPts,))
    jjj = 0
    timeDiff = timeStamps[1:] - timeStamps[0:-1]

    for iii in range(len(timeDiff) - 1):
        checkNum = 1
        buffTime = timeDiff[iii]
        if ene[iii] == 0:
            break
        while (buffTime <= cutOffClk) & ((checkNum + iii) < len(timeStamps) - 3):
            eneSum = ene[iii] + ene[checkNum + iii]
            if (eneSum <= peakMax) & (eneSum >= peakMin):
                deltaT[jjj] = buffTime
                mod1[jjj] = modList[iii] + 1
                mod2[jjj] = modList[checkNum + iii] + 1
                ed1[jjj] = ene[iii]
                ed2[jjj] = ene[iii + checkNum]
                yList1[jjj] = yList[iii]
                yList2[jjj] = yList[checkNum + iii]
                xList1[jjj] = xList[iii]
                xList2[jjj] = xList[checkNum + iii]
                zList1[jjj] = zList[iii]
                zList2[jjj] = zList[checkNum + iii]
                jjj += 1

            checkNum = checkNum + 1
            buffTime = buffTime + timeDiff[iii + checkNum]

        if jjj > maxPts - 2:

            break
    mod1 = numpy.trim_zeros(mod1) - 1
    mod2 = numpy.trim_zeros(mod2) - 1
    doubs = numpy.array((numpy.trim_zeros(ed1), numpy.trim_zeros(xList1), numpy.trim_zeros(yList1), numpy.trim_zeros(zList1), numpy.trim_zeros(ed2), numpy.trim_zeros(xList2), numpy.trim_zeros(yList2), numpy.trim_zeros(zList2))).T

    return numpy.trim_zeros(deltaT), mod1, mod2, doubs


#  updated code to check PET coincidences [modifed from original code by paul maggi]
def coincCheckPET(timeStamps, ene, modList, xList, yList, zList, cutOffClk, peakMin, peakMax):
    '''
    deltaT, mod1, mod2, doubs, eng1, eng2 = coincCheckPET(timeStamps, ene, modList, xList, yList, zList, cutOffClk, peakMin, peakMax)

    this function should only be passed events that are single pixel events, e.g. ene[npx==1]. For each interaction, it looks up to cutOffClk beyond that point
	for a combination of events that are both within the energy window specified by peakMin and peakMax.

    input:
        timeStamps, ene, modList, xList, yList, zList:
            these are the t, edep, mdl, x, y, and z array from the output of readInC, after only selecting one pixel events
        cutOffClk:
            the number of 10 ns clock cycles beyond a point to look for coincidence events.
				example:
					cutOffSec = 1E-6 %this is a 1 us time window
					cutOffClk = cutOffSec/(10E-9) %this gives a value of 100
    output:
        deltaT:
			time difference (in clk cycles) between the two accepted coincidence points
		mod1, mod2:
			module numbers of the first and second interactions of a coincidence pair, respectively
		doubs:
			list of doubles in the standard format (time1, x1, y1, z1, time2, x2, y2, z2)
		eng1, eng2:
			energy deposited in the first and second interactions of a coincidence pair, respectively

    '''
    #  note: maxPts is the max number of coincidence Singles->doubles this code looks for.
    maxPts = 100000

    print ('  -- Length of timeStamps = {} / maxPts = {}'.format(len(timeStamps), maxPts))

    deltaT = numpy.zeros((maxPts,))
    mod1 = numpy.zeros((maxPts,))
    mod2 = numpy.zeros((maxPts,))
    yList1 = numpy.zeros((maxPts,))
    yList2 = numpy.zeros((maxPts,))
    ed1 = numpy.zeros((maxPts,))
    ed2 = numpy.zeros((maxPts,))
    t1 = numpy.zeros((maxPts,))
    t2 = numpy.zeros((maxPts,))
    xList1 = numpy.zeros((maxPts,))
    xList2 = numpy.zeros((maxPts,))
    zList1 = numpy.zeros((maxPts,))
    zList2 = numpy.zeros((maxPts,))
    jjj = 0

    #  timeDiff is the next time step minus the current time stamp [always positive]
    timeDiff = timeStamps[1:] - timeStamps[0:-1]

    #  loop through all of the time steps
    for iii in range(len(timeDiff) - 1):
        checkNum = 1
        buffTime = timeDiff[iii]

        #  kills loop if deposited energy is zero
        if ene[iii] == 0:
            break

        #  loops successive time steps (using checkNum) until outside cutOff time (cutOffClk)
        while (buffTime <= cutOffClk) & ((checkNum + iii) < len(timeStamps) - 3):

            #  saves a data point if the energy value for both time steps (iii & checkNum + iii) are within energy window
            if (ene[iii] <= peakMax) & (ene[iii] >= peakMin) & (ene[checkNum + iii] <= peakMax) & (ene[checkNum + iii] >= peakMin):

                deltaT[jjj] = buffTime
                mod1[jjj] = modList[iii] + 1
                mod2[jjj] = modList[checkNum + iii] + 1
                ed1[jjj] = ene[iii]
                ed2[jjj] = ene[iii + checkNum]
                # convert time into microseconds
                t1[jjj] = timeStamps[iii] * (10E-9) / (1E-6)
                t2[jjj] = timeStamps[iii + checkNum] * (10E-9) / (1E-6)
                yList1[jjj] = yList[iii]
                yList2[jjj] = yList[checkNum + iii]
                xList1[jjj] = xList[iii]
                xList2[jjj] = xList[checkNum + iii]
                zList1[jjj] = zList[iii]
                zList2[jjj] = zList[checkNum + iii]
                jjj += 1

            #  increment checkNum and time difference (buffTime)
            checkNum = checkNum + 1
            buffTime = buffTime + timeDiff[iii + checkNum]

        #  kill loop if max number reached
        if jjj > maxPts - 2:
            break

    #  clean up output arrays (remove zeros from the end)
    mod1 = numpy.trim_zeros(mod1) - 1
    mod2 = numpy.trim_zeros(mod2) - 1
    doubs = numpy.array((numpy.trim_zeros(t1), numpy.trim_zeros(xList1), numpy.trim_zeros(yList1), numpy.trim_zeros(zList1), numpy.trim_zeros(t2), numpy.trim_zeros(xList2), numpy.trim_zeros(yList2), numpy.trim_zeros(zList2))).T
    eng1 = numpy.trim_zeros(ed1)
    eng2 = numpy.trim_zeros(ed2)

    return numpy.trim_zeros(deltaT), mod1, mod2, doubs, eng1, eng2



#  filtering raw data to produce list of back to back coincident 511 keV gammas
#   - function pulls timing/energy settings from main program
def grouping_backtoback_coincidence_events(df, timing, tWindow, eWindow):

    #  converting data from pandas to numpy arrays (not really necessary)
    npx = numpy.array(df['scatters'])
    mdl = numpy.array(df['detector'])
    edep = numpy.array(df['energy'])
    x = numpy.array(df['x'])
    y = numpy.array(df['y'])
    z = numpy.array(df['z'])
    t = numpy.array(df['time'])

    #  calculating number of 10 ns clock cycles for grouping (based on timing setting)
    number_of_clock_cycles = timing/(10E-9)
    print ('   - Using timing window of {} s, thus looking at {} 10 ns clock cycles for coincidences'.format(timing, number_of_clock_cycles))

    #  convert time into seconds (t given in units of 10 ns clock cycles), so tS = t * 10E-9
    tS = t * 10E-9
    tMin = numpy.amin(tWindow[0])
    tMax = numpy.amax(tWindow[1])

    #  filtering out single pixel events (npx == 1) within time window
    checkInds = numpy.logical_and(npx == 1, tS > tMin)
    checkInds = numpy.logical_and(checkInds, tS < tMax)
    #  print time filtered data  to screen
    numFiltered = numpy.sum(checkInds)
    perFiltered = numpy.sum(checkInds) / float(npx.size)
    print ('    - data filter: single pixel events within time window: {0} to {1} s / number of filtered events: {2} ({3:.2f}%)'.format(tMin, tMax, numFiltered, perFiltered * 100.0))

    #  setting energy window (in keV)
    engMin = eWindow[0] * 1000.0
    engMax = eWindow[1] * 1000.0
    #  print energy filtered data  to screen
    engFiltered = numpy.sum(numpy.logical_and(edep > engMin, edep < engMax))
    print ('    - number of events between {} and {} kev is {}'.format(engMin, engMax, engFiltered))

    #  grouping singles into PET coincidences
    #   - looks at single pixel events in different detectors with the given time window with each gamma having a value within the given energy window
    deltaTM, modL1, modL2, time_position_doubles, E1, E2 = coincCheckPET(t[checkInds], edep[checkInds], mdl[checkInds], x[checkInds], y[checkInds], z[checkInds], number_of_clock_cycles, engMin, engMax)

    #  print counts to screen
    print ('    - results -> number of coincidence events returned: {}'.format(deltaTM.size))

    #  return two arrays
    #   - time_position_doubles: time1, x1, y1, z1, time2, x2, y2, z2
    #   - coincidence_data: deltaT, module1, module2, eng1, eng2
    return time_position_doubles, numpy.array( (deltaTM, modL1, modL2, E1, E2) )


#  updated code to check PET coincidences [modifed from original code by paul maggi]
def coincCheckMS(timeStamps, ene, modList, xList, yList, zList, cutOffClk, peakMin, peakMax):
    '''
    deltaT, mod1, mod2, doubs, eng1, eng2 = coincCheckMS(timeStamps, ene, modList, xList, yList, zList, cutOffClk, peakMin, peakMax)

    this function should only be passed events that are single pixel events, e.g. ene[npx==1]. For each interaction, it looks up to cutOffClk beyond that point
	for a combination of events that are both within the energy window specified by peakMin and peakMax.

    input:
        timeStamps, ene, modList, xList, yList, zList:
            these are the t, edep, mdl, x, y, and z array from the output of readInC, after only selecting one pixel events
        cutOffClk:
            the number of 10 ns clock cycles beyond a point to look for coincidence events.
				example:
					cutOffSec = 1E-6 %this is a 1 us time window
					cutOffClk = cutOffSec/(10E-9) %this gives a value of 100
    output:
        deltaT:
			time difference (in clk cycles) between the two accepted coincidence points
		mod1, mod2:
			module numbers of the first and second interactions of a coincidence pair, respectively
		doubs:
			list of doubles in the standard format (eng1, x1, y1, z1, eng2, x2, y2, z2)
		t1, t2:
			time stamp of the first and second interactions of a coincidence pair, respectively

    '''
    #  note: maxPts is the max number of coincidence Singles->doubles this code looks for.
    maxPts = 100000

    print ('  -- Length of timeStamps = {} / maxPts = {}'.format(len(timeStamps), maxPts))

    deltaT = numpy.zeros((maxPts,))
    mod1 = numpy.zeros((maxPts,))
    mod2 = numpy.zeros((maxPts,))
    yList1 = numpy.zeros((maxPts,))
    yList2 = numpy.zeros((maxPts,))
    ed1 = numpy.zeros((maxPts,))
    ed2 = numpy.zeros((maxPts,))
    t1 = numpy.zeros((maxPts,))
    t2 = numpy.zeros((maxPts,))
    xList1 = numpy.zeros((maxPts,))
    xList2 = numpy.zeros((maxPts,))
    zList1 = numpy.zeros((maxPts,))
    zList2 = numpy.zeros((maxPts,))
    jjj = 0

    #  timeDiff is the next time step minus the current time stamp [always positive]
    timeDiff = timeStamps[1:] - timeStamps[0:-1]

    #  loop through all of the time steps
    for iii in range(len(timeDiff) - 1):
        checkNum = 1
        buffTime = timeDiff[iii]

        #  kills loop if deposited energy is zero
        if ene[iii] == 0:
            break

        #  loops successive time steps (using checkNum) until outside cutOff time (cutOffClk)
        while (buffTime <= cutOffClk) & ((checkNum + iii) < len(timeStamps) - 3):

            #  saves a data point if the energy value for both time steps (iii & checkNum + iii) are within energy window
            if (ene[iii] <= peakMax) & (ene[iii] >= peakMin) & (ene[checkNum + iii] <= peakMax) & (ene[checkNum + iii] >= peakMin):

                deltaT[jjj] = buffTime
                mod1[jjj] = modList[iii] + 1
                mod2[jjj] = modList[checkNum + iii] + 1
                ed1[jjj] = ene[iii]
                ed2[jjj] = ene[iii + checkNum]
                # convert time into microseconds
                t1[jjj] = timeStamps[iii] * (10E-9) / (1E-6)
                t2[jjj] = timeStamps[iii + checkNum] * (10E-9) / (1E-6)
                yList1[jjj] = yList[iii]
                yList2[jjj] = yList[checkNum + iii]
                xList1[jjj] = xList[iii]
                xList2[jjj] = xList[checkNum + iii]
                zList1[jjj] = zList[iii]
                zList2[jjj] = zList[checkNum + iii]
                jjj += 1

            #  increment checkNum and time difference (buffTime)
            checkNum = checkNum + 1
            buffTime = buffTime + timeDiff[iii + checkNum]

        #  kill loop if max number reached
        if jjj > maxPts - 2:
            break

    #  clean up output arrays (remove zeros from the end)
    mod1 = numpy.trim_zeros(mod1) - 1
    mod2 = numpy.trim_zeros(mod2) - 1
    doubs = numpy.array((numpy.trim_zeros(ed1), numpy.trim_zeros(xList1), numpy.trim_zeros(yList1), numpy.trim_zeros(zList1), numpy.trim_zeros(ed2), numpy.trim_zeros(xList2), numpy.trim_zeros(yList2), numpy.trim_zeros(zList2))).T
    eng1 = numpy.trim_zeros(ed1)
    eng2 = numpy.trim_zeros(ed2)
    t1 = numpy.trim_zeros(t1)
    t2 = numpy.trim_zeros(t2)

    return numpy.trim_zeros(deltaT), mod1, mod2, doubs, eng1, eng2, t1, t2


#  filtering raw data to produce list of multistage coincident gammas (any energy)
#   - function pulls timing/energy settings from main program
def grouping_multistage_coincidence_events(df, timing, tWindow, eWindow):

    #  converting data from pandas to numpy arrays (not really necessary)
    npx = numpy.array(df['scatters'])
    mdl = numpy.array(df['detector'])
    edep = numpy.array(df['energy'])
    x = numpy.array(df['x'])
    y = numpy.array(df['y'])
    z = numpy.array(df['z'])
    t = numpy.array(df['time'])

    #  calculating number of 10 ns clock cycles for grouping (based on timing setting)
    number_of_clock_cycles = timing/(10E-9)
    print ('   - Using timing window of {} s, thus looking at {} 10 ns clock cycles for coincidences'.format(timing, number_of_clock_cycles))

    #  convert time into seconds (t given in units of 10 ns clock cycles), so tS = t * 10E-9
    tS = t * 10E-9
    tMin = numpy.amin(tWindow[0])
    tMax = numpy.amax(tWindow[1])

    #  filtering out single pixel events (npx == 1) within time window
    checkInds = numpy.logical_and(npx == 1, tS > tMin)
    checkInds = numpy.logical_and(checkInds, tS < tMax)
    #  print time filtered data  to screen
    numFiltered = numpy.sum(checkInds)
    perFiltered = numpy.sum(checkInds) / float(npx.size)
    print ('    - data filter: single pixel events within time window: {0} to {1} s / number of filtered events: {2} ({3:.2f}%)'.format(tMin, tMax, numFiltered, perFiltered * 100.0))

    #  setting energy window (in keV)
    engMin = eWindow[0] * 1000.0
    engMax = eWindow[1] * 1000.0
    #  print energy filtered data  to screen
    engFiltered = numpy.sum(numpy.logical_and(edep > engMin, edep < engMax))
    print ('    - number of events between {} and {} kev is {}'.format(engMin, engMax, engFiltered))

    #  grouping singles into PET coincidences
    #   - looks at single pixel events in different detectors with the given time window with each gamma having a value within the given energy window
    deltaTM, modL1, modL2, eng_position_doubles, E1, E2, T1, T2 = coincCheckMS(t[checkInds], edep[checkInds], mdl[checkInds], x[checkInds], y[checkInds], z[checkInds], number_of_clock_cycles, engMin, engMax)

    #  print counts to screen
    print ('    - results -> number of coincidence events returned: {}'.format(deltaTM.size))

    #  return two arrays
    #   - time_position_doubles: eng1, x1, y1, z1, eng2, x2, y2, z2
    #   - coincidence_data: deltaT, module1, module2, eng1, eng2, time1, time2
    return eng_position_doubles, numpy.array( (deltaTM, modL1, modL2, E1, E2, T1, T2) )



### End of Functions ###
