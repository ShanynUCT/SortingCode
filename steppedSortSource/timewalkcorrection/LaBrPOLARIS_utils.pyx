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
  
    detector_numbers = set(df.detector.values)
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
        #xyz = detector_slices[i][:, 4:6]

        #  the POLARIS detectors have left handed Coordinate system (so we first flip y-axis)
        xyz[:,1] *= -1  # flip y-axis
        print ('   - detector {} - converted to right handed coordinates (flip y-axis)'.format(i))
        xyz = xyz.T
        xyz_four_length = numpy.append(xyz, [numpy.ones(xyz.shape[1])], axis=0)

        #  the detector to the reconstruction transformation is defined as the inverse transformation
        M = numpy.linalg.inv(transformation_matrices[n])
        xyz_transformed = numpy.dot(M, xyz_four_length)
        detector_slices[i][:, 3:6] = xyz_transformed[:-1, :].T
        #detector_slices[i][:, 4:6] = xyz_transformed[:-1, :].T
        print ('   - detector {} - successful coordinate transformation!'.format(i))

    #  recombines transformed coordinates back into full dataFrame
    new_df = pandas.DataFrame(numpy.concatenate(detector_slices), columns=df.columns)

    #  re-sort data by time stamp [important for coincidence grouping]
    new_df = new_df.sort_values(by = ['time', 'y'], ascending = [True, True])
    #  not sure if this is important
    new_df = new_df.reset_index(drop = True)

    return new_df

