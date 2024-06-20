# **********************************************************
#  replaces the MATLAB script Calibration.m
# __author__ = "Shanyn Hart"
# __date__ = "2022-08-23"
# **********************************************************
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
import glob
import os
from scipy.optimize import curve_fit

# **************************************************************************************************
# ****************************************** POLYNOMIALS *******************************************
# **************************************************************************************************
def first_order_polynomial(x, a, b):
    return a*x + b

def second_order_polynomial(x, a, b, c):
    return a*x**2 + b*x + c

def third_order_polynomial(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

# **************************************************************************************************
# First 11 peaks correspond to 152Eu R84
# Last peak corresponds to AmBe R86
PeakEnergyValuesFullRange = np.array([121.8, 244.7, 344.3, 411.1, 444, 511, 778.9, 867.4, 964.1, 1408, 1470])

# **************************************************************************************************
# **************************************** LaBr3 Detector 0 ****************************************
# **************************************************************************************************
PeakChannelValues0FullRange = np.array([])
PeakFWHM0 = np.array([])
PeakChannelError0FullRange = np.array([])


# --------------------------------------------------------------------------------------------------
Calibration1stOrder_FullRange_L0 = curve_fit(first_order_polynomial, PeakChannelValues0FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError0FullRange, absolute_sigma=True)
a1st_FullRange_L0 = Calibration1stOrder_FullRange_L0[0][0]
b1st_FullRange_L0 = Calibration1stOrder_FullRange_L0[0][1]
'''
plt.errorbar(PeakChannelValues0FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError0FullRange, fmt='o', label='L0')
plt.plot(PeakChannelValues0FullRange, first_order_polynomial(PeakChannelValues0FullRange, a1st_FullRange_L0, b1st_FullRange_L0), label='L0 1st order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'First Order Calibration of L0 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.6, r'$E = a_1 \cdot C + b_1$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif') 
plt.text(0.5, 0.55, r'$a_1 = $' + str(round(a1st_FullRange_L0, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.5, r'$b_1 = $' + str(round(b1st_FullRange_L0, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''


# --------------------------------------------------------------------------------------------------
Calibration2ndOrder_FullRange_L0 = curve_fit(second_order_polynomial, PeakChannelValues0FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError0FullRange, absolute_sigma=True)
a2nd_FullRange_L0 = Calibration2ndOrder_FullRange_L0[0][0]
b2nd_FullRange_L0 = Calibration2ndOrder_FullRange_L0[0][1]
c2nd_FullRange_L0 = Calibration2ndOrder_FullRange_L0[0][2]

'''
plt.errorbar(PeakChannelValues0FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError0FullRange, fmt='o', label='L0')
plt.plot(PeakChannelValues0FullRange, second_order_polynomial(PeakChannelValues0FullRange, a2nd_FullRange_L0, b2nd_FullRange_L0, c2nd_FullRange_L0), label='L0 2nd order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'Second Order Calibration of L0 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.65, r'$E = a_2 \cdot C^2 + b_2 \cdot C + c_2$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.6, r'$a_2 = $' + str(round(a2nd_FullRange_L0, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.55, r'$b_2 = $' + str(round(b2nd_FullRange_L0, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.5, r'$c_2 = $' + str(round(c2nd_FullRange_L0, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

# **************************************************************************************************
# **************************************** LaBr3 Detector 1 ****************************************
# **************************************************************************************************
PeakChannelValues1FullRange = np.array([])
PeakFWHM1 = np.array([])
PeakChannelError1FullRange = np.array([])

# --------------------------------------------------------------------------------------------------
Calibration1stOrder_FullRange_L1 = curve_fit(first_order_polynomial, PeakChannelValues1FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError1FullRange, absolute_sigma=True)
a1st_FullRange_L1 = Calibration1stOrder_FullRange_L1[0][0]
b1st_FullRange_L1 = Calibration1stOrder_FullRange_L1[0][1]

'''
plt.errorbar(PeakChannelValues1FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError1FullRange, fmt='o', label='L1')
plt.plot(PeakChannelValues1FullRange, first_order_polynomial(PeakChannelValues1FullRange, a1st_FullRange_L1, b1st_FullRange_L1), label='L1 1st order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'First Order Calibration of L1 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.8, r'$E = a_1 \cdot C + b_1$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.75, r'$a_1 = $' + str(round(a1st_FullRange_L1, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.7, r'$b_1 = $' + str(round(b1st_FullRange_L1, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

# --------------------------------------------------------------------------------------------------
Calibration2ndOrder_FullRange_L1 = curve_fit(second_order_polynomial, PeakChannelValues1FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError1FullRange, absolute_sigma=True)
a2nd_FullRange_L1 = Calibration2ndOrder_FullRange_L1[0][0]
b2nd_FullRange_L1 = Calibration2ndOrder_FullRange_L1[0][1]
c2nd_FullRange_L1 = Calibration2ndOrder_FullRange_L1[0][2]

'''
plt.errorbar(PeakChannelValues1FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError1FullRange, fmt='o', label='L1')
plt.plot(PeakChannelValues1FullRange, second_order_polynomial(PeakChannelValues1FullRange, a2nd_FullRange_L1, b2nd_FullRange_L1, c2nd_FullRange_L1), label='L1 2nd order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'Second Order Calibration of L1 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.8, r'$E = a_2 \cdot C^2 + b_2 \cdot C + c_2$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.75, r'$a_2 = $' + str(round(a2nd_FullRange_L1, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.7, r'$b_2 = $' + str(round(b2nd_FullRange_L1, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.65, r'$c_2 = $' + str(round(c2nd_FullRange_L1, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

# **************************************************************************************************
# **************************************** LaBr3 Detector 2 ****************************************
# **************************************************************************************************
PeakChannelValues2FullRange = np.array([])
PeakFWHM2 = np.array([])
PeakChannelError2FullRange = np.array([])

# --------------------------------------------------------------------------------------------------
Calibration1stOrder_FullRange_L2 = curve_fit(first_order_polynomial, PeakChannelValues2FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError2FullRange, absolute_sigma=True)
a1st_FullRange_L2 = Calibration1stOrder_FullRange_L2[0][0]
b1st_FullRange_L2 = Calibration1stOrder_FullRange_L2[0][1]

'''
plt.errorbar(PeakChannelValues2FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError2FullRange, fmt='o', label='L2')
plt.plot(PeakChannelValues2FullRange, first_order_polynomial(PeakChannelValues2FullRange, a1st_FullRange_L2, b1st_FullRange_L2), label='L2 1st order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'First Order Calibration of L2 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.8, r'$E = a_1 \cdot C + b_1$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.75, r'$a_1 = $' + str(round(a1st_FullRange_L2, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.7, r'$b_1 = $' + str(round(b1st_FullRange_L2, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

# --------------------------------------------------------------------------------------------------
Calibration2ndOrder_FullRange_L2 = curve_fit(second_order_polynomial, PeakChannelValues2FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError2FullRange, absolute_sigma=True)
a2nd_FullRange_L2 = Calibration2ndOrder_FullRange_L2[0][0]
b2nd_FullRange_L2 = Calibration2ndOrder_FullRange_L2[0][1]
c2nd_FullRange_L2 = Calibration2ndOrder_FullRange_L2[0][2]

'''
plt.errorbar(PeakChannelValues2FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError2FullRange, fmt='o', label='L2')
plt.plot(PeakChannelValues2FullRange, second_order_polynomial(PeakChannelValues2FullRange, a2nd_FullRange_L2, b2nd_FullRange_L2, c2nd_FullRange_L2), label='L2 2nd order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'Second Order Calibration of L2 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.8, r'$E = a_2 \cdot C^2 + b_2 \cdot C + c_2$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.75, r'$a_2 = $' + str(round(a2nd_FullRange_L2, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.7, r'$b_2 = $' + str(round(b2nd_FullRange_L2, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.65, r'$c_2 = $' + str(round(c2nd_FullRange_L2, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

# **************************************************************************************************
# **************************************** LaBr3 Detector 3 ****************************************
# **************************************************************************************************
PeakChannelValues3FullRange = np.array([])
PeakFWHM3 = np.array([])
PeakChannelError3FullRange = np.array([])

# --------------------------------------------------------------------------------------------------
Calibration1stOrder_FullRange_L3 = curve_fit(first_order_polynomial, PeakChannelValues3FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError3FullRange, absolute_sigma=True)
a1st_FullRange_L3 = Calibration1stOrder_FullRange_L3[0][0]
b1st_FullRange_L3 = Calibration1stOrder_FullRange_L3[0][1]

'''
plt.errorbar(PeakChannelValues3FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError3FullRange, fmt='o', label='L3')
plt.plot(PeakChannelValues3FullRange, first_order_polynomial(PeakChannelValues3FullRange, a1st_FullRange_L3, b1st_FullRange_L3), label='L3 1st order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'First Order Calibration of L3 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.8, r'$E = a_1 \cdot C + b_1$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.75, r'$a_1 = $' + str(round(a1st_FullRange_L3, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.7, r'$b_1 = $' + str(round(b1st_FullRange_L3, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

# --------------------------------------------------------------------------------------------------
Calibration2ndOrder_FullRange_L3 = curve_fit(second_order_polynomial, PeakChannelValues3FullRange, PeakEnergyValuesFullRange, sigma=PeakChannelError3FullRange, absolute_sigma=True)
a2nd_FullRange_L3 = Calibration2ndOrder_FullRange_L3[0][0]
b2nd_FullRange_L3 = Calibration2ndOrder_FullRange_L3[0][1]
c2nd_FullRange_L3 = Calibration2ndOrder_FullRange_L3[0][2]

'''
plt.errorbar(PeakChannelValues3FullRange, PeakEnergyValuesFullRange, yerr=PeakChannelError3FullRange, fmt='o', label='L3')
plt.plot(PeakChannelValues3FullRange, second_order_polynomial(PeakChannelValues3FullRange, a2nd_FullRange_L3, b2nd_FullRange_L3, c2nd_FullRange_L3), label='L3 2nd order fit')
plt.xlabel('Channel Number', font = 'serif', fontsize=18)
plt.ylabel('Energy (keV)', font = 'serif', fontsize=18)
plt.title(r'Second Order Calibration of L3 - Peaks Obtained from 152Eu and AmBe', fontsize=22, fontweight='bold', font = 'serif')
plt.text(0.5, 0.8, r'$E = a_2 \cdot C^2 + b_2 \cdot C + c_2$', fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.75, r'$a_2 = $' + str(round(a2nd_FullRange_L3, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.7, r'$b_2 = $' + str(round(b2nd_FullRange_L3, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.text(0.5, 0.65, r'$c_2 = $' + str(round(c2nd_FullRange_L3, 4)), fontsize=18, transform=plt.gcf().transFigure, font = 'serif')
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
plt.legend(prop={'family':'serif', 'size':18})
plt.show()
'''

