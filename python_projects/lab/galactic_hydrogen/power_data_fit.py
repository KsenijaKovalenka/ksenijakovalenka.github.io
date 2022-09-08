# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 19:10:11 2020

Gaussian fiiting algorithm.

This code fits a gaussian curve to the input data, and plots the data together
with the fit line.

It was used in the galactic hydrogen experiment, inputting a data of the scan
over strong sources from the 7m radio telescope.

by Ksenija Kovalenka
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit


# gaussian curve with x as independent variable and others
# as parameters of fit
# fit parameters = vertical_offset, normalisation_const, mean, std

def gauss(x, vertical_offset, normalisation_const, mean, std):
    return vertical_offset + normalisation_const * np.exp(-(x - mean) ** 2
                                                          / (2 * std ** 2))


# an initial quess of fitting parameters
# in the save order as defined in function of parameters
INITIAL_GUESS = np.array([1, 1, 1, 1])
TEMPERATURE_DATA = SCAN_OFFSET_DATA = np.array([])

DATA_FILE = open('data_taurus_A.txt', 'r')

for line in DATA_FILE:
    split_data = line.split(',')

    SCAN_OFFSET_DATA = np.append(SCAN_OFFSET_DATA, float(split_data[0]))
    TEMPERATURE_DATA = np.append(TEMPERATURE_DATA, float(split_data[1]))

DATA_FILE.close()

PARAMETERS = curve_fit(gauss, SCAN_OFFSET_DATA, TEMPERATURE_DATA, INITIAL_GUESS)

CHI_SQUARED_REDUCED = stats.chisquare(TEMPERATURE_DATA, 
                                      gauss(SCAN_OFFSET_DATA, 
                                            PARAMETERS[0][0],
                                            PARAMETERS[0][1], 
                                            PARAMETERS[0][2],
                                            PARAMETERS[0][3]), 
                                        (len(TEMPERATURE_DATA) - 4))[0] / (len(TEMPERATURE_DATA) - 4)


FWHM = PARAMETERS[0][3] * np.sqrt(8* np.log(2))
PARAMETER_ERRORS = np.sqrt(np.diag(PARAMETERS[1]))
PERCENTAGE_STD_ERROR = PARAMETER_ERRORS[-1] / PARAMETERS[0][3] * 100

# max power  = normalisation constant + vertical offset
MAX_POWER = PARAMETERS[0][1] + PARAMETERS[0][0]
# absolute error in max power
MAX_POWER_ERROR = np.sqrt(PARAMETER_ERRORS[0]**2 + PARAMETER_ERRORS[1]**2)


# visualising data and fit
plt.title('Taurus-A power scan')
plt.xlabel('Scan Offset/ Degrees')
plt.ylabel('Temperature/ K')
plt.plot(SCAN_OFFSET_DATA, TEMPERATURE_DATA)
plt.plot(SCAN_OFFSET_DATA, gauss(SCAN_OFFSET_DATA, PARAMETERS[0][0],
                                 PARAMETERS[0][1], PARAMETERS[0][2],
                                 PARAMETERS[0][3]))

# checking for goodnes of fit. Red - bad, green - good.
if 1 - np.sqrt(8 / (len(TEMPERATURE_DATA) - 4)) < CHI_SQUARED_REDUCED < 1 + np.sqrt((8 / len(TEMPERATURE_DATA))):
    plt.annotate('reduced $\chi^2$ = {0:.3f}'.format(CHI_SQUARED_REDUCED),
                 xy=(1.9, 15), c='g')
else:
    plt.annotate('reduced $\chi^2$ = {0:.3f}'.format(CHI_SQUARED_REDUCED),
                 xy=(1.9, 15), c='r')

plt.savefig('Taurus_A.png', dpi=300, bbox_inches='tight')
plt.show()

print('FWHM = {0:.4f} +- {1:.2f} %'.format(FWHM, PERCENTAGE_STD_ERROR))
print('Maximum power = {0:2.2f} +- {1:1.2f}'.format(MAX_POWER, MAX_POWER_ERROR))