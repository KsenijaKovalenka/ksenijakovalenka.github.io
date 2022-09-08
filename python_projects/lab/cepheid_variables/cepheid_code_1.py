# -*- coding: utf-8 -*-
"""
Cepheid Variables experiment.

Reads in the data file with counts, area in pixels over which counts were taken,
and times of observation (start exposure time). PLots the graph of counts
against time and makes a phase folding over a given set of peiods. String length
is calculated as a measure of fit, but the optimal period has to be picked
manually because of the quality of data. To make that easier, a plot of string
length against period is plotted. Best period minimises string length.

Final plot of the phasefolded counts against time over two periods is created
for the value of the best period (has to be entered manually) as well. 

Created on Fri Feb 12 15:40:19 2021

author: Ksenija Kovalenka
"""

import numpy as np
import matplotlib.pyplot as plt

#  read in the data
data_set = np.genfromtxt('data1.txt', delimiter=',', comments='%')
times = data_set[:, 4] - data_set[0, 4]

data_set_other = np.genfromtxt('data14.txt', delimiter=',', comments='%')
counts_1 = data_set_other[:, 0]
counts_2 = data_set_other[:, 1]
pixels_1 = data_set_other[:, 2]
pixels_2 = data_set_other[:, 3]

#  subtract backround count
counts = ((counts_1 / pixels_1) - ((counts_2 - counts_1) / (pixels_2 - pixels_1))) * np.mean(pixels_1)
counts_unc = np.sqrt(((1 + pixels_1 / (pixels_2 - pixels_1)) 
                      * np.sqrt(counts_1))**2 + ((pixels_1 / (pixels_2 - pixels_1) * np.sqrt(counts_1))**2))

#  calculate apparent magnitude (weighted mean)
magnitudes = 22.57 - 2.5 * np.log10(counts / 1000)
magnitudes_error = (2.5 * counts_unc) / (counts * np.log(10))

w_i = 1 / magnitudes_error ** 2
w_mean = np.sum(magnitudes * w_i) / np.sum(w_i)
w_unc = 1/ np.sqrt(np.sum(w_i))

print('m = {}'.format(w_mean))
print('m uncertainty = {}'.format(w_unc))

#  first plot
plt.scatter(times, counts)
plt.errorbar(times, counts, counts_unc, linestyle=' ')
plt.title('Cepheid 13')
plt.xlabel('time (days)')
plt.ylabel('counts')

plt.savefig('cepheid13-1.png', dpi=300, bbox_inches='tight')
plt.show()


#  finding optimal period
trial_periods = np.arange(21, 25, 0.1)
strings = np.array([])

for period in trial_periods:
    
    #  convert to phases
    phases = times / period
    
    for index in range(len(phases)):
        if phases[index] > 1:
            phases[index] -= 1
    
    #  crate new data set with phases (sort from smallest to greathest)       
    combined_data = np.vstack((counts, phases, counts_unc))
    combined_data = combined_data[:, combined_data[1].argsort()]
    
    #  calculate string
    total_string = 0
    
    for index in range(len(combined_data[1]) - 1):
        string = ((combined_data[1][index + 1] - combined_data[1][index])**2 +
                  (combined_data[0][index + 1] - combined_data[0][index])**2)
        total_string += string
    
    strings = np.append(strings, total_string)
    
#  print what computer thinks is a minimum
print('min string: {}'.format(np.min(strings)))
print('optimal period: {}'.format(trial_periods[np.where(strings == np.min(strings))]))

#  supporting plot to make a better estimate
plt.scatter(trial_periods, strings)
plt.show()

# the nicest graph x2 periods (same rountine but for only one of the periods)

period = 24

phases = times / period

for index in range(len(phases)):
    if phases[index] > 1:
        phases[index] -= 1

combined_data = np.vstack((counts, phases, counts_unc))
combined_data = combined_data[:, combined_data[1].argsort()]

counts = np.append(combined_data[0], combined_data[0])
times = np.append(combined_data[1], combined_data[1] + 1)
uncertainties = np.append(combined_data[2], combined_data[2])

plt.scatter(times, counts)
plt.errorbar(times, counts, uncertainties, linestyle=' ')
plt.title('Cepheid 13, period {}'.format(period))
plt.xlabel('time (days)')
plt.ylabel('counts')
plt.savefig('Cepheid13', dpi=300, bbox_inches='tight')
plt.show()
