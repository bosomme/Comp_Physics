# -*- coding: utf-8 -*-
"""
Program that plots from .txt-file
First column is x-value
Second column is y-value

Put name of txt-file in "filename", and edit title and labels

Always reset filename = ".txt"
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

#Change this name!#
filename = "L=100.txt"
#Change this name!#

if filename == ".txt":
   sys.exit("Error: please put in filename, title and labels")

file = open(filename, 'r')
file.readline()
file.readline()

data = file.readlines()
"""MC_cycles = []; E = []; M = []

for line in data:
    MC_cycles.append(line.split()[0])
    E.append(line.split()[1])
    M.append(line.split()[2])
"""
T = np.zeros(len(data)); E = np.zeros(len(data)); M = np.zeros(len(data))
CV = np.zeros(len(data)); X = np.zeros(len(data))

for i in range(len(data)):
    T[i] = data[i].split()[0]
    E[i] = data[i].split()[1]
    M[i] = data[i].split()[2]
    CV[i] = data[i].split()[3]
    X[i] = data[i].split()[4]

file.close()
"""
E_hist = E[3000:]

plt.hist(E_hist,bins=np.linspace(min(E_hist), max(E_hist), 100))
plt.xlabel('Mean Energy')
plt.ylabel('Number of Occurances')
"""
plt.subplot(411)
plt.ylabel('Mean Energy')
plt.plot(T, E)

plt.subplot(412)
plt.ylabel('Mean Magnetization')
plt.plot(T, M)

plt.subplot(413)
plt.ylabel('Heat Capacity')
plt.plot(T, CV)

plt.subplot(414)
plt.xlabel('Temperature')
plt.ylabel('Susceptibility')
plt.plot(T, X)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #Scientific x-axis

plt.show()
