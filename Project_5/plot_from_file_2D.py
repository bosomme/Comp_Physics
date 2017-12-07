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

def read_from_file(filename):
    file = open(filename, 'r')
    #file.readline()
    #file.readline()

    data = file.readlines()

    X = np.zeros(len(data))
    U = np.zeros((len(data), len(data)))

    for i in range(len(data)):
        X[i] = data[i].split()[0]
        U[i,:] = data[i].split()[1:] #U[x, y]
    
    file.close()
    return X, U


def analytical_solution(t):
    nx = 100
    X = np.linspace(0, 1, nx)
    x, y = np.meshgrid(X, X)
    
    pi = np.pi
    U = np.sin(pi*x)*np.sin(pi*y)*np.exp(-2*pi**2*t)
    
    return X, U


#Change this name!#
filenames = ("2-D_t=0.050000_dx=0.100000.txt", "2-D_t=0.050000_dx=0.010000.txt", "2-D_t=0.500000_dx=0.100000.txt", "2-D_t=0.500000_dx=0.010000.txt")
#Change this name!#

#legends = ('Forward Euler t=0.05','Forward Euler t=0.5','Backward Euler t=0.05','Backward Euler t=0.5','Crank-Nicolson t=0.05','Crank-Nicolson t=0.5', 'Analytical t=0.05', 'Analytical t=0.5')

"""
plt.subplot(111)
for i in range(len(filenames)):
    X, U = read_from_file(filenames[i])
    plt.imshow(U, extent=(0,1,0,1), cmap="gist_gray", vmin=0, vmax=1)
    plt.colorbar()
    plt.show()
"""




for t in (0, 0.05, 0.5):
    X, U = analytical_solution(t)
    plt.imshow(U, extent=(0,1,0,1), cmap="gist_gray", vmin=0, vmax=1)
    plt.colorbar()
    plt.show()







