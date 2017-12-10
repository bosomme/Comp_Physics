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


def analytical_solution(t, dx):
    nx = 1./dx + 1
    X = np.linspace(0, 1, nx)
    x, y = np.meshgrid(X, X)
    
    pi = np.pi
    U = np.sin(pi*x)*np.sin(pi*y)*np.exp(-2*pi**2*t)
    
    return X, U


#Change this name!#
filenames = ("2-D_t=0.000000_dx=0.100000.txt", "2-D_t=0.000000_dx=0.010000.txt", "2-D_t=0.050000_dx=0.100000.txt", "2-D_t=0.050000_dx=0.010000.txt", "2-D_t=0.500000_dx=0.100000.txt", "2-D_t=0.500000_dx=0.010000.txt")
#Change this name!#


for i in range(len(filenames)):
    plt.figure(figsize=(5,4))
    X, U = read_from_file(filenames[i])
    if i == 0:
        X1, U1 = analytical_solution(0.0, 0.1)
        UU = U - U1
    if i == 1:
        X1, U1 = analytical_solution(0.0, 0.01)
        UU = U - U1
    if i == 2:
        X1, U1 = analytical_solution(0.05, 0.1)
        UU = U - U1
    if i == 3:
        X1, U1 = analytical_solution(0.05, 0.01)
        UU = U - U1
    if i == 4:
        X1, U1 = analytical_solution(0.5, 0.1)
        UU = U - U1
    if i == 5:
        X1, U1 = analytical_solution(0.5, 0.01)
        UU = U - U1
    
    plt.imshow(UU, extent=(0,1,0,1), cmap="gist_gray") #, vmin=0, vmax=1)
    plt.colorbar()
    plt.xlabel('x'); plt.ylabel('y')
    plt.savefig('2D_error(%d).png' %i, dpi=300)
    plt.show()




"""
for t in (0, 0.05, 0.5):
    X, U = analytical_solution(t, 0.01)
    plt.imshow(U, extent=(0,1,0,1), cmap="gist_gray") #, vmin=0, vmax=1
    plt.colorbar()
    plt.show()
"""






