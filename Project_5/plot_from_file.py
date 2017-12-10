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

    X = np.zeros(len(data)); U = np.zeros(len(data));

    for i in range(len(data)):
        X[i] = data[i].split()[0]
        U[i] = data[i].split()[1]
    
    file.close()
    return X, U


def analytical_solution(t, dx):
    nx = 1./dx + 1
    Xx = np.linspace(0, 1, nx)
    Uu = np.zeros(len(Xx))
    pi = np.pi
    for x in range(len(Xx)):
        for k in range(1, len(Xx)):
            k_pi = k*pi
            Uu[x] += (2*(-1)**k)/k_pi * np.exp(-(k_pi**2)*t) * np.sin(Xx[x]*k_pi)
        Uu[x] += Xx[x]
    return Xx, Uu


#Change this name!#
# filenames = ("F_E_t=0.050000_dx=0.100000.txt", "F_E_t=0.500000_dx=0.100000.txt", "B_E_t=0.050000_dx=0.100000.txt", "B_E_t=0.500000_dx=0.100000.txt", "C-N_t=0.050000_dx=0.100000.txt", "C-N_t=0.500000_dx=0.100000.txt")
filenames = ("F_E_t=0.050000_dx=0.100000.txt", "F_E_t=0.500000_dx=0.100000.txt", "B_E_t=0.050000_dx=0.100000.txt", "B_E_t=0.500000_dx=0.100000.txt", "C-N_t=0.050000_dx=0.100000.txt", "C-N_t=0.500000_dx=0.100000.txt")
#Change this name!#

legends = ('Forward Euler t=0.05','Forward Euler t=0.5','Backward Euler t=0.05','Backward Euler t=0.5','Crank-Nicolson t=0.05','Crank-Nicolson t=0.5', 'Analytical t=0.05', 'Analytical t=0.5')

"""
if filename == ".txt":
   sys.exit("Error: please put in filename, title and labels")
"""

"""
plt.figure(figsize=(7,3.6))
for i in range(len(filenames)):
    X, U = read_from_file(filenames[i])
    plt.plot(X, U, label=legends[i])

for t in (0.05, 0.5):
    X, U = analytical_solution(t)
    plt.plot(X, U)

plt.legend(legends)

plt.xlabel('x')
plt.ylabel('U(x,t)')

plt.savefig('dx=0.01.png', dpi=300) #saving figure
"""


# plot error
plt.figure(figsize=(7,3.6))
legends = ('Forward Euler', 'Backward Euler', 'Crank-Nicolson')

X, U = read_from_file(filenames[0])
analX, analU = analytical_solution(0.05,0.1)
plt.plot(X, abs(U-analU), label=legends[0])


X, U = read_from_file(filenames[2])
analX, analU = analytical_solution(0.05,0.1)
plt.plot(X, abs(U-analU), label=legends[1])

X, U = read_from_file(filenames[4])
analX, analU = analytical_solution(0.05,0.1)
plt.plot(X, abs(U-analU), label=legends[2])

plt.legend(legends)

plt.xlabel('x')
plt.ylabel('U(x,t)')

plt.savefig('Error_dx=0.1.png', dpi=300) #saving figure


plt.show()



