#Calculate isotopic distribution of macromolecules (high approximation) using array manipulation and FFT based algorithm by Rockwood et al.(1995)

import math
import re
import numpy as np
import matplotlib.pyplot as plt

formula = raw_input("Enter the chemical formula: ") #Input formula, eg: C6H12O6
res = raw_input("Enter required resolution: ")      #Resolution


elem_limit = 5
mass_limit = 5000000

ele = re.findall('[A-Z]', formula)
elem_num = re.findall('[0-9]+', formula)

Z = []
S = []
M = np.zeros(5)

if 'H' in ele:
    M[0] = elem_num[ele.index('H')]

if 'C' in ele:
    M[1] = elem_num[ele.index('C')]

if 'N' in ele:
    M[2] = elem_num[ele.index('N')]

if 'O' in ele:
    M[3] = elem_num[ele.index('O')]

if 'S' in ele:
    M[4] = elem_num[ele.index('S')]  #M contains empirical formula of the macromolecule [H C O N S]

integer_tot_mass = M[0]*1 + M[1]*12 + M[2]*14 + M[3]*16 + M[4]*32

A = np.zeros((elem_limit, mass_limit))   #Isotopic abundance matrix

A[0][10078] = 0.9998443     #H1
A[0][20141] = 0.0001557     #H2
A[1][120000] = 0.98889      #C12
A[1][130033] = 0.01111      #C13
A[2][140030] = 0.99634      #N14
A[2][150001] = 0.00366      #N15
A[3][159949] = 0.997628     #O16
A[3][169991] = 0.000372     #O17
A[3][179991] = 0.002000     #O18
A[4][319720] = 0.95018      #S32
A[4][329714] = 0.00750      #S33
A[4][339678] = 0.04215      #S34
A[4][359670] = 0.00017      #S36

ele_fft = np.fft.fft(A)     #Element wise FFT of the isotopic abundance matrix

multi_fft = np.ones((1, mass_limit))

for i in range(elem_limit):
    dA = ele_fft[i,:]
    prdt = np.power(dA, M[i])
    multi_fft = np.multiply(multi_fft, prdt)     #Multiplying transforms

inv_fft = np.fft.ifft(multi_fft)    #Inverse FFT of the multiplied transforms to get convolutions
real_fft = np.real(inv_fft)         #Real part of inverse FFT

iso_abundance = np.zeros((1, mass_limit))

for s in range(mass_limit-1):
   iso_abundance[0][s] = real_fft[0][s+1]

mass = np.arange(1, mass_limit+1)

j = np.where(iso_abundance[0] == np.max(iso_abundance[0]))  #Relative abundance %, highest abundance to 100
iso_abundance[0][j] = 1

for d in range(len(mass)):
	if iso_abundance[0][d]*100 > float(res):
		Z.append(float(mass[d])/10000)
		S.append(100*iso_abundance[0][d])
		print(float(mass[d])/10000, 100*iso_abundance[0][d])

plt.bar(Z, S, width=0.02)
plt.xlabel('Integer Mass')
plt.ylabel('Relative Abundance')
plt.show()