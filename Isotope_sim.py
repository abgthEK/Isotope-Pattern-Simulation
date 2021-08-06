
import math
import re
import numpy as np
import matplotlib.pyplot as plt

formula = raw_input("Enter the chemical formula: ")

max_ele = 5
max_mass = 2^10

ele = re.findall('[A-Z]', formula)
elem_num = re.findall('[0-9]+', formula)

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
    M[4] = elem_num[ele.index('S')]

tot_mass = M[0]*1 + M[1]*12 + M[2]*14 + M[3]*16 + M[4]*32

A = np.zeros((max_ele, max_mass))

if M[0] > 0:
    A[0][1] = 0.998443
    A[0][2] = 0.0001557

if M[1] > 0:
    A[1][12] = 0.98889
    A[1][13] = 0.01111

if M[2] > 0:
    A[2][14] = 0.99634
    A[2][15] = 0.00366

if M[3] > 0:
    A[3][16] = 0.997628
    A[3][17] = 0.000372
    A[3][18] = 0.002000

if M[4] > 0:
    A[4][32] = 0.95018
    A[4][33] = 0.00750  
    A[4][34] = 0.04215
    A[4][35] = 0
    A[4][36] = 0.00017

tA = np.fft.fft(A)

ptA = np.ones((1, max_mass))

for i in range(max_ele):
    dA = tA[i,:]
    prd = np.power(dA, M[i])
    ptA = np.multiply(ptA, prd)

sTA = np.fft.ifft(ptA)
riptA = np.real(sTA)

id = np.zeros((1, max_mass))

for s in range(max_mass-1):
   id[0][s] = riptA[0][s+1]

l = np.arange(1, max_mass+1)

j = np.where(id[0] == np.max(id[0]))
id[0][j] = 1

plt.bar(l, 100*id[0], width=0.1)
plt.xlim([tot_mass-2, tot_mass+10])
plt.ylim([0, 105])
plt.xlabel('Integer Mass')
plt.ylabel('Relative Abundance')
plt.xticks(np.arange(tot_mass-2, tot_mass+10, 1))
plt.yticks(np.arange(0, 105, 10))
plt.show()
