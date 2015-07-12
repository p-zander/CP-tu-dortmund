# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

# ––– 1 ––––––––––
data = np.genfromtxt("mc_1.txt")

plt.plot(data[:, 0], np.tanh(data[:, 0]), 'b-')
plt.plot(-data[:, 0], np.tanh(-data[:, 0]), 'b-')
plt.plot(data[:, 0], data[:, 1], 'r.')
plt.plot(-data[:, 0], -data[:, 1], 'r.')

plt.xlabel(r"$H$")
plt.ylabel(r"$M$")
plt.savefig("mc_1.pdf")

# ––– 2 ––––––––––
# ––– a ––––––––––
# plotting snapshots of the system
data = np.genfromtxt("mc_2_init.txt")

plt.matshow(data, vmin=0, vmax=1, cmap='Blues')
plt.savefig("mc_2_init.pdf")

data = np.genfromtxt("mc_2_snap_1.txt")

plt.matshow(data, vmin=0, vmax=1, cmap='Blues')
plt.savefig("mc_2_snap_1.pdf")

data = np.genfromtxt("mc_2_snap_3.txt")

plt.matshow(data, vmin=0, vmax=1, cmap='Blues')
plt.savefig("mc_2_snap_3.pdf")

# ––– b ––––––––––
# testing the equilibration

data_1 = np.genfromtxt("mc_2_b_100.txt")
data_2 = np.genfromtxt("mc_2_b_225.txt")
data_3 = np.genfromtxt("mc_2_b_300.txt")

plt.figure()

plt.plot(data_1[:, 0], data_1[:, 1], linewidth=1,
         label=r"$k_BT = 1,00$, random")
plt.plot(data_2[:, 0], data_2[:, 1], linewidth=1,
         label=r"$k_BT = 2,25$, random")
plt.plot(data_3[:, 0], data_3[:, 1], linewidth=1,
         label=r"$k_BT = 3,00$, random")
plt.plot(data_1[:, 0], data_1[:, 2], linewidth=1,
         label=r"$k_BT = 1,00$, const")
plt.plot(data_2[:, 0], data_2[:, 2], linewidth=1,
         label=r"$k_BT = 2,25$, const")
plt.plot(data_3[:, 0], data_3[:, 2], linewidth=1,
         label=r"$k_BT = 3,00$, const")

plt.ylim(-2, 0)
plt.xlabel(r"steps")
plt.ylabel(r"$E$")
plt.legend(ncol=2, loc='best')
plt.savefig("mc_2_b.pdf")

# ––– cde ––––––––––


def C_T(data, temp):
    return np.var(data) / temp ** 2


def U_L(data):
    return 1 - (len(data) / 3) * np.sum(data ** 4) / (np.sum(data ** 2)) ** 2


def M_T(data):
    return np.abs(np.mean(data))

# plot the energy and magnetisation for long simulation durations
data_100 = np.genfromtxt("mc_2_cde_100.txt")
data_225 = np.genfromtxt("mc_2_cde_225.txt")
data_300 = np.genfromtxt("mc_2_cde_300.txt")

plt.figure()

plt.plot(data_100[:, 0], data_100[:, 1], linewidth=1, label=r"$k_BT = 1,00$")
plt.plot(data_225[:, 0], data_225[:, 1], linewidth=1, label=r"$k_BT = 2,25$")
plt.plot(data_300[:, 0], data_300[:, 1], linewidth=1, label=r"$k_BT = 3,00$")

plt.legend(loc='best')
plt.xlabel(r"steps")
plt.ylabel(r"$E$")
plt.savefig("mc_2_c_E.pdf")

plt.figure()

plt.plot(data_100[:, 0], np.abs(data_100[:, 2]),
         linewidth=1, label=r"$k_BT = 1,00$")
plt.plot(data_225[:, 0], np.abs(data_225[:, 2]),
         linewidth=1, label=r"$k_BT = 2,25$")
plt.plot(data_300[:, 0], np.abs(data_300[:, 2]),
         linewidth=1, label=r"$k_BT = 3,00$")

plt.legend(loc='best')
plt.xlabel(r"steps")
plt.ylabel(r"$M$")
plt.savefig("mc_2_c_M.pdf")

# calculate the binder-cumulant for the first system
files = np.array([100, 150, 200, 225, 229, 250, 300, 350])
values = len(data_100[:, 0])
temps = files / 100.0

C = np.zeros(len(files))
U1 = np.zeros(len(files))
M = np.zeros(len(files))

for i in range(0, len(files)):
    data = np.genfromtxt("mc_2_cde_" + str(files[i]) + ".txt")
    C[i] = C_T(data[:, 1], temps[i])
    U1[i] = U_L(data[:, 2])
    M[i] = M_T(data[:, 2])

# plot m
plt.figure()

plt.plot(temps, M, 'ro')

plt.xlim(0.5, 4)
plt.xlabel(r"$T$")
plt.ylabel(r"$M$")
plt.savefig("M_T.pdf")

# plot c
plt.figure()

plt.plot(temps, C, 'r.')

plt.xlim(0.5, 4)
plt.xlabel(r"$T$")
plt.ylabel(r"$\frac{c}{k_B}$")
plt.savefig("C_T.pdf")

# calculate all the other binder-cumulants

U2 = np.zeros(len(files))
for i in range(0, len(files)):
    data = np.genfromtxt("mc_2_U1_" + str(files[i]) + ".txt")
    U2[i] = U_L(data[:, 1])

U3 = np.zeros(len(files))
for i in range(0, len(files)):
    data = np.genfromtxt("mc_2_U2_" + str(files[i]) + ".txt")
    U3[i] = U_L(data[:, 1])

U4 = np.zeros(len(files))
for i in range(0, len(files)):
    data = np.genfromtxt("mc_2_U3_" + str(files[i]) + ".txt")
    U4[i] = U_L(data[:, 1])

U5 = np.zeros(len(files))
for i in range(0, len(files)):
    data = np.genfromtxt("mc_2_U4_" + str(files[i]) + ".txt")
    U5[i] = U_L(data[:, 1])

U6 = np.zeros(len(files))
for i in range(0, len(files)):
    data = np.genfromtxt("mc_2_U5_" + str(files[i]) + ".txt")
    U6[i] = U_L(data[:, 1])

# plot all binder-cumulants
plt.figure()

plt.plot(temps, U1, 'b-', label=r"$N=100$")
plt.plot(temps, U2, 'g--', label=r"$N=50$")
plt.plot(temps, U3, 'r-.', label=r"$N=25$")
plt.plot(temps, U4, 'c:', label=r"$N=20$")
plt.plot(temps, U5, 'm-', label=r"$N=10$")
plt.plot(temps, U6, 'y--', label=r"$N=5$")

plt.xlim(0.5, 4)
plt.ylim(0, 1)
plt.xlabel(r"$T$")
plt.ylabel(r"$U_L$")
plt.legend(loc='best')
plt.savefig("U_T.pdf")