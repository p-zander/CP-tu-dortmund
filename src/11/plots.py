import numpy as np
from matplotlib import pyplot as plt

data = np.genfromtxt("mc_2_init.txt")

plt.matshow(data, vmin=0, vmax=1, cmap='Blues')
plt.savefig("mc_2_init.pdf")

data = np.genfromtxt("mc_2_snap_1.txt")

plt.matshow(data, vmin=0, vmax=1, cmap='Blues')
plt.savefig("mc_2_snap_1.pdf")

data = np.genfromtxt("mc_2_snap_3.txt")

plt.matshow(data, vmin=0, vmax=1, cmap='Blues')
plt.savefig("mc_2_snap_3.pdf")

# testing the equilibration

data_1 = np.genfromtxt("mc_2_b_1_00.txt")
data_2 = np.genfromtxt("mc_2_b_2_25.txt")
data_3 = np.genfromtxt("mc_2_b_3_00.txt")

plt.figure()

plt.plot(data_1[:, 0], data_1[:, 1], linewidth=1, label=r"$k_B T = 1,00$, random")
plt.plot(data_2[:, 0], data_2[:, 1], linewidth=1, label=r"$k_B T = 2,25$, random")
plt.plot(data_3[:, 0], data_3[:, 1], linewidth=1, label=r"$k_B T = 3,00$, random")
plt.plot(data_1[:, 0], data_1[:, 2], linewidth=1, label=r"$k_B T = 1,00$, const")
plt.plot(data_2[:, 0], data_2[:, 2], linewidth=1, label=r"$k_B T = 2,25$, const")
plt.plot(data_3[:, 0], data_3[:, 2], linewidth=1, label=r"$k_B T = 3,00$, const")

plt.ylim(-2, 0)
plt.legend(ncol=2, loc='best')
plt.savefig("mc_2_b.pdf")

# samples of energy and magnetisation

data_1 = np.genfromtxt("mc_2_c_1_00.txt")
data_2 = np.genfromtxt("mc_2_c_2_25.txt")
data_3 = np.genfromtxt("mc_2_c_3_00.txt")

plt.figure()

plt.plot(data_1[:, 0], data_1[:, 1], linewidth=1, label=r"$k_B T = 1,00$")
plt.plot(data_2[:, 0], data_2[:, 1], linewidth=1, label=r"$k_B T = 2,25$")
plt.plot(data_3[:, 0], data_3[:, 1], linewidth=1, label=r"$k_B T = 3,00$")

plt.legend(loc='best')
plt.savefig("mc_2_c_E.pdf")

plt.figure()

plt.plot(data_1[:, 0], data_1[:, 2], linewidth=1, label=r"$k_B T = 1,00$")
plt.plot(data_2[:, 0], data_2[:, 2], linewidth=1, label=r"$k_B T = 2,25$")
plt.plot(data_3[:, 0], data_3[:, 2], linewidth=1, label=r"$k_B T = 3,00$")

plt.legend(loc='best')
plt.savefig("mc_2_c_M.pdf")

# magnetisation, energy and binder cumulant depending on temperature

data = np.genfromtxt("mc_2_d.txt")

plt.figure()

plt.plot(data[:, 0], data[:, 1], linewidth=1, label=r"$|M(T)|$")
plt.plot(data[:, 0], data[:, 3], linewidth=1, label=r"$U(T)$")

plt.legend(loc='best')
plt.savefig("mc_2_d_MU.pdf")

plt.figure()

plt.plot(data[:, 0], data[:, 2], linewidth=1, label=r"$c(T)$")

plt.legend(loc='best')
plt.savefig("mc_2_d_c.pdf")
