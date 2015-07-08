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

plt.plot(data_1[:, 0], data_1[:, 1], linewidth=1,
         label=r"$k_B T = 1,00$, random")
plt.plot(data_2[:, 0], data_2[:, 1], linewidth=1,
         label=r"$k_B T = 2,25$, random")
plt.plot(data_3[:, 0], data_3[:, 1], linewidth=1,
         label=r"$k_B T = 3,00$, random")
plt.plot(data_1[:, 0], data_1[:, 2], linewidth=1,
         label=r"$k_B T = 1,00$, const")
plt.plot(data_2[:, 0], data_2[:, 2], linewidth=1,
         label=r"$k_B T = 2,25$, const")
plt.plot(data_3[:, 0], data_3[:, 2], linewidth=1,
         label=r"$k_B T = 3,00$, const")

plt.ylim(-2, 0)
plt.legend(ncol=2, loc='best')
plt.savefig("mc_2_b.pdf")

# samples of energy and magnetisation

data_100 = np.genfromtxt("mc_2_cde_100.txt")
data_225 = np.genfromtxt("mc_2_cde_225.txt")
data_300 = np.genfromtxt("mc_2_cde_300.txt")

plt.figure()

plt.plot(data_100[:, 0], data_100[:, 1], linewidth=1, label=r"$k_B T = 1,00$")
plt.plot(data_225[:, 0], data_225[:, 1], linewidth=1, label=r"$k_B T = 2,25$")
plt.plot(data_300[:, 0], data_300[:, 1], linewidth=1, label=r"$k_B T = 3,00$")

plt.legend(loc='best')
plt.savefig("mc_2_c_E.pdf")

plt.figure()

plt.plot(data_100[:, 0], np.abs(data_100[:, 2]),
         linewidth=1, label=r"$k_B T = 1,00$")
plt.plot(data_225[:, 0], np.abs(data_225[:, 2]),
         linewidth=1, label=r"$k_B T = 2,25$")
plt.plot(data_300[:, 0], np.abs(data_300[:, 2]),
         linewidth=1, label=r"$k_B T = 3,00$")

plt.legend(loc='best')
plt.savefig("mc_2_c_M.pdf")


files = np.array([100, 150, 200, 225, 250, 300, 350])
values = len(data_100[:, 0])
temps = files / 100.0

C_T = np.zeros(7)
U_T1 = np.zeros(7)
M_T = np.zeros(7)
for i in range(0, 7):
    data = np.genfromtxt("mc_2_cde_" + str(files[i]) + ".txt")
    C_T[i] = np.std(data[:, 1])**2 / temps[i]**2
    U_T1[i] = 1 - (np.sum(data[:, 2]**4) / values) / \
        (np.sum(data[:, 2]**2) / values)**2
    M_T[i] = np.mean(data[:,2])

plt.figure()
plt.plot(temps, M_T, 'ro', label=r"M(T)")
plt.legend(loc='best')
plt.savefig("M_T.pdf")

plt.figure()
plt.plot(temps, C_T, 'r.', label=r"C(T)")
plt.legend(loc='best')
plt.savefig("C_T.pdf")

U_T2 = np.zeros(7)
for i in range(0, 7):
    data = np.genfromtxt("mc_2_U1_" + str(files[i]) + ".txt")
    U_T2[i] = 1 - (np.sum(data[:, 1]**4) / values) / \
        (np.sum(data[:, 1]**2) / values)**2

U_T3 = np.zeros(7)
for i in range(0, 7):
    data = np.genfromtxt("mc_2_U2_" + str(files[i]) + ".txt")
    U_T3[i] = 1 - (np.sum(data[:, 1]**4) / values) / \
        (np.sum(data[:, 1]**2) / values)**2

plt.figure()
plt.plot(temps, U_T1, 'ro', label=r"U(T), N=100")
plt.plot(temps, U_T2, 'kx', label=r"U(T), N=50")
plt.plot(temps, U_T3, 'b.', label=r"U(T), N=25")
plt.xlim(0.5, 4)
plt.ylim(-2.5, 0.5)
plt.legend(loc='best')
plt.savefig("U_T.pdf")

# magnetisation, energy and binder cumulant depending on temperature

# data = np.genfromtxt("mc_2_d.txt")
# data = np.sort(data, axis=0)

# plt.figure()

# plt.plot(data[:, 0], data[:, 1], linewidth=1, label=r"$|M(T)|$")
# plt.plot(data[:, 0], data[:, 3], linewidth=1, label=r"$U(T)$")

# plt.legend(loc='best')
# plt.savefig("mc_2_d_MU.pdf")

# plt.figure()

# plt.plot(data[:, 0], data[:, 2], linewidth=1, label=r"$c(T)$")

# plt.legend(loc='best')
# plt.savefig("mc_2_d_c.pdf")
