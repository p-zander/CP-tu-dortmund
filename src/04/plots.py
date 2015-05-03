
# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize

plt.style.use('ggplot')

data = np.genfromtxt("1_ab.txt")

plt.plot(data[:, 0], 1000*data[:, 1], label=r"$(x,y) = (18,1)$cm")
plt.plot(data[:, 0], 1000*data[:, 2], label=r"$(x,y) = (16,0)$cm")
plt.plot(data[:, 0], 1000*data[:, 3], label=r"$(x,y) = (16,2)$cm")
plt.xlim(0,2*np.pi)
plt.ylabel(r"$V(\theta)$ in mJ")
plt.xlabel(r"$\theta$")
plt.legend(loc='best')
plt.savefig("1_ab_V.pdf")

plt.figure()
plt.plot(data[:, 0], 1000*data[:, 4], "r-", label=r"$F$ analytisch")
plt.plot(data[:, 0], 1000*data[:, 5], "b--", label=r"$F$ numerisch")
plt.xlim(0,2*np.pi)
plt.ylabel(r"$F(\theta)$ in mN")
plt.xlabel(r"$\theta$")
plt.legend(loc='best')
plt.savefig("1_ab_F.pdf")

# plt.figure()
# diff = data[:, 4]-data[:, 5]
# plt.plot(data[:, 0], diff)
# plt.ylim(np.min(diff), np.max(diff))
# plt.xlim(0,2*np.pi)
# plt.ylabel(u"Differenz numerische-analytische Lösung")
# plt.xlabel(r"$\theta$")
# plt.savefig("1_ab_F_diff.pdf")

del data#, diff


data1602 = np.genfromtxt("1_b_16_0_2.txt")
data1604 = np.genfromtxt("1_b_16_0_4.txt")
data1622 = np.genfromtxt("1_b_16_2_2.txt")
data1624 = np.genfromtxt("1_b_16_2_4.txt")

plt.figure()
plt.plot(data1602[:, 0], data1602[:, 1], label=r"$(x,y) = (16,0)$cm")#$\theta_0 = 2$
plt.plot(data1604[:, 0], data1604[:, 1], label=r"$(x,y) = (16,0)$cm")#$\theta_0 = 4$
plt.plot(data1622[:, 0], data1622[:, 1], label=r"$(x,y) = (16,2)$cm")#$\theta_0 = 2$
plt.plot(data1624[:, 0], data1624[:, 1], label=r"$(x,y) = (16,2)$cm")#$\theta_0 = 4$

plt.ylim(1.5,4.5)
plt.ylabel(r"$\theta(t)$")
plt.xlabel(r"Zeit $t$ in s")
plt.legend(ncol=2, loc='best')
plt.savefig("1_b.pdf")

del data1602, data1604, data1622, data1624


data = np.genfromtxt("1_c.txt")

plt.figure()
plt.plot(data[:len(data)/2, 0], data[:len(data)/2, 1], label=r"$(x,y)(t) = \left(16, 2\frac{t}{100}\right)$cm")
plt.plot(data[len(data)/2:, 0], data[len(data)/2:, 1], label=r"$(x,y)(t) = \left(16, 2-\frac{2t}{100}\right)$cm")
plt.ylabel(r"$\theta(t)$")
plt.xlabel(r"Zeit $t$ in s")
plt.legend(loc='best')
plt.savefig("1_c.pdf")

del data


data = np.genfromtxt("2_a.txt")

x = [21.800, 21.900, 21.925, 21.950, 22.00, 22.250]

for i in np.arange(2, 13, 2):
    plt.figure()
    plt.plot(data[:,i-1], data[:,i], label="Phasenraum "+`x[int(i/2-1)]`+" cm")
    plt.xlim(-2, 2)
    plt.ylim(-0.6, 0.8)
    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$\omega$ in Hz")
    plt.legend(loc='upper left')
    plt.savefig("2_a_"+`1000*x[int(i/2-1)]`+".pdf")
    
del data, x


data = np.genfromtxt("2_b.txt")

plt.figure()
plt.plot(data[:,1], data[:,2], '.', label=u"Poincaré-Schnitt an "+r"$x = 22$cm")
plt.xlabel(r"$\theta$")
plt.ylabel(r"$\omega$ in Hz")
plt.legend(loc='best')

del data


data = np.genfromtxt("2_c.txt")

plt.figure(figsize=[8,6])
plt.plot(data[:,0]*100, data[:,1], '.', alpha=0.4)
plt.xlabel(r"$x$ in cm")
plt.xlim(21.8,22)
plt.ylabel(r"$\theta$")
#plt.legend(loc='best')
plt.savefig("2_c.pdf")

del data

