# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('ggplot')

# ––– 1 –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


def f(t, m, b):
    return 4 * m * t + b

# a = [0.01, 0.1, 1, 10]
a = [1]
for a in a:
    data = np.genfromtxt('rnd_a_' + str(a) + '.txt')
    t = np.linspace(1, len(data[:, 0]), len(data[:, 0]))
    rxm = data[:,0]
    rym= data[:,1]
    r2 = data[:,2]
    plt.figure()
    plt.plot(t, r2 / (a**2), 'r,')
    plt.xlabel(r'$\langle r(t)^2\rangle/a^2$')
    plt.ylabel(r't/$\tau$')
    plt.savefig("r2_a_" + str(a) + ".png")
    plt.figure()
    plt.plot(rxm / a, rym / a, 'r,', alpha=0.7)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig("r_a_" + str(a) + ".png")

    par, _ = opt.curve_fit(f, t, r2)
    if (a==1):
        D = par[0]
    plt.figure()
    plt.plot(t, f(t, par[0], par[1]), '--k',label="Ausgleichsgerade")
    plt.plot(t, r2, 'r,',label="Datenpunkte")
    plt.xlabel(r'$\langle r(t)^2\rangle$')
    plt.ylabel(r't/$\tau$')
    plt.title("D = " + str(par[0]))
    plt.legend(loc='best')
    plt.savefig("r2_fit_" + str(a) + ".png")

def P(x, y, t):
    return 1 / (4*D * np.pi * t) * np.exp(-(x**2 + y**2) / (4*D * t))
ranges = [25, 35, 75]
bins = [50,70,75]
ti = [10, 100, 1000]
data = np.genfromtxt("rnd_a_1_hist.txt")
for i in range(0, 3):
    if (a == 1):
        x = data[i, ::2]
        y = data[i, 1::2]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        hist, xedges, yedges = np.histogram2d(x, y, range=[
                                              [-ranges[i], ranges[i]], [-ranges[i], ranges[i]]], bins=bins[i], normed=True)

        elements = (len(xedges) - 1) * (len(yedges) - 1)
        xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)

        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(elements)
        dx = 0.5 * np.ones_like(zpos)
        dy = dx.copy()
        dz = hist.flatten()

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,
                 zsort='max')
        xy = np.arange(-ranges[i], ranges[i], 0.1)
        ax.plot(xy, xy, P(xy, xy, ti[i]), 'r-')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig("hist" + str(ti[i]) + ".png")


# ––– 2 –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
data = np.genfromtxt('mc_int.txt')
data_2 = np.genfromtxt('mc_int_rndwalk.txt')

plt.figure()
plt.loglog(data[:, 0], np.abs(data[:, 2]), label=u"direct sampling")
plt.loglog(data_2[:, 0], np.abs(data_2[:, 1]), label=u"Liegenlassen")
plt.loglog(data_2[:, 0], np.abs(data_2[:, 2]), "y--", label=u"Reflektion")
plt.xlim(10**1, 10**6)
plt.xlabel(r"$N$")
plt.ylabel(r"$\delta$")
plt.legend(loc='best')
plt.savefig("mc_int_pi.pdf")

plt.figure()
plt.plot(np.arange(1, 10), data[-1:, 3:-1].T, 'ro')
plt.xlim(0, 10)
plt.xlabel(r"$ab$")
plt.ylabel(r"$A$")
plt.savefig("mc_int_area_ellipse.pdf")

data = np.genfromtxt('mc_int_hist.txt')

plt.figure()
_, bins, _ = plt.hist(data, range=[2.95, 3.35], bins=40)
plt.xlabel(r"approx $\pi$")
plt.ylabel(u"Einträge pro " + '%.3f' % (bins[1] - bins[0]))
plt.savefig("mc_int_pi_hist.pdf")
