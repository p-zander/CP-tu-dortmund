# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('ggplot')

# ––– 1 –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


def P(x, y, t):
    return 1 / (0.65 * np.pi * t) * np.exp(-(x**2 + y**2) / (0.65 * t))


def f(t, m, b):
    return 4 * m * t + b

# a = [0.01, 0.1, 1, 10]
a = [1]
for a in a:
    data = np.genfromtxt('rnd_a_' + str(a) + '.txt')
    N = len(data[0, 0::2])
    t = np.linspace(1, len(data[:, 0]) + 1, len(data[:, 0]))
    rx = data[:, ::2]
    ry = data[:, 1::2]
    rxm = np.sum(rx, axis=1) / N
    rym = np.sum(ry, axis=1) / N
    r2 = np.sum((rx**2 + ry**2), axis=1) / N
    plt.figure()
    plt.plot(t, r2 / (a**2), 'r,')
    plt.savefig("r2_a_" + str(a) + ".png")
    plt.figure()
    plt.plot(rxm / a, rym / a, 'r,', alpha=0.7)
    plt.savefig("r_a_" + str(a) + ".png")

    par, _ = opt.curve_fit(f, t, r2)
    plt.figure()
    plt.plot(t, f(t, par[0], par[1]), '-k')
    plt.plot(t, r2, 'r,')
    plt.title("D = " + str(par[0]))
    plt.savefig("r2_fit_" + str(a) + ".png")

    ranges = [25, 35, 50]
    ti = [9, 99, 999]
    for i in range(0, 3):
        if (a == 1):
            x = data[ti[i], ::2]
            y = data[ti[i], 1::2]
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            hist, xedges, yedges = np.histogram2d(x, y, range=[
                                                  [-ranges[i], ranges[i]], [-ranges[i], ranges[i]]], bins=2 * ranges[i], normed=True)

            elements = (len(xedges) - 1) * (len(yedges) - 1)
            xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)

            xpos = xpos.flatten()
            ypos = ypos.flatten()
            zpos = np.zeros(elements)
            dx = 0.5 * np.ones_like(zpos)
            dy = dx.copy()
            dz = hist.flatten()

            ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='r',
                     zsort='average')
            xy = np.arange(-ranges[i], ranges[i], 0.1)
            ax.plot(xy, xy, P(xy, xy, ti[i]+1), 'r-')
            plt.savefig("hist" + str(ti[i]+1) + ".png")


# ––– 2 –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
data = np.genfromtxt('mc_int.txt')

plt.figure()
plt.loglog(data[:, 0], np.abs(data[:, 2]))
plt.xlim(10**1, 10**6)
plt.savefig("mc_int_pi.pdf")

data = np.genfromtxt('mc_int_hist.txt')

plt.figure()
plt.hist(data, bins=25)
plt.savefig("mc_int_pi_hist.pdf")

data = np.genfromtxt('mc_int_rndwalk.txt')

plt.figure()
plt.loglog(data[:, 0], np.abs(data[:, 1]))
plt.loglog(data[:, 0], np.abs(data[:, 2]))
plt.savefig("mc_int_rndwalk.pdf")
