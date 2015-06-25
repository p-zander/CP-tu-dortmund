# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

plt.style.use('ggplot')

# ––– 1 –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


def f(d, t):
    return 4 * d * t

data = np.genfromtxt('rnd.txt')

plt.figure()
plt.plot(data[:, 0], data[:, 1], 'r,', alpha=0.2)
plt.savefig("rw.png")

timesteps = [10, 50, 100, 500, 1000]
steps = len(data)
for ts in timesteps:
    t = np.linspace(0, steps / ts, steps / ts)
    par, _ = opt.curve_fit(f, t, data[::ts, 2])
    plt.figure()
    plt.plot(t, data[::ts, 2], 'r,')
    plt.plot(t, f(par[0], t), '-k')
    plt.title("D= " + str(par[0] / ts))
    plt.savefig("r2_" + str(ts) + ".png")

timesteps = [10, 100, 1000]
for ts in timesteps:
    plt.figure()
    plt.hist2d(
        data[::ts, 0], data[::ts, 1], bins=150, normed=True, cmap="Reds")
    plt.savefig("hist_" + str(ts) + ".png")

# ––– 2 –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
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
plt.loglog(data[:,0], np.abs(data[:,1]))
plt.loglog(data[:,0], np.abs(data[:,2]))
plt.savefig("mc_int_rndwalk.pdf")
