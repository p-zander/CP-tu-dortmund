# coding: utf-8

from matplotlib import pyplot as plt
import numpy as np

plt.style.use('ggplot')
data = np.genfromtxt("rand.txt")


def normal(x):
    return 1.0 / (2 * np.pi)**0.5 * np.exp(-x**2 / 2)

fig = plt.figure()
plt.subplots_adjust(hspace=0.3)

fig.add_subplot(221, title='Box Muller')
plt.hist(data[:, 0], bins=50, normed=True)
plt.plot(np.linspace(-4, 4, 100), normal(np.linspace(-4, 4, 100)))

ax = fig.add_subplot(222, title='Zentraler Grenzwertsatz')
plt.hist(data[:, 1], bins=50, normed=True)
plt.plot(np.linspace(-4, 4, 100), normal(np.linspace(-4, 4, 100)))

ax = fig.add_subplot(223, title='Rueckweisungsmethode')
plt.hist(data[:, 2], bins=50, normed=True)
plt.plot(np.linspace(0, np.pi, 100), 0.5 * np.sin(np.linspace(0, np.pi, 100)))

ax = fig.add_subplot(224, title='Inversionsmethode')
plt.hist(data[:, 3], bins=50, normed=True)
plt.plot(np.linspace(-4, 4, 100), 3 * np.linspace(-4, 4, 100)**2)
plt.xlim(0, 1)
plt.ylim(0, 3)

plt.savefig("dists.png", dpi=400)
