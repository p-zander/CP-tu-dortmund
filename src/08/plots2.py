# coding: utf-8

from matplotlib import pyplot as plt
import numpy as np

plt.style.use('ggplot')
data = np.genfromtxt("rand.txt")


def normal(x):
    return 1.0 / (2 * np.pi)**0.5 * np.exp(-x**2 / 2)

plt.figure()
plt.subplots_adjust(hspace=0.3)
plt.subplot(221)
plt.hist(data[:, 0], bins=100, normed=True)
plt.plot(np.linspace(-4, 4, 100), normal(np.linspace(-4, 4, 100)), 'r-')
plt.title('Box Mueller')

plt.subplot(222)
plt.hist(data[:, 1], bins=100, normed=True)
plt.plot(np.linspace(-4, 4, 100), normal(np.linspace(-4, 4, 100)), 'r-')
plt.title('Zentraler Grenzwertsatz')

plt.subplot(223)
plt.hist(data[:, 2], bins=100, normed=True)
plt.plot(np.linspace(0, np.pi, 100), 0.5 * np.sin(np.linspace(0, np.pi, 100)))
plt.title('Rueckweisungsmethode')

plt.subplot(224)
plt.hist(data[:, 3], bins=100, normed=True)
plt.plot(np.linspace(-4, 4, 100), 3 * np.linspace(-4, 4, 100)**2, 'r-')
plt.xlim(0, 1)
plt.ylim(0, 3)
plt.title('Inversionsmethode')

plt.savefig("dists.pdf")
