# coding: utf-8

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.style.use('ggplot')
plt.switch_backend('AGG')

data = np.genfromtxt("lcg.txt")

fig = plt.figure()
plt.subplots_adjust(hspace=0.3)

ax = fig.add_subplot(221, title='LCG1')
plt.hist(data[:, 0], normed=True, range=(0, 1))

ax = fig.add_subplot(222, title='LCG2')
plt.hist(data[:, 1], normed=True, range=(0, 1))

ax = fig.add_subplot(223, title='LCG3')
plt.hist(data[:1e4, 2], normed=True, range=(0, 1))

ax = fig.add_subplot(224, title='LCG4')
plt.hist(data[:, 3], normed=True, range=(0, 1))

plt.savefig('lcg_hist.png', dpi=200)

fig = plt.figure()
ax = fig.add_subplot(221, title='LCG1')
ax.plot(data[::2, 0], data[1::2, 0], 'r,')

ax = fig.add_subplot(222, title='LCG2')
ax.plot(data[::2, 1], data[1::2, 1], 'r,')

ax = fig.add_subplot(223, title='LCG3')
ax.plot(data[::2, 2], data[1::2, 2], 'r,')

ax = fig.add_subplot(224, title='LCG4')
ax.plot(data[::2, 3], data[1::2, 3], 'r,')

plt.savefig('lcg_corr.png', dpi=100)


fig = plt.figure()

x = data[:-1:3, :].T
y = data[1::3, :].T
z = data[2::3, :].T

ax = fig.add_subplot(221, projection='3d', title='LCG1')
ax.scatter(y[0], z[0], x[0], c='r', marker='.', alpha=0.1)

ax = fig.add_subplot(222, projection='3d', title='LCG2')
ax.scatter(x[1], y[1], z[1], c='r', marker='.', alpha=0.1)

ax = fig.add_subplot(223, projection='3d', title='LCG3')
ax.scatter(y[2], z[2], x[2], c='r', marker='.', alpha=0.1)

ax = fig.add_subplot(224, projection='3d', title='LCG4')
ax.scatter(y[3], z[3], x[3], c='r', marker='.', alpha=0.1)

fig.savefig('lcg_corr_3.png', dpi=900)
