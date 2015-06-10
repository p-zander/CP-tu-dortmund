# coding: utf-8

from matplotlib import pyplot as plt
import numpy as np

plt.style.use('ggplot')

data = np.genfromtxt("lcg.txt")

plt.figure()
plt.subplots_adjust(hspace=0.3)
plt.subplot(221)
plt.hist(data[:1e4, 0], normed=True, range=(0, 1))
plt.title('LCG 1')

plt.subplot(222)
plt.hist(data[:1e4, 1], normed=True, range=(0, 1))
plt.title('LCG 2')

plt.subplot(223)
plt.hist(data[:1e4, 2], normed=True, range=(0, 1))
plt.title('LCG 3')

plt.subplot(224)
plt.hist(data[:1e4, 3], normed=True, range=(0, 1))
plt.title('LCG 4')

plt.savefig("lcg_hist.pdf")

plt.figure()
plt.subplot(221)
plt.plot(data[:6074:2, 0], data[1:6075:2, 0], 'b,')
plt.title('LCG1')

plt.subplot(222)
plt.plot(data[:256:2, 1], data[1:256:2, 1], 'b,')
plt.title('LCG2')

plt.subplot(223)
plt.plot(data[::2, 2], data[1::2, 2], 'b,')
plt.title('LCG3')

plt.subplot(224)
plt.plot(data[::2, 3], data[1::2, 3], 'b,')
plt.title('LCG4')

plt.savefig('lcg_corr.png')
