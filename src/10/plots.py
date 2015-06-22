import matplotlib.pyplot as plt
import numpy as np

plt.style.use('ggplot')

data = np.genfromtxt('rnd.txt')

plt.figure()
plt.plot(data[:, 0], data[:, 1], 'r,', alpha=0.2)
plt.savefig("rw.png")

plt.figure()
plt.plot(np.linspace(0,1000,1000),data[::100, 2], 'r,')
plt.savefig("r2.png")
