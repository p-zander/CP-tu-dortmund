
import matplotlib.pyplot as plt

plt.figure()
plt.plot(np.tile(np.linspace(rmin1 + dr1, rmin1 + dr1 * (1 + steps), steps),
                 steady), data1, 'b,', alpha=0.6, label=r'$x_{n+1} = rx_n(1-x_n)$')
plt.legend(loc=0)
plt.savefig("bifurk1.pdf")

plt.figure()
plt.plot(np.tile(np.linspace(rmin2 + dr2, rmin2 + dr2 * (1 + steps), steps),
                 steady), data2, 'b,', alpha=0.6, label=r'$x_{n+1} = rx_n - x_n^3$')
plt.legend(loc=0)
plt.savefig("bifurk2.pdf")

print "done!"