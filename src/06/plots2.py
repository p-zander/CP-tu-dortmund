import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
line, = plt.plot(np.linspace(-10, 10, 201), psi2[:, 0], 'b-')


def animate(i):
	line.set_xdata(np.linspace(-10, 10, 201))
	line.set_ydata(psi2[:, i])
	return line,


def init():
	line.set_ydata(np.ma.array(psi2[:, 0], mask=True))
	plt.xlim(-10, 10)
	plt.ylim(0,.04)
	return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0,499,1), interval=25, init_func=init)
ani.save("Welle.gif", writer='imagemagick', fps=30)

norm = np.sum(psi2, axis=0)
mean = np.sum( (np.linspace(-10, 10, N)*psi2.T).T , axis=0)
var  = np.sum( (np.linspace(-10, 10, N)**2 * psi2.T).T , axis=0) - mean

plt.figure()
plt.plot(np.linspace(-10, 10, N), psi2[:, 0], 'b-', label=r'$\vert\psi_n(t=0)\vert^2$')
plt.plot(np.linspace(-10, 10, N), psi2[:, t-1], 'r-', label=r'$\vert\psi_n(t=10)\vert^2$')
plt.legend()
plt.savefig("t0_10.pdf")

plt.figure()
plt.plot(np.linspace(0, t/0.02, t), norm, label=r'$\Sigma_n(\Delta\xi)\vert\psi_n(t)\vert^2$')
plt.ylim(0.9,1.1)
plt.legend()
plt.savefig("norm.pdf")

plt.figure()
plt.plot(np.linspace(0, t/0.02, t), var, label=r'$\langle\xi^2\rangle(t) - \langle\xi\rangle^2(t)$')
# plt.ylim(0.9,1.1)
plt.legend()
plt.savefig("var.pdf")

plt.figure()
plt.plot(np.linspace(0, t/0.02, t), mean, label=r'$\Sigma_n(\Delta\xi)\xi_n\vert\psi_n(t)\vert^2$')
plt.legend()
plt.savefig("mean.pdf")
