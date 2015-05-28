import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
line, = plt.plot(np.linspace(-10, 10, 201), psi2[:, 0], 'b-')
ttl = plt.title("t=0")

def animate(i):
	line.set_xdata(np.linspace(-10, 10, 201))
	line.set_ydata(psi2[:, i])
	ttl = plt.title("t=" + str(i*d_t))
	return line,


def init():
	line.set_ydata(np.ma.array(psi2[:, 0], mask=True))
	plt.xlim(-10, 10)
	plt.ylim(0,.8)
	return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0,t-1,1), interval=200, init_func=init)
ani.save("Welle.gif", writer='imagemagick', fps=30)

norm = d_l * np.sum(psi2, axis=0)
mean = d_l * np.sum( (np.linspace(-10, 10, N)*psi2.T).T , axis=0)
var  = d_l * np.sum( (np.linspace(-10, 10, N)**2 * psi2.T).T , axis=0) - mean

plt.figure()
plt.plot(np.linspace(-10, 10, N), psi2[:, 0], 'b-', label=r'$\vert\psi_n(t=0)\vert^2$')
plt.plot(np.linspace(-10, 10, N), psi2[:, t-1], 'r-', label=r'$\vert\psi_n(t=' + str(t*d_t) + r')\vert^2$')
plt.legend()
plt.savefig("t0_10.pdf")

plt.figure()
plt.plot(np.linspace(0, t*d_t, t), norm, label=r'$\Sigma_n(\Delta\xi)\vert\psi_n(t)\vert^2$')
plt.ylim(0.9999,1.0001)
plt.legend()
plt.savefig("norm.pdf")

plt.figure()
plt.plot(np.linspace(0, t*d_t, t), var, label=r'$\langle\xi^2\rangle(t) - \langle\xi\rangle^2(t)$')
# plt.ylim(0.9,1.1)
plt.legend()
plt.savefig("var.pdf")

plt.figure()
plt.plot(np.linspace(0, t*d_t, t), mean, label=r'$\Sigma_n(\Delta\xi)\xi_n\vert\psi_n(t)\vert^2$')
plt.legend()
plt.savefig("mean.pdf")
