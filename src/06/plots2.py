import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
line, = plt.plot(np.linspace(-10,10,201), psim[:,0], 'b-')


def animate(i):
    line.set_xdata(np.linspace(-10,10,201))  # update the data
    line.set_ydata(psim[:,i])
    return line,


def init():
    line.set_ydata(np.ma.array(psim[:, 0], mask=True))
    plt.xlim(-10,10)
    plt.ylim(0,.1)
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(
    0, 499, 1), init_func=init, interval=25)
ani.save("Welle.gif", writer='imagemagick', fps=10)
