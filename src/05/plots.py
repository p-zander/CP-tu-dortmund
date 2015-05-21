# coding: utf-8

plt.plot(loc_sample[0,:, 0],   loc_sample[0,:, 1], 'bo') 

plt.quiver(loc_sample[0,:, 0], loc_sample[0,:, 1], vel_sample[0,:, 0], vel_sample[0,:, 1], pivot='middle', headlength=2, headaxislength=1) 

plt.xlim(0, L)
plt.ylim(0, L)

plt.savefig('init' + file_ext)

# for i in range(N):
#     plt.plot(loc_sample[:, i, 0], loc_sample[:, i, 1], ',')
# plt.savefig('time' + file_ext)

# fig = plt.figure()
# line, = plt.plot(loc_sample[0, :, 0], loc_sample[0,:, 1], 'bo')

# def animate(i):
# line.set_xdata(loc_sample[i, :, 0])  # update the data
#     line.set_ydata(loc_sample[i, :, 1])
#     return line,

# def init():
#     line.set_ydata(np.ma.array(loc_sample[0, :, 0], mask=True))
#     plt.xlim(0, 8)
#     plt.ylim(0, 8)
#     return line,

# ani = animation.FuncAnimation(fig, animate, np.arange(0,len(loc_sample),samples/100), init_func=init, interval=25)
# ani.save("positions.gif", writer='imagemagick')

plt.figure()

time = np.linspace(0, t_max, samples)

plt.plot(time, E_kin_sample, label=u'$E_{kin}$')
plt.plot(time, E_pot_sample, label=u'$E_{pot}$')
plt.plot(time, E_pot_sample + E_kin_sample, label=u'$E_{tot}$')
plt.plot(time, 1 / 2 * cm_vel_sample, label=u'$E_{kin,CM}$')

plt.vlines(t_aqui, -30, 30)

plt.legend(loc='best')

plt.savefig('energy' + file_ext)

plt.figure()
plt.bar(np.arange(0, cutoff, cutoff / bins), radial_sample, cutoff / bins)

plt.savefig('g_' + file_ext)

plt.figure()

plt.plot(time, E_kin_sample/N, label=u'$T$')
plt.vlines(t_aqui, 0, 2)

plt.legend(loc='best')

plt.savefig('temperature' + file_ext)

print 'mean temperature after equilibration was: ', np.mean(E_kin_sample[int(samples*t_aqui/t_max + 1):] / N)

# T_init = np.array([1.000, 0.800, 0.600, 0.500, 0.300, 0.250, 0.001, 1e-05, 1e-10])
# T_mean = np.array([1.4467, 1.2614, 1.0957, 1.0072, 0.8508, 0.8183, 0.6588, 0.6519, 0.6474])

# plt.figure()
# plt.semilogx(T_init, T_mean)

# plt.savefig('T_init_T_mean_curve' + file_ext)
