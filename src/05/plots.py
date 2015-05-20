# coding: utf-8

plt.plot(loc_sample[0, :, 0],   loc_sample[0,:, 1], 'bo') 

plt.quiver(loc_sample[0, :, 0], loc_sample[0,:, 1], vel_sample[0,:, 0], vel_sample[0,:, 1], scale=velocity_scale, pivot='middle', headlength=2, headaxislength=1) 

plt.xlim(0, L)
plt.ylim(0, L)

plt.savefig('init' + file_ext)

# for i in range(N):
#     plt.plot(loc_sample[:, i, 0], loc_sample[:, i, 1], ',')
# plt.savefig('time' + file_ext)

# fig = plt.figure()
# line, = plt.plot(loc_sample[0, :, 0], loc_sample[0,:, 1], 'bo')

# def animate(i):
#     line.set_xdata(loc_sample[i, :, 0])  # update the data
#     line.set_ydata(loc_sample[i, :, 1])
#     return line,

# def init():
#     line.set_ydata(np.ma.array(loc_sample[0, :, 0], mask=True))
#     plt.xlim(0, 8)
#     plt.ylim(0, 8)
#     return line,

# ani = animation.FuncAnimation(fig, animate, np.arange(0,len(loc_sample),samples/100), init_func=init, interval=25)
# ani.save("time.gif", writer='imagemagick')

plt.figure()

time = np.arange(0, t_max, sbs)

plt.plot(E_kin_sample, label=u'$E_{kin}$')
plt.plot(E_pot_sample, label=u'$E_{pot}$')
plt.plot(E_pot_sample+E_kin_sample, label=u'$E_{tot}$')
plt.plot(cm_vel_sample, label=u'$v_{CM}$')
plt.ylim(-40, 40)
plt.legend(loc='best')

# plt.plot(cm_vel_sample[:len(E_kin_sample)])

plt.savefig('energy' + file_ext)
