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

plt.plot(E_kin_sample)
plt.plot(E_pot_sample)
plt.plot(cm_vel_sample)
plt.plot(E_kin_sample[0:-1] + E_pot_sample[1:])

plt.savefig('energy' + file_ext)
