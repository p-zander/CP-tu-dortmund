
# coding: utf-8

from scipy import ndimage
print 'importing matplotlib...'
from matplotlib import pyplot as plt
from matplotlib import cm
plt.style.use('ggplot')
print 'done'

print 'plotting...'

# –– gradient field –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
dx = ndimage.sobel(A, 0)  # horizontal derivative
dy = ndimage.sobel(A, 1)  # vertical derivative

# –– angles at boundaries –––––––––––––––––––––––––––––––––––––––––––––––––––––
plt.figure()

plt.plot(np.abs(np.arctan2(dx[:, 0], dy[:, 0])), label=u'at $y = 0$')
plt.plot(np.abs(np.arctan2(dx[0, :], dy[0, :])), label=u'at $x = 0$')
plt.plot(np.abs(np.arctan2(dx[:, -1], dy[:, -1])), label=u'at $y = 1$')
plt.plot(np.abs(np.arctan2(dx[-1, :], dy[-1, :])), label=u'at $x = 1$')

print 'Charge along x=0 axis: ', sum(dx[0, :])
print 'Charge along x=1 axis: ', sum(dx[-1, :])
print 'Charge along y=0 axis: ', sum(dy[:, 0])
print 'Charge along y=1 axis: ', sum(dy[:, -1])

plt.legend(loc='best')
plt.savefig('angles.pdf')

# –– plot potential and field –––––––––––––––––––––––––––––––––––––––––––––––––
fig = plt.figure(figsize=[8, 8])

dx = np.reshape(dx, N * N)
dy = np.reshape(dy, N * N)

points = np.array(list(np.ndindex(N, N)))

ax = fig.gca()
ax.imshow(A.T, cmap='seismic')
ax.quiver(points[:, 0], points[:, 1], -dx, -dy, angles='xy', pivot='middle')

fig.savefig('A.pdf')
