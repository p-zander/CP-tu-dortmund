
# coding: utf-8
from matplotlib import pyplot as plt

plt.figure(figsize=[8, 8])
plt.plot(x, y, 'bo')
circle = plt.Circle((x[0], y[0]), cutoff, color='k', fill=False, linestyle='dashed')
plt.quiver(x, y, v_x, v_y, pivot='middle', headlength=2, headaxislength=1)
plt.quiver(x, y, F_x, F_y, color='red', pivot='middle',
           headlength=2, headaxislength=1)

plt.xlim(-L, 2 * L)
plt.ylim(-L, 2 * L)

plt.hlines([0, L], -L, 2 * L, linestyles='dashed')
plt.vlines([0, L], -L, 2 * L, linestyles='dashed')

for m in [-L, 0, L]:
    for n in [-L, 0, L]:
        if n is not 0 or m is not 0:
            plt.plot(x + n, y + m, 'bo', alpha=0.5)
            plt.quiver(x + n, y + m, v_x, v_y, pivot='middle',
                       alpha=0.5, headlength=2, headaxislength=1)

plt.gca().add_artist(circle)

plt.savefig('init.pdf')
