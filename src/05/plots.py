
# coding: utf-8
from matplotlib import pyplot as plt

plt.figure(figsize=[8, 8])
plt.plot(x_0, y_0, 'bo')
circle = plt.Circle((x_0[0], y_0[0]), cutoff, color='k', fill=False, linestyle='dashed')
plt.quiver(x_0, y_0, v_x_0, v_y_0, pivot='middle', headlength=2, headaxislength=1)
plt.quiver(x_0, y_0, F_x_0, F_y_0, color='red', pivot='middle',
           headlength=2, headaxislength=1)

plt.xlim(-L, 2 * L)
plt.ylim(-L, 2 * L)

plt.hlines([0, L], -L, 2 * L, linestyles='dashed')
plt.vlines([0, L], -L, 2 * L, linestyles='dashed')

for m in [-L, 0, L]:
    for n in [-L, 0, L]:
        if n is not 0 or m is not 0:
            plt.plot(x_0 + n, y_0 + m, 'bo', alpha=0.5)
            plt.quiver(x_0 + n, y_0 + m, v_x_0, v_y_0, pivot='middle',
                       alpha=0.5, headlength=2, headaxislength=1)

plt.gca().add_artist(circle)

plt.savefig('init.pdf')

plt.figure(figsize=[8, 8])
plt.plot(x_0, y_0, 'bo')
plt.plot(x_1, y_1, 'bo')
plt.quiver(x_1, y_1, v_x_1, v_y_1, pivot='middle', headlength=2, headaxislength=1)
plt.quiver(x_1, y_1, F_x_1, F_y_1, color='red', pivot='middle',
           headlength=2, headaxislength=1)

plt.xlim(0, L)
plt.ylim(0, L)

# plt.hlines([0, L], -L, 2 * L, linestyles='dashed')
# plt.vlines([0, L], -L, 2 * L, linestyles='dashed')

# for m in [-L, 0, L]:
#     for n in [-L, 0, L]:
#         if n is not 0 or m is not 0:
#             plt.plot(x_1 + n, y_1 + m, 'bo', alpha=0.5)
#             plt.quiver(x_1 + n, y_1 + m, v_x_1, v_y_1, pivot='middle',
#                        alpha=0.5, headlength=2, headaxislength=1)


plt.savefig('first_step.pdf')
