
# coding: utf-8

plt.plot(x, y, 'bo')
plt.quiver(x, y, v_x, v_y, pivot='middle', linestyle='dashed')

plt.xlim(-L, 2*L)
plt.ylim(-L, 2*L)

plt.hlines([0, L], -L, 2*L, linestyles='dashed')
plt.vlines([0, L], -L, 2*L, linestyles='dashed')

for m in [-L,0,L]:
    for n in [-L,0,L]:
        if [n,m] is not [0,0]:
            plt.plot(x+n, y+m, 'bo', alpha=0.5)
            plt.quiver(x+n, y+m, v_x, v_y, pivot='middle', alpha=0.5, linestyle='dashed')

plt.savefig('init.pdf')
