#! /usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt("2a.txt")
b = np.genfromtxt("2b.txt")

plt.plot(a[:, 0], a[:, 1], label=u'$x_0$')
plt.plot(b[:, 0], b[:, 1], label=u'$x_1$')
plt.plot(b[:, 0], b[:, 2], label=u'$y_1$')
plt.xlabel("t")
plt.ylabel("x/y(t)")
plt.legend(loc='best', fontsize='large')
plt.savefig("2a.pdf")
plt.clf()

plt.plot(b[:, 0], 1/2 * np.sum(b[:, 4:7]**2, axis=1), label=u'$E_{kin}$')
plt.plot(b[:, 0], 1/2 * np.sum(b[:, 1:4]**2, axis=1), label=u'$E_{pot}$')
plt.plot(b[:, 0], 1/2 * np.sum(b[:, 4:7]**2, axis=1) + 1/2 * np.sum(b[:, 1:4]**2, axis=1), label=u'$E_{ges}$')
plt.xlabel("t")
plt.ylabel("E")
plt.ylim(0, 1)
plt.legend(loc='best', fontsize='large')
plt.savefig("2b.pdf")
plt.clf()

E_a = 1/2 * np.sum(a[:, 4:7]**2, axis=1) + 1/2 * np.sum(a[:, 1:4]**2, axis=1)
E_b = 1/2 * np.sum(b[:, 4:7]**2, axis=1) + 1/2 * np.sum(b[:, 1:4]**2, axis=1)
plt.semilogy(b[:, 0], np.abs(E_b - E_b[0]), label=u'$E_{ges}$')
plt.semilogy(a[:, 0], np.abs(E_a - E_a[0]), label=u'$E_{ges}$')
plt.xlabel("t")
plt.ylabel(u"$\Delta E$")
plt.ylim(1e-6, 1)
plt.legend(loc='upper left', fontsize='large')
plt.savefig("2b_2.pdf")
plt.clf()

del a, b, E_a, E_b

a1 = np.genfromtxt("a_tn_37_N_25_vy_1.3.txt")
a2 = np.genfromtxt("a_tn_37_N_100_vy_1.3.txt")
a3 = np.genfromtxt("a_tn_37_N_100_vy_0.88.txt")
bcd = np.genfromtxt("bcd.txt")
e1 = np.genfromtxt("e_alpha_0.9.txt")
e2 = np.genfromtxt("e_alpha_1.1.txt")

# Three Plots with different N and v_y for problem a.
plt.plot(a1[:, 1], a1[:, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-4, 2)
plt.ylim(-3, 3)
plt.savefig("a1.pdf")
plt.clf()

plt.plot(a2[:, 1], a2[:, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-6, 1.2)
plt.ylim(-3.5, 3.7)
plt.savefig("a2.pdf")
plt.clf()

plt.plot(a3[:, 1], a3[:, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-0.8, 1.2)
plt.ylim(-1, 1)
plt.savefig("a3.pdf")
plt.clf()

# All plots for b) with two separate legends
plt.figure(figsize=[10, 6])
Ldot, = plt.plot(bcd[:, 0], bcd[:, 7], label=r'$\dot{\mathrm{L}}$')
Lx, = plt.plot(bcd[:, 0], bcd[:, 8], label=r'$\mathrm{L}_x$')
Ly, = plt.plot(bcd[:, 0], bcd[:, 9], label=r'$\mathrm{L}_y$')
Lz, = plt.plot(bcd[:, 0], bcd[:, 10], label=r'$\mathrm{L}_z$')
Ek, = plt.plot(
    bcd[:, 0], .5 * np.sum(bcd[:, 4:7]**2, axis=1), label=r'$\mathrm{E}_{kin}$')
V, = plt.plot(
    bcd[:, 0], -1 / np.sqrt(np.sum(bcd[:, 1:4]**2, axis=1)), label=r'V')
E, = plt.plot(bcd[:, 0], .5 * np.sum(bcd[:, 4:7]**2, axis=1) - 1 /
              np.sqrt(np.sum(bcd[:, 1:4]**2, axis=1)), label=r'$\mathrm{E}_{ges}$')

legend2 = plt.legend(handles=[V, Ek, E], loc=1)
ax = plt.gca().add_artist(legend2)

plt.legend(handles=[Ldot, Lx, Ly, Lz], loc=2)
plt.ylim(-1.1, 2.5)
plt.xlabel('t')
plt.savefig("b.pdf")
plt.clf()

# First Plot: Lenz-Runge components. second plot: ellipse witch LR-vector
plt.plot(bcd[:, 0], bcd[:, 11], label=r'$\mathrm{LR}_x$')
plt.plot(bcd[:, 0], bcd[:, 12], label=r'$\mathrm{LR}_y$')
plt.plot(bcd[:, 0], bcd[:, 13], label=r'$\mathrm{LR}_z$')
plt.legend(loc=6)
plt.xlabel('t')
plt.ylim(-0.01, 0.71)
plt.savefig("c1.pdf")
plt.clf()

plt.gca().set_aspect('equal', adjustable='box')
plt.plot(bcd[:, 1], bcd[:, 2])
arr = plt.arrow(0.31, 0, 0.69, 0, head_width=0.05,
                head_length=0.1, fc='r', ec='r', length_includes_head=True)
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-6, 1.2)
plt.ylim(-3.5, 3.7)
plt.savefig("c2.pdf")
plt.clf()

# Plots the difference ru-r and vu-v to show equivalency
plt.plot(bcd[:, 0], bcd[:, 1] - bcd[::-1, -6], label=r'$r_x - ru_x$')
plt.plot(bcd[:, 0], bcd[:, 2] - bcd[::-1, -5], label=r'$r_y - ru_y$')
plt.plot(bcd[:, 0], bcd[:, 3] - bcd[::-1, -4], label=r'$r_z- ru_z$')
plt.plot(bcd[:, 0], np.abs(bcd[:, 4]) - np.abs(bcd[::-1, -3]),
         label=r'$\vert v_x\vert -\vert vu_x\vert$')
plt.plot(bcd[:, 0], np.abs(bcd[:, 5]) - np.abs(bcd[::-1, -2]),
         label=r'$\vert v_y\vert -\vert vu_y\vert$')
plt.plot(bcd[:, 0], np.abs(bcd[:, 6]) - np.abs(bcd[::-1, -1]),
         label=r'$\vert v_z\vert -\vert vu_z\vert$')
plt.legend(loc='best')
plt.xlabel('t')
plt.savefig('d.pdf')
plt.clf()

# Plots for e) LR-vector and ellipse alpha 0.9
plt.plot(e1[:, 0], e1[:, 4], label=r'$\mathrm{LR}_x$')
plt.plot(e1[:, 0], e1[:, 5], label=r'$\mathrm{LR}_y$')
plt.plot(e1[:, 0], e1[:, 6], label=r'$\mathrm{LR}_z$')
plt.plot(e1[:, 0], np.sqrt(e1[:, 4]**2 + e1[:, 5]**2 + e1[:, 6]**2),
         label=r'$\vert \mathrm{LR} \vert$')
plt.legend(loc=5)
plt.xlabel('t')
plt.ylim(-0.7, 0.85)
plt.savefig("e09LR.pdf")
plt.clf()

plt.plot(e1[:, 1], e1[:, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-8, 2)
plt.ylim(-4, 6)
plt.savefig("e09EL.pdf")
plt.clf()

# Plots for e) LR-vector and ellipse alpha 1.1
plt.plot(e2[:, 0], e2[:, 4], label=r'$\mathrm{LR}_x$')
plt.plot(e2[:, 0], e2[:, 5], label=r'$\mathrm{LR}_y$')
plt.plot(e2[:, 0], e2[:, 6], label=r'$\mathrm{LR}_z$')
plt.plot(e2[:, 0], np.sqrt(e2[:, 4]**2 + e2[:, 5]**2 + e2[:, 6]**2),
         label=r'$\vert \mathrm{LR} \vert$')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylim(-0.7, 0.85)
plt.savefig("e11LR.pdf")
plt.clf()

plt.plot(e2[:, 1], e2[:, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-6, 2)
plt.ylim(-5, 3)
plt.savefig("e11EL.pdf")
plt.clf()
