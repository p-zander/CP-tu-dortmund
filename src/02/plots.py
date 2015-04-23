#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

V1 = np.genfromtxt("V1.txt", skip_header=1)
V2 = np.genfromtxt("V2.txt", skip_header=1)

plt.plot(V1[:, 0], V1[:, 1], 'r-', label='V1 numerisch')
plt.plot(V1[:, 0], V1[:, 2], 'g-', label='V1 asymptotisch (Monopol)')
plt.legend(loc='best')
plt.xlabel(r'Distance x [$\frac{x}{a}$] ')
plt.ylabel(r'[$V\,(4\pi\varepsilon_0)\,(\rho_0 a^2)^{-1}$]')
plt.ylim(0, 10)
plt.title('Potential 1')
plt.savefig('V1.pdf')
plt.clf()

plt.plot(V2[:, 0], V2[:, 1], 'b-', label='V2 numerisch')
plt.plot(V2[:, 0], V2[:, 2], 'g-', label='V2 asymptotisch (Dipol)')
plt.legend(loc='best')
plt.xlabel(r'Distance x [$\frac{x}{a}$] ')
plt.ylabel(r'[$V\,(4\pi\varepsilon_0)\,(\rho_0 a^2)^{-1}$]')
plt.ylim(0, 2)
plt.title('Potential 2')
plt.savefig('V2.pdf')
plt.clf()
