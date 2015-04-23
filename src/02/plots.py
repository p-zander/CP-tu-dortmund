#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np


V1_in = np.genfromtxt("V1_in.txt", skip_header=1)
V1_out = np.genfromtxt("V1_out.txt", skip_header=1)
V2_out = np.genfromtxt("V2_out.txt", skip_header=1)


plt.plot(V1_in[:, 0], V1_in[:, 1], 'b-', label='V1 inside')
plt.plot(V1_out[:, 0], V1_out[:, 1], 'r-', label='V1 outside')
plt.plot(V1_out[:, 0], V1_out[:, 2], 'g-', label='V1 monopol')
plt.legend(loc='best')
plt.xlabel(r'Distance x [$\frac{x}{a}$] ')
plt.ylabel(r'V1')
plt.title('Potential 1')
plt.savefig('V1.pdf')
plt.clf()

plt.plot(V2_out[:, 0], V2_out[:, 1], 'b-', label='V2 outside')
plt.plot(V2_out[:, 0], V2_out[:, 2], 'g-', label='V2 Dipol')
plt.legend(loc='best')
plt.xlabel(r'Distance x [$\frac{x}{a}$] ')
plt.ylabel(r'V2')
plt.title('Potential 2')
plt.savefig('V2.pdf')
plt.clf()
