#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np


EF2 = np.genfromtxt("E_Ferro_2.txt", skip_header=1)
EF5 = np.genfromtxt("E_Ferro_5.txt", skip_header=1)
EF10 = np.genfromtxt("E_Ferro_10.txt", skip_header=1)

EAF2 = np.genfromtxt("E_Antiferro_2.txt", skip_header=1)
EAF5 = np.genfromtxt("E_Antiferro_5.txt", skip_header=1)
EAF10 = np.genfromtxt("E_Antiferro_10.txt", skip_header=1)

TF2 = np.genfromtxt("T_Ferro_2.txt", skip_header=1)
TF5 = np.genfromtxt("T_Ferro_5.txt", skip_header=1)
TF10 = np.genfromtxt("T_Ferro_10.txt", skip_header=1)

TAF2 = np.genfromtxt("T_Antiferro_2.txt", skip_header=1)
TAF5 = np.genfromtxt("T_Antiferro_5.txt", skip_header=1)
TAF10 = np.genfromtxt("T_Antiferro_10.txt", skip_header=1)


plt.plot(EF2[:, 0], EF2[:, 1], 'b-', label='N=2')
plt.plot(EF5[:, 0], EF5[:, 1], 'r-', label='N=5')
plt.plot(EF10[:, 0], EF10[:, 1], 'g-', label='N=5')
plt.legend(loc='upper right')
plt.xlabel(r'$\theta$ [$^\circ$]')
plt.ylabel(r'E($\theta$)')
plt.title('Energy ferromagnetic')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.savefig('E_Ferro.pdf')
plt.clf()

plt.plot(EAF2[:, 0], EAF2[:, 1], 'b-', label='N=2')
plt.plot(EAF5[:, 0], EAF5[:, 1], 'r-', label='N=5')
plt.plot(EAF10[:, 0], EAF10[:, 1], 'g-', label='N=10')
plt.legend(loc='upper right')
plt.xlabel(r'$\theta$ [$^\circ$]')
plt.ylabel(r'E($\theta$)')
plt.title('Energy antiferromagnetic')
plt.subplots_adjust(left=0.14)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.savefig('E_Antiferro.pdf')
plt.clf()

plt.plot(TF2[:, 0], TF2[:, 1], 'b-', label='N=2')
plt.plot(TF5[:, 0], TF5[:, 1], 'r-', label='N=5')
plt.plot(TF10[:, 0], TF10[:, 1], 'g-', label='N=10')
plt.legend(loc='upper right')
plt.xlabel(r'$\theta$ [$^\circ$]')
plt.ylabel(r'T($\theta$)')
plt.title('Momentum ferromagnetic')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.savefig('T_Ferro.pdf')
plt.clf()

plt.plot(TAF2[:, 0], TAF2[:, 1], 'b-', label='N=2')
plt.plot(TAF5[:, 0], TAF5[:, 1], 'r-', label='N=5')
plt.plot(TAF10[:, 0], TAF10[:, 1], 'g-', label='N=10')
plt.legend(loc='upper right')
plt.xlabel(r'$\theta$ [$^\circ$]')
plt.ylabel(r'T($\theta$)')
plt.title('Momentum antiferromagnetic')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.savefig('T_Antiferro.pdf')
plt.clf()
