#!/usr/bin/env python3

#This file is used to fitting the Birch-Murnaghan equation.
#Dr. Xianzhong Zhou
#11.08.2018

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
plt.style.use('seaborn-paper')
mpl.rcParams['figure.figsize'] = (4.8, 3.2)
rc('text', usetex=True)
b2a = 0.52918
ry2ev = 13.6056923
lattice, energy = np.genfromtxt('./lattice.txt')
lattice = lattice * b2a
energy = energy * ry2ev
plt.plot(lattice, energy, '*')
plt.xlabel('Lattice Constant (\AA)')
plt.ylabel('Energy (eV)')

f = lambda a, E0, a0, B0, B1 : E0+9*a0**3*B0/16*(((a0/a)**(2)-1)**3*B1+((a0/a)**2-1)**2*(6-4*(a0/a)**(2)))
E0 = min(energy)
a0 = lattice[np.argmin(energy)]
print(E0, a0)

quadratic = lambda a, E0, a0, B0 : E0+9/2*a0*B0*(a-a0)**2

xid = np.argmin(energy)
offset = 3
x = lattice[xid-offset:xid+offset+1]
print(x)
y = energy[xid-offset:xid+offset+1]
popt, pcov = curve_fit(quadratic, x, y, p0=(E0, a0, 1))
print(popt)
#  print(pcov)
E0, a0, B0 = popt
print(E0, a0, B0)
x = np.linspace(min(lattice), max(lattice), 100)
y = quadratic(x, *popt)
#plt.plot(x,y, '--')
B1 = 0
x = np.linspace(min(lattice), max(lattice), 100)
y = f(x, E0, a0, B0, B1)
#plt.plot(x,y, '-')
popt, pcov = curve_fit(f, lattice, energy, p0=(E0, a0, B0, B1))
print(popt)
print(pcov)
x = np.linspace(5.8, 6.4, 100)
y = f(x, *popt)
plt.plot(x,y)
plt.tight_layout()
#plt.savefig('../lattice2energy.pdf')

plt.figure()
plt.plot(lattice, energy, '*')
plt.plot(x, y, label="with SOC")
lattice, energy = np.genfromtxt('./lattice_no_spin.txt')
lattice = lattice * b2a
energy = energy * ry2ev
popt, pcov = curve_fit(f, lattice, energy, p0=(E0, a0, B0, B1))
print(popt)
print(pcov)
y = f(x, *popt)
plt.plot(x, y, label="without SOC")
plt.plot(lattice, energy, '^')
plt.legend()
plt.show()
