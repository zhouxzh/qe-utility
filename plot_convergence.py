#!/opt/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os

cutoff=[]
energy=[]
for out_file in os.listdir('./'):
    if out_file.find('out')>=0:
        print(out_file)
        with open(out_file) as f:
            for i in f:
                if i.find('kinetic-energy cutoff')>=0:
                    cutoff.append(float(i.split()[3]))
                if i.find('!')>=0:
                    energy.append(float(i.split()[4]))

dtype=[('cutoff',float),('energy',float)]
cutoff=np.array(cutoff)
energy=np.array(energy)
energy=energy[np.argsort(cutoff)]
cutoff=cutoff[np.argsort(cutoff)]
plt.plot(cutoff, energy)
plt.xlabel('ecutwfc (Ry) (ecutrho=4*ecutwfc)')
plt.ylabel('Total Energy (Ry)')
plt.tight_layout()
plt.savefig('convergence.pdf')
for i in range(len(cutoff)):
    print(cutoff[i], energy[i])
