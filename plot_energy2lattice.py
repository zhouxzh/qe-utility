#!/usr/bin/env python3

import os
import re
import numpy as np
import matplotlib.pyplot as plt

Emax = 5
Emin = -5

def filter_folder():
    my_folder = []
    lattice = []
    for i in os.listdir('.'):
        if os.path.isdir(i):
            if re.search('\d+(\.\d*)?', i):
                my_folder.append(i)
    my_folder.sort()
    return my_folder
           

def find_energy(folder):
    for my_file in os.listdir(folder):
        if 'relax.out' in my_file:
            with open(folder+'/'+my_file) as f:
                for line in f:
                    if '!' in line:
                        print(line)
                        energy = re.findall('-?\d+\.\d*', line)
                f.close()
    return float(energy[0])



my_folder = filter_folder()
print(my_folder)
energy = []
if my_folder:
    for i in my_folder:
        energy.append(find_energy(i))


lattice = np.zeros(len(my_folder))
energy = np.array(energy)
for i in range(len(my_folder)):
    lattice[i] = float(my_folder[i])
    print(lattice[i], energy[i])
print(lattice, energy)
plt.figure()
plt.plot(lattice, energy)
plt.xlabel('Lattice Constant (bohr)')
plt.ylabel('Energy (Ry)')
plt.tight_layout()
plt.savefig('lattce2energy.pdf')
