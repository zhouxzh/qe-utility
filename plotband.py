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
           

def band(folder = '.'):
    for i in os.listdir(folder):
        if 'dat.gnu' in i:
            j = 0
            with open(folder+'/'+i) as f:
                for line in f:
                    if not line.strip():
                        break
                    j = j+1
                f.close()
            print(j)
            data = np.genfromtxt(folder+'/'+i)
            m, n = data.shape
            k = m//j
            data = data.reshape(k, j, n)
            x = data[0,:,0]
            plt.figure()
            fermi = find_fermi(folder)
            if not fermi:
                fermi = 0
            print('The Fermi energy is ',fermi, ' eV.')
            for line in range(k):
                plt.plot(x, data[line,:,1]-fermi)
            plt.xlim(min(x), max(x))
            plt.ylim(Emin, Emax)
            plt.ylabel('Energy (eV)')
            value, point = find_ticker(folder)
            print(value, point)
            plt.xticks(value, point)
            plt.tight_layout()
            plt.savefig(folder+'/bands.pdf')
            

def find_ticker(folder = '.'):
    print(folder)
    for my_file in os.listdir(folder):
        if 'bands.in' in my_file:
            print(my_file)
            point = []
            with open(folder+'/'+my_file) as f:
                for line in f:
                    if 'crystal_b' in line:
                        point_nu = int(f.readline().strip())
                        print(point_nu)
                        for i in range(point_nu):
                            point.append(f.readline().split()[0])
                        print(point)
                f.close()
        if 'plotband.out' in my_file:
            print(my_file)
            value = []
            with open(folder+'/'+my_file) as f:
                for line in f:
                    if 'symmetry' in line:
                        value.append(float(line.split()[-1]))
                f.close()
            print(value)
    return value, point


def find_fermi(folder):
    for my_file in os.listdir(folder):
        if 'scf.out' in my_file:
            print(my_file)
            fermi = []
            with open(folder+'/'+my_file) as f:
                for line in f:
                    if 'highest occupied' in line:
                        m = re.findall('\d+\.\d*', line)
                        if len(m) == 1:
                            fermi = float(m[0])
                        if len(m) == 2:
                            fermi = (float(m[0])+float(m[1]))/2
                        return fermi



my_folder = filter_folder()


if my_folder:
    for i in my_folder:
        band(i)
else:
    band()
