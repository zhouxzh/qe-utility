#!/opt/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
from io import BytesIO
from io import StringIO
from scipy.interpolate import CubicSpline
from scipy.optimize import fmin
import re

b2a = 0.52917706


def plot_energy_step(energy_step):
    plt.figure()
    plt.plot(energy_step,label='Lattice Constant = '+my_dir+' Bohr')
    plt.xlabel('Iteration Step')
    plt.ylabel('Energy (Ry)')
    plt.tight_layout()
    plt.legend()
    plt.savefig('relax_energy.pdf')

def plot_lattice_energy(lattice, energy):
    lattice = np.array(lattice)
    energy = np.array(energy)
    energy = energy[np.argsort(lattice)]
    lattice = np.sort(lattice)
    plt.figure()
    plt.plot(lattice*b2a, energy, '*')
    cs = CubicSpline(lattice, energy)
    x = np.linspace(min(lattice),max(lattice), 100)
    plt.plot(x*b2a, cs(x))
    plt.xlabel(r'Lattice Constant ($\AA$)')
    plt.ylabel('Energy (Ry)')
    plt.tight_layout()
    plt.savefig('lattice_energy.pdf')

def search_energy():
    lattice=[]
    energy=[]
    for my_dir in os.listdir('./'):
        if os.path.isdir(my_dir):
            if re.search('\d+(\.\d*)?', my_dir):
                lattice.append(float(my_dir))
            else:
                continue
            os.chdir(my_dir)
            for out_file in os.listdir('./'):
                if out_file.find('relax.out')>=0:
                    print(my_dir, out_file)
                    with open(out_file) as f:
                        energy_step=[]
                        for i in f:
                            if i.find('!')>=0:
                                energy_step.append(float(i.split()[4]))
                        #plot_energy_step(energy_step)
                        energy.append(energy_step[-1])
            os.chdir('../')
    return lattice, energy

def plot_sigle_band(data):
    plt.figure()
    circle = np.argwhere(data[:,0] == data[0,0])[:,0]
    for i in range(len(circle)-1):
        plt.plot(data[circle[i]:circle[i+1],0], data[circle[i]:circle[i+1],1])
    plt.ylim(-5, 5)
    plt.xlim(0, 2.5231)
    positions = [0, 0.8660, 1.366, 1.816, 2.5231]
    labels = [r'$\Gamma$', 'R', 'M', 'X', 'R']
    plt.xticks(positions, labels)
    plt.ylabel('Energy (eV)')
    plt.savefig('bands.pdf')
    plt.close()

def plot_bm_band(data):
    circle = np.argwhere(data[:,0] == data[0,0])[:,0]
    vbm_up = data[circle[48]:circle[49],:]
    vbm_down = data[circle[49]:circle[50],:]
    cbm_up = data[circle[50]:circle[51],:]
    cbm_down = data[circle[51]:circle[52],:]
    plt.figure()
    plt.plot(vbm_up[:,0],vbm_up[:,1],label='VBM up')
    plt.plot(vbm_down[:,0],vbm_down[:,1], label='VBM down')
    plt.plot(cbm_up[:,0],cbm_up[:,1], label='CBM up')
    plt.plot(cbm_down[:,0],cbm_down[:,1], label='CBM down')
    plt.xlim(0, 2.5231)
    positions = [0, 0.8660, 1.366, 1.816, 2.5231]
    labels = [r'$\Gamma$', 'R', 'M', 'X', 'R']
    plt.xticks(positions, labels)
    plt.ylabel('Energy (eV)')
    plt.legend()
    plt.savefig('bm.pdf')
    plt.close()

def plot_all_band():
    lattice = []
    for my_dir in os.listdir():
        if os.path.isdir(my_dir):
            os.chdir(my_dir)
            lattice.append(float(my_dir))
            print(my_dir)
            for out_file in os.listdir():
                if out_file.find('relax.out')>=0:
                    print(out_file)
                    with open(out_file) as f:
                        for i in f:
                            if i.find('highest')>=0:
                                homo, lomo = i.split()[6::]
                                fermi = (float(homo)+float(lomo))/2
                    print(fermi)
            for out_file in os.listdir():
                if out_file.find('bands.dat.gnu')>=0:
                    print(out_file)
                    data = np.genfromtxt(out_file)
                    data[:,1] = data[:,1]-fermi
                    plot_sigle_band(data)
                    plot_bm_band(data)
            os.chdir('../')

def rashba(lower, upper):
    lattice = []
    ER = []
    k0 = []
    for my_dir in os.listdir():
        if os.path.isdir(my_dir):
            if re.search('\d+(\.\d*)?', my_dir):
                lattice.append(float(my_dir))
            else:
                continue
            os.chdir(my_dir)
            print(my_dir)
            for out_file in os.listdir():
                if out_file.find('relax.out')>=0:
                    print(out_file)
                    with open(out_file) as f:
                        for i in f:
                            if i.find('highest')>=0:
                                homo, lomo = i.split()[6::]
                                fermi = (float(homo)+float(lomo))/2
                    print(fermi)
            for out_file in os.listdir():
                if out_file.find('bands.dat.gnu')>=0:
                    print(out_file)
                    data = np.genfromtxt(out_file)
                    data[:,1] = data[:,1]-fermi
                    circle = np.argwhere(data[:,0] == data[0,0])[:,0]
                    vbm_up = data[circle[48]:circle[49],1]
                    vbm_down = data[circle[49]:circle[50],1]
                    cbm_up = data[circle[50]:circle[51],1]
                    cbm_down = data[circle[51]:circle[52],1]
                    k = data[circle[51]:circle[52],0]
                    k_range = np.bitwise_and( k>= lower, k<= upper )
                    vbm_up = vbm_up[k_range]
                    vbm_down = vbm_down[k_range]
                    cbm_up = cbm_up[k_range]
                    cbm_down = cbm_down[k_range]
                    k = k[k_range]
                    print(k)
                    ER_vbm_up = max(vbm_up)-vbm_up[k==0.866]
                    ER_vbm_down = max(vbm_down)-vbm_down[k==0.866]
                    ER_cbm_up = min(cbm_up)-cbm_up[k==0.866]
                    ER_cbm_down = min(cbm_down)-cbm_down[k==0.866]
                    k0_vbm_up = k[np.argmax(vbm_up)]-0.866
                    k0_vbm_down = k[np.argmax(vbm_down)]-0.866
                    k0_cbm_up = k[np.argmin(cbm_up)]-0.866
                    k0_cbm_down = k[np.argmin(cbm_down)]-0.866
                    print(ER_vbm_up, ER_vbm_down, ER_cbm_up, ER_cbm_down)
                    print(k0_vbm_up, k0_vbm_down, k0_cbm_up, k0_cbm_down)
                    cs = CubicSpline(k, cbm_up)
                    minimum = fmin(cs, 0.866)
                    k0.append(minimum[0]-0.866)
                    ER.append(cs(minimum[0])-cs(0.866))
            os.chdir('../')
    lattice = np.array(lattice)
    ER = np.array(ER)
    k0 = np.array(k0)
    ER = ER[np.argsort(lattice)]
    k0 = k0[np.argsort(lattice)]
    lattice = np.sort(lattice)
    return lattice, np.abs(k0), np.abs(ER)

def plot_rashba(lower, upper, prefix):
    lattice, k0, ER = rashba(lower, upper)
    plt.figure()
    lattice = lattice * b2a
    k0 = k0*2*np.pi/lattice
    plt.plot(lattice, ER)
    plt.xlabel('Lattice Constant ($\AA$)')
    plt.ylabel('Rashba Energy (eV)')
    plt.tight_layout()
    plt.savefig(prefix+'_ER.pdf')
    plt.close()
    plt.figure()
    plt.xlabel('Lattice Constant ($\AA$)')
    plt.ylabel('Momentum offset (${\AA}^{-1}$)')
    plt.plot(lattice, k0, label='VBM up')
    plt.tight_layout()
    plt.savefig(prefix+'_k0.pdf')
    plt.close()
    plt.figure()
    plt.plot(lattice, ER/(2*k0))
    plt.xlabel('Lattice Constant ($\AA$)')
    plt.ylabel('Rashba coefficient (eV$\AA$)')
    plt.tight_layout()
    plt.savefig(prefix+'_alfa.pdf')
    plt.close()




#lattice, energy = search_energy()
#plot_lattice_energy(lattice, energy)
#plot_all_band()



##### search gG-R #####
lower = 0
upper = 0.866
prefix = 'GR'
plot_rashba(lower, upper, prefix)

##### search R-M #####
lower = 0.866
upper = 1.36
prefix = 'RM'
plot_rashba(lower, upper, prefix)
