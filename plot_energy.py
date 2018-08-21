#!/opt/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.optimize import fmin, fmin_l_bfgs_b
import argparse
import re

plt.style.use('seaborn-paper')
mpl.rcParams['figure.figsize'] = (4.8, 3.2)
from matplotlib import rc
rc('text', usetex=True)

b2a = 0.52917706
Emax = 5
Emin = -5


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
    #plt.plot(lattice*b2a, energy, '*')
    cs = CubicSpline(lattice, energy)
    x = np.linspace(min(lattice),max(lattice), 100)
    plt.plot(x*b2a, cs(x))
    plt.xlabel(r'Lattice Constant ($\mathrm{\AA}$)')
    plt.ylabel('Energy (Ry)')
    plt.tight_layout()
    plt.savefig('lattice2energy.pdf')

def find_energy(folder):
    os.chdir(folder)
    for out_file in os.listdir():
        if re.search('^scf.out', out_file):
            print(folder, out_file)
            with open(out_file) as f:
                energy_step=[]
                for i in f:
                    if i.find('!')>=0:
                        energy_step.append(float(i.split()[4]))
            os.chdir('../')
            energy = energy_step[-1]
            return energy
        elif re.search('^relax.out', out_file):
            print(folder, out_file)
            with open(out_file) as f:
                energy_step=[]
                for i in f:
                    if i.find('!')>=0:
                        energy_step.append(float(i.split()[4]))
            os.chdir('../')
            energy = energy_step[-1]
            return energy


def filter_folder():
    my_folder = []
    lattice = []
    for i in os.listdir('.'):
        if os.path.isdir(i):
            if re.search('^\d+(\.\d*)?', i):
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
            for i in value:
                plt.plot([i,i], [Emin, Emax], color='0.75', linewidth=0.75)
            plt.tight_layout()
            plt.savefig(folder+'/bands.pdf')
            

def find_ticker(folder = '.'):
    #print(folder)
    for my_file in os.listdir(folder):
        if 'bands.in' in my_file:
            #print(my_file)
            point = []
            with open(folder+'/'+my_file) as f:
                for line in f:
                    if 'crystal_b' in line:
                        point_nu = int(f.readline().strip())
                        #print(point_nu)
                        for i in range(point_nu):
                            point.append(f.readline().split()[0])
                            if point[-1]=='gG':
                                point[-1] = '$\Gamma$'
                        #print(point)
                f.close()
        if 'plotband.out' in my_file:
            #print(my_file)
            value = []
            with open(folder+'/'+my_file) as f:
                for line in f:
                    if 'symmetry' in line:
                        value.append(float(line.split()[-1]))
                f.close()
            #print(value)
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

def find_electrons(folder = '.'):
    for i in os.listdir(folder):
        if 'scf.out' == i:
            scf = os.path.join(folder, i)
            with open(scf) as f:
                for line in f:
                    if 'number of electrons' in line:
                        print(line.split())
                        return int(float(line.split()[4]))

def find_bm_band(folder = '.'):
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
            electrons = find_electrons(folder)
            print('The number of electrons: ', electrons)
            bm = data[electrons-2:electrons+2,:,1]-fermi
            x = data[0,:,0]
            return x, bm



def plot_bm_band(folder, x, bm, value, point):
    band = ['VBM', 'CBM']
    spin = ['up', 'down']
    plt.figure()
    for i in range(2):
        for j in range(2):
            plt.plot(x, bm[i*2+j,:], label=band[i]+' spin '+spin[j])
    plt.legend()
    plt.xlim(min(x), max(x))
    plt.ylabel('Energy (eV)')
    plt.xticks(value, point)
    plt.tight_layout()
    plt.savefig(os.path.join(folder, 'bm_band.pdf'))
    plt.close()

def find_lattice_constant(folder='.'):
    for my_file in os.listdir(folder):
        if 'in' in my_file:
            with open(my_file) as f:
                for line in f:
                    if 'celldm(1)' in line:
                        print(line)
                        m = re.findall('\d+\.?\d*', line)
                        if m:
                            return float(m[1])


def find_rashba(lattice, x, bm, value, point):
    rashba = []
    band = ['VBM', 'CBM']
    spin = ['up', 'down']
    for i in range(2):
        for j in range(2):
            energy = bm[i*2+j,:]
            R = point.index('R')
            if R > 0 and R < len(point):
                line_name = [point[R-1]+point[R], point[R]+point[R+1]]
            for m in range(len(line_name)):
                E0 = energy[x==value[i]]
                if m == 0:
                    bounds = [value[R-1], value[R]]
                if m == 1:
                    bounds = [value[R], value[R+1]]
                k_slice = np.bitwise_and(bounds[0]<=x, x<=bounds[1])
                k = x[k_slice]
                e = energy[k_slice]
                if i==0:
                    e = -e
                #  if np.min(e) == e[k==value[R]]:
                    #  k0 = 0
                    #  ER = 0
                    #  alpha = 0
                #  else:
                    #cs = CubicSpline(k, e)
                sp = UnivariateSpline(k, e, k=4)
                sp1 = sp.derivative()
                print(band[i], spin[j], line_name[m], sp1.roots(), sp1.roots().shape)
                spdr = sp1.roots()
                k0 = []
                if len(spdr) == 0:
                    k0 = 0
                else:
                    spdr = spdr[np.abs(spdr-value[R])<0.1]
                    print(spdr)
                    if len(spdr) == 0:
                        k0 = 0
                    elif len(spdr) > 1:
                        print('error')
                        exit()
                    else:
                        spdr = spdr[0]
                if k0 == 0:
                    ER = 0
                    alpha = 0
                else:
                    print(spdr, value[R], lattice)
                    k0 = np.abs((spdr-value[R])*2*np.pi/(lattice*b2a))
                    ER = np.abs(sp(spdr)-sp(value[R]))
                    #results = fmin(cs, value[R])
                    #k0 = (results[0]-value[R])*2*np.pi/(lattice*b2a)
                    #ER = cs(results)[0]-cs(value[R])
                    #ER = e[k==value[R]]-np.min(e)
                    #ER = ER[0]
                    #k0 = np.abs(value[R]-k[np.argmin(e)])
                    #k0 = k0*2*np.pi/(lattice*b2a)
                    alpha = 2*ER/k0
                rashba.append([lattice, band[i], spin[j], line_name[m], k0, ER, alpha])
    return rashba

def plot_rashba(rashba_lattice):
    m, n, o = np.shape(rashba_lattice)
    print(m, n ,o)
    lattice = []
    for i in range(m):
        lattice.append(rashba_lattice[i][0][0])
    lattice = np.array(lattice)*b2a
    print(lattice)
    label_name = []
    prefix = ['k0', 'ER', 'alpha']
    for i in range(n):
        label_name = ''.join(rashba_lattice[0][i][1:4])
        print(label_name)
        for j in range(3):
            print(prefix[j])
            y = []
            for k in range(m):
                y.append(rashba_lattice[k][i][j+4])
            plt.figure()
            plt.plot(lattice, y)
            if j == 0:
                plt.xlabel('Lattice Constant ($\mathrm{\AA}$)')
                plt.ylabel('Momentum offset (${\mathrm{\AA}}^{-1}$)')
            elif j == 1:
                plt.xlabel('Lattice Constant ($\mathrm{\AA}$)')
                plt.ylabel('Rashba Energy (eV)')
            else:
                plt.xlabel('Lattice Constant ($\mathrm{\AA}$)')
                plt.ylabel('Rashba coefficient (eV$\mathrm{\AA}$)')
            plt.tight_layout()
            plt.savefig(prefix[j]+''+label_name+'.pdf')
            plt.close()

def plot_convergence():
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
    plt.xlabel('ecutwfc (Ry)')
    plt.ylabel('Total Energy (Ry)')
    plt.tight_layout()
    plt.savefig('convergence.pdf')
    for i in range(len(cutoff)):
        print(cutoff[i], energy[i])
    np.savetxt('cutoff.txt',(cutoff, energy))




parser = argparse.ArgumentParser()
parser.add_argument("out", choices=['step','convergence', 'lattice', 'band', 'rashba'], help="The kind of file want to plot")
args = parser.parse_args()

if args.out == 'lattice':
    my_folder = filter_folder()
    energy = list(map(find_energy, my_folder))
    lattice = list(map(float, my_folder))
    plot_lattice_energy(lattice, energy)
    np.savetxt('lattice.txt',(lattice, energy))


if args.out == 'convergence':
    plot_convergence()


if args.out == 'band':
    my_folder = filter_folder()
    if my_folder:
        for i in my_folder:
            band(i)
    else:
        band()


if args.out == 'rashba':
    my_folder = filter_folder()
    if my_folder:
        rashba_lattice = []
        for i in my_folder:
            print('The lattice constant is: ',i)
            try:
                x, bm = find_bm_band(i)
            except:
                continue
            value, point = find_ticker(i)
            print(value, point)
            plot_bm_band(i, x, bm, value, point)
            rashba = find_rashba(float(i), x, bm, value, point)
            rashba_lattice.append(rashba)
        plot_rashba(rashba_lattice)
    else:
        x, bm = find_bm_band()
        value, point = find_ticker()
        print(value, point)
        plot_bm_band('.',x,bm,value,point)
        lattice = find_lattice_constant()
        print('The lattice constant is ', lattice*b2a)
        rashba = find_rashba(lattice, x, bm, value, point)
        print(rashba)

