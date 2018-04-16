#!/opt/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt

data=np.genfromtxt('energy.txt')
plt.plot(data[:,0], data[:,1])
plt.xlabel('ecutwfc (Ry) (ecutrho=8*ecutwfc)')
plt.ylabel('Total Energy (Ry)')
plt.tight_layout()
plt.savefig('convergence.pdf')
plt.show()
