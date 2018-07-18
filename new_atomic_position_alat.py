#!/usr/bin/env python3

import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("out",help="The vc-relax output file of qe")
args = parser.parse_args()
print(args.out)

with open(args.out) as f:
    for line in f:
        if 'Begin final' in line:
            #print(line)
            while 'End final' not in line:
                if 'CELL' in line:
                    m = re.search('[0-9]+\.[0-9]+',line)
                    alat = float(m.group())
                    cell = np.zeros((3,3))
                    for i in range(3):
                        a = f.readline().split()
                        for j in range(3):
                            cell[i,j]=float(a[j])
                    print('new alat is ', alat*cell[0,0])
                    print('new cell is \n', cell/cell[0,0])
                    new_cell=cell[0,0]
                m = re.search(r'^[A-Z][a-z]?\b',line)
                if m:
                    atomic=line.split()
                    for i in range(3):
                        new_atomic = float(atomic[i+1])/new_cell
                        if new_atomic < 0:
                            new_atomic = new_atomic + 1
                        if new_atomic > 1:
                            new_atomic = new_atomic - 1
                        atomic[i+1]='%10.6f'%(new_atomic)
                    print(' '.join(atomic))
                line=f.readline()
