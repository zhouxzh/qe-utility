#!/usr/bin/env python3

import argparse
import re


parser = argparse.ArgumentParser()
parser.add_argument("out", help="The qe input file")
parser.add_argument("-a", "--atom", help="The original atom you want to move")
parser.add_argument("-p", "--point", nargs=3, help="The point you want the original atom to move")
args = parser.parse_args()
if args.point:
    point = [float(i) for i in args.point]
else:
    point = [0, 0, 0]

atomic_position = []
with open(args.out) as f:
    for line in f:
        if 'ATOMIC_POSITION' in line:
            print(line, end='')
        m = re.search(r'^[A-Z][a-z]?\b',line)
        if m:
            atomic_position.append(line)


for i in atomic_position:
    if 'Pb' in i:
        origin = i.split()
        x0 = []
        for j in range(3):
            x0.append(float(origin[j+1]))

for i in atomic_position:
    atomic = i.split()
    for j in range(3):
        new_atomic = float(atomic[j+1])
        new_atomic = new_atomic - x0[j] + point[j]
        if new_atomic < 0:
            new_atomic = new_atomic + 1
        atomic[j+1] = '%10.6f'%(new_atomic)
    print(' '.join(atomic))

