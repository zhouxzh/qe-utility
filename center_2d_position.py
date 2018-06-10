#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='make the atomic position at the center position')
parser.add_argument('name', help='The input qe file')
args = parser.parse_args()
print(args.name)

position=[]
with open(args.name) as f:
    for i in f:
        if i.find('ATOMIC_POSITIONS')>=0:
            j=f.readline()
            print(j)

