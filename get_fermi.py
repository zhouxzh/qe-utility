#!/usr/bin/env python3

import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("out",help="The output file of qe")
args = parser.parse_args()
print(args.out)

with open(args.out) as f:
    for line in f:
        if 'highest' in line:
            print(line)
            if 'lowest' not in line:
                print("Error, no lowest unoccopied level")
                exit()
            m = re.findall('\d+\.\d+', line)
            fermi = (float(m[0])+float(m[1]))/2
            print('Fermi Energy is ', fermi, ' eV')


