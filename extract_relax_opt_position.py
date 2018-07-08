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
                m = re.search(r'^[A-Z][a-z]?\b',line)
                if m:
                    print(line, end='')
                line=f.readline()
