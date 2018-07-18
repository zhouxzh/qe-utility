#!/usr/bin/env python3

import argparse
import re
import os
import os.path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("out", help="The file you want to copy, such as scf file, relax file")
args = parser.parse_args()
print('The ', args.out, ' results will be copy.')

results_dir=args.out+'_results'
if os.path.isdir(results_dir):
    print('The results dir have been made: '+results_dir)
else:
    os.mkdir(results_dir)
    print('The results dir will be made: '+results_dir)


for my_dir in os.listdir():
    if os.path.isdir(my_dir):
        if re.search('\d+(\.\d*)?', my_dir):
            print(my_dir)
            dst_dir = os.path.join(results_dir, my_dir)
            if not os.path.isdir(dst_dir):
                os.mkdir(dst_dir)
                print('destination dir is making: ',dst_dir)
            for my_file in os.listdir(my_dir):
                src_file = os.path.join(my_dir, my_file)
                if args.out in my_file:
                    print(my_file)
                    dst_file = os.path.join(results_dir, my_dir, my_file)
                    shutil.copy(src_file, dst_file)




