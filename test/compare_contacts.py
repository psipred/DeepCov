#!/usr/bin/env python

import sys
from numpy import argsort

def get_ordering(filename=""):
    scores = list()
    f = open(filename, 'r')
    f = f.readlines()
    for line in f:
        l = line.split()
        scores.append(l[4])
        
    return argsort(scores)

def main():
    if len(sys.argv) < 3:
        print('usage: python compare_contacts.py file1 file2')
        sys.exit(2)

    o1 = get_ordering(sys.argv[1])
    o2 = get_ordering(sys.argv[2])

    ## Testing: manually create a mismatch by swapping two entries in one of the lists
#    o1[1], o1[2] = o1[2], o1[1]

    if all(o1 == o2):
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
