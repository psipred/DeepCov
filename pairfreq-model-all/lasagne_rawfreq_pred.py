#!/usr/bin/env python

# This script runs a trained DeepCov model using pair frequency data

from __future__ import print_function

import sys
import os
import time

from math import sqrt

import numpy as np
import theano
import theano.tensor as T

import lasagne

# This assumes nndef.py is in pwd
from nndef import build_cnn

def main():

    input_var = T.ftensor4('inputs')

    # Create neural network model
    network = build_cnn(input_var)
    
    # Load parameters; the file must be in pwd
    with np.load('FINAL_fullmap_metapsicov_model.npz') as f:
        param_values = [f['arr_%d' % i] for i in range(len(f.files))]
    
    lasagne.layers.set_all_param_values(network, param_values)

    # Load input data (.21m file)
    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)

    length = int(sqrt(mapdata.shape[0]/21/21))
    inputs = mapdata.reshape(1,441,length,length)

    predict = lasagne.layers.get_output(network, input_var, deterministic=True)
    predict_fn = theano.function([input_var], predict)

    # Make the predictions
    result = predict_fn(inputs)

    # Write output to stdout; average values for residue pairs i,j and j,i
    for wi in range(0, length-1):
        for wj in range(wi+1, length):
            print("{} {} 0 8 {}".format(wi+1, wj+1, 0.5 * (result[0,0,wi,wj] + result[0,0,wj,wi])))

if __name__=="__main__":
    main()
