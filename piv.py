#!/usr/bin/env python2
import numpy as np

import glob

def piv_to_im(piv):
    im1 = np.empty([256, 256], dtype=np.float32)
    im2 = np.empty([256, 256], dtype=np.float32)

    for i in xrange(256):
        im1[i,:] = piv[256*i:256*(i+1),2]
        im2[i,:] = piv[256*i:256*(i+1),3]

    return im1, im2

def piv_mean_rms(fn, debug=False):
    fp = sorted(glob.glob(fn))

    piv1 = np.empty([len(fp), 65536])
    piv2 = np.empty([len(fp), 65536])

    for i in xrange(len(fp)):
        if debug:
            if (i-1) % 10 == 0:
                print fp[i]
        piv1[i,:], piv2[i,:] = np.loadtxt(fp[i], dtype=np.float32, usecols=(2,3), unpack=True)

    pivm1 = piv1.mean(axis=0)
    pivm2 = piv2.mean(axis=0)

    for i in xrange(len(fp)):
        piv1[i,:] = piv1[i,:] - pivm1
        piv2[i,:] = piv2[i,:] - pivm2

    pivrms1 = np.sqrt(np.mean(piv1**2, axis=0))
    pivrms2 = np.sqrt(np.mean(piv2**2, axis=0))

    pos = np.loadtxt(fp[0], dtype=np.float32, usecols=(0,1))
    repivm = np.concatenate([pos, pivm1.reshape((65536,1)), pivm2.reshape((65536,1))], axis=1)
    repivrms = np.concatenate([pos, pivrms1.reshape((65536,1)), pivrms2.reshape((65536,1))], axis=1)

    return repivm, repivrms

# vim:set ts=4 sw=4:
