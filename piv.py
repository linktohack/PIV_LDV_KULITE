#!/usr/bin/env python2
import numpy as np

import glob

import utils

def piv_to_im(piv):
    im1 = np.empty([256, 256], dtype=np.float32)
    im2 = np.empty([256, 256], dtype=np.float32)

    for i in xrange(256):
        # Last 2 colums are u & v, postitions are optional
        im1[i,:] = piv[256*i:256*(i+1),-2]
        im2[i,:] = piv[256*i:256*(i+1),-1]

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

def piv_filter_n_sigma(piv, pm, pr, n=3, debug=False):
    """Filter the signal by n(three) sigma"""
    errors_count = np.zeros(piv.shape[0], dtype=np.int32)

    for j in xrange(piv.shape[1]):
        if debug:
            if j%1000 == 0:
                print 'loop %d/65536, errors count %d' % (j, errors_count.sum())
        save_end = ()
        ends = []
        middle_X = []
        middle_Y = []

        for i in xrange(piv.shape[0]):
            finished = False
            if i == piv.shape[0] - 1:
                finish = True

            X = (i, piv[i,j])
            # bad point
            if np.abs(piv[i,j] - pm[j]) > 3*np.abs(pr[j]):
                middle_X.append(X)
            # good point, keep only one good point, save the other end
            # in case of there is no other good point at the end
            else:
                ends.append(X)
                if len(middle_X) == 0 and len(ends) > 1:
                    save_end = ends[-2]
                    ends = [ends[-1]]

            # restore the saved end
            if finished:
                ends = save_end + ends

            # interpolate when there are 2 good points and a list of bad points
            # TODO: Raise exeption when there is only one good point. NO WAY!
            if len(middle_X) > 0 and len(ends) == 2:
                middle_Y = interpolate(ends[-2], ends[-1], [x[0] for x in middle_X])
                errors_count[i] = errors_count[i] + len(middle_X)
                for (x, _), y in zip(middle_X, middle_Y):
                    piv[x,j] = y
                middle_X = []

    return piv, errors_count

# vim:set ts=5 sw=4 tw=78:
