#!/usr/bin/env python2
import os
import numpy as np

def cali_coeffs(fpgen, fp, fpamb, rn=range(1, 11)):
    pgen = np.loadtxt(fpgen)

    p = []
    for i in xrange(10):
        pamb = np.fromfile(fpamb % (i+1, i+1), dtype=np.float32)
        pi = [pamb.mean()]
        for j in [800, 900, 1000, 1100, 1200]:
            pj = np.fromfile(fp % (i+1, j, i+1), dtype=np.float32)
            pi.append(pj.mean())
        p.append(pi)
    p = np.array(p).T

    fit = np.empty(p.shape[1])
    for i in xrange(len(fit)):
        fit[i] = np.polyfit(p[:,i], pgen[:,i], 1)[0]

    return fit

# vim:set ts=4 sw=4 et:
