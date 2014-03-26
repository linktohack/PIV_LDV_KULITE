#!/usr/bin/env python2
import os
import numpy as np

import glob

from lib import waveletrec

def datarate(fn, debug=False):
    """Calculate datarate of file(s) with pattern fn"""
    def _dr(fn):
        """Calculate datarate of one file"""
        with open(fn, 'r') as f:
            header = [f.next() for _ in xrange(6)]
            pos = np.array(header[3].strip()[:-4].split(' mm;')).astype(np.float64)
            firstline =  f.next()

            for lastline in f:
                pass

            first = np.array(firstline.strip().split('\t')).astype(np.float64)
            last = np.array(lastline.strip().split('\t')).astype(np.float64)

            dr = (last[0] - first[0])/(last[1] - first[1])*1000
            return pos, dr

    fl = sorted(glob.glob(fn))
    if not fl:
        raise IOError('File not found! %s' % fn)
    if len(fl) < 2:
        if debug:
            print fl[0]
        return _dr(fl[0])
    else:
        pos = []
        dr = []

        for fn_ in fl:
            if debug:
                print fn_
            pos_, dr_ = _dr(fn_)
            pos.append(pos_)
            dr.append(dr_)

        return pos, dr


def quantities(fi, fo=None, rot=0, off=(), debug=False):
    """ Calculate statistical quantities of file(s) with pattern fi
    Optional ouput file can be set to export the data"""
    def _qt(fn, rot, off, debug):
        """ Calculate statistical quantities of one file"""
        with open(fn, 'r') as f:
            header = [f.next() for _ in xrange(6)]
            pos = np.array(header[3].strip()[:-4].split(' mm;')).astype(np.float64)
            if debug:
                print '[x y z] = ', pos

            lda = np.loadtxt(fn, skiprows=6, usecols=(3,4), dtype=np.float64)
            if debug:
                print 'Number of samples = ', lda.shape[0]

            a1 = np.cos((1./4 - rot/180.)*np.pi)
            a2 = np.sin((1./4 - rot/180.)*np.pi)
            b1 = np.cos((-1./4 - rot/180.)*np.pi)
            b2 = np.sin((-1./4 - rot/180.)*np.pi)

            u = lda[:,0]*a1 + lda[:,1]*b1
            v = lda[:,0]*a2 + lda[:,1]*b2
            U = np.mean((u**2 + v**2)**0.5)

            u = u - u.mean() # ubar, vbar
            v = v - v.mean()

            k = (2*(u*u).mean()**0.5 + (v*v).mean()**0.5)/2
            uvbar = np.abs((u*v).mean())**0.5

            return pos, U, k, uvbar

    fl = sorted(glob.glob(fi))
    if not fl:
        raise IOError(2, 'No such file or directory' % fi)
    if len(fl) < 2:
        if debug:
            print fl[0]
        pos, U, k, uvbar = _qt(fl[0], rot, off, debug)
    else:
        pos = []
        U = []
        k = []
        uvbar = []

        for fn_ in fl:
            if debug:
                print fn_
            pos_, U_, k_, uvbar_ = _qt(fn_, rot, off, debug)
            pos.append(pos_)
            U.append(U_)
            k.append(k_)
            uvbar.append(uvbar_)

    if fo is not None:
        po = fo[:fo.rfind('/')]
        if not os.path.exists(po):
            os.makedirs(po)
        np.save(fo, [pos, U, k, uvbar])

    return pos, U, k, uvbar

# vim:set ts=4 sw=4 tw=78:
