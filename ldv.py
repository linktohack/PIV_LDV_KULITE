#!/usr/bin/env python2
import os
import numpy as np

from scipy.optimize import curve_fit
from scipy.special import erf

import glob

def datarate(fn, cb=None):
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
            return { 'pos': pos, 'dr': dr }

    fl = sorted(glob.glob(fn))
    if not fl:
        raise IOError('File not found! %s' % fn)
    if len(fl) < 2:
        if cb:
            cb('File name', fl[0])
        return _dr(fl[0])
    else:
        dr = {}
        for fn_ in fl:
            if cb:
                cb('File name', fn_)
            d = _dr(fn_)
            for k in d.iterkeys():
                dr.setdefault(k, []).append(d[k])

        return dr


def quantities(fn, rot=0, fit=[450, 50, 25, 10], cb=None):
    """ Calculate statistical quantities of file(s) with pattern fn

    Also try to fit the curve to the erf function with `fit=[ua, ub, yref, dw]`
    """
    def _qt(fn, rot, cb):
        """ Calculate statistical quantities of one file"""
        with open(fn, 'r') as f:
            header = [f.next() for _ in xrange(6)]
            pos = np.array(header[3].strip()[:-4].split(' mm;')).astype(np.float64)
            if cb:
                cb('Position', pos)

            lda = np.loadtxt(fn, skiprows=6, usecols=(3,4), dtype=np.float64)
            if cb:
                cb('Number of samples', lda.shape[0])

            a1 = np.cos((1./4 - rot/180.)*np.pi)
            a2 = np.sin((1./4 - rot/180.)*np.pi)
            b1 = np.cos((-1./4 - rot/180.)*np.pi)
            b2 = np.sin((-1./4 - rot/180.)*np.pi)

            u = lda[:,0]*a1 + lda[:,1]*b1
            v = lda[:,0]*a2 + lda[:,1]*b2
            U = np.mean((u**2 + v**2)**0.5)

            # std-n (rms), oh fluid dynamic!
            su = u.std()
            sv = v.std()

            # u', v'
            u = u - u.mean()
            v = v - v.mean()


            k = (2*np.mean(u**2) + np.mean(v**2))/2

            uv = np.abs(np.mean(u*v))
            u_v = np.abs(np.mean(u/v))

            return {
                      'pos': pos,
                      'U': U,
                      'su': sv,
                      'sv': sv,
                      'k': k,
                      "u'v'": uv,
                      "u'/v'": u_v
                    }


    fl = sorted(glob.glob(fn))
    if not fl:
        raise IOError(2, 'No such file or directory' % fn)
    if len(fl) < 2:
        if cb:
            cb('File name', fl[0])
        qt = _qt(fl[0], rot, cb)
    else:
        qt = {}
        for fn_ in fl:
            if cb:
                cb('File name', fn_)
            q = _qt(fn_, rot, cb)
            for k in q.iterkeys():
                qt.setdefault(k, []).append(q[k])

        # fit the curve
        f = lambda x, ua, ub, y0, dw: 0.5*(1+erf((np.sqrt(np.pi)/dw)*(-x+y0)))*(ua-ub)+ub
        # f2 = lambda x, ub, y0, dw: f(x, ua, ub, y0, dw)
        z = np.array([p[2] for p in qt['pos']])
        U = np.array(qt['U'])
        popt, pcov = curve_fit(f, z, U, [450, 50, 25, 10])
        qt['fit'] = popt
    return qt

# vim:set ts=4 sw=4 tw=78:
