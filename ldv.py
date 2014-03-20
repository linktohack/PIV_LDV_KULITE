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

def ldv_wavelet_rec(t, sig, fs=1e5, NFFT=4096, stw=0.3, err=1e-3):
    """Reconstruct LDV signal by Wavelet Technique
    Module used by JAUNET"""
    nsig1 = sig.shape[1]
    n1rec = NFFT
    nmin = 2
    ferec = np.float64(fs)
    erconv = np.float64(err)

    no_part = int(t[-1]/(n1rec-1)*ferec)

    s1rec_all = np.empty([n1rec*no_part,nsig1], dtype=np.float64, order='F')
    at1rec_all = np.arange(n1rec*no_part)/ferec

    part = 0
    n0 = 0
    n1 = 0
    for i in xrange(len(t)):
        if t[i] > at1rec_all[n1rec*(part+1)-1]:
            n0 = n1
            n1 = i + 1

            at1rec = at1rec_all[n1rec*part:n1rec*(part+1)]
            at1rec = at1rec - at1rec_all[n1rec*part]
            s1rec = np.empty([n1rec,nsig1], dtype=np.float64, order='F')

            at1 = t[n0:n1] - at1rec_all[n1rec*part]
            s1 = sig[n0:n1,:]
            s1[:,0] = s1[:,0] - np.mean(s1[:,0])

            dt1 = at1[1:] - at1[:-1]
            mask = np.where(dt1 > 1/ferec)[0]
            nfil = len(mask)

            at2 = at1[mask]
            s2 = s1[mask,:]


            waveletrec.waveletrecldv(at2, s2, at1rec, s1rec, erconv, stw, nmin,
                    nn=nfil, nsig=nsig1, nn2=n1rec)


            s1rec_all[n1rec*part:n1rec*(part+1),:] = s1rec

            # Next part
            part += 1
            if n1rec*(part+1) > len(s1rec_all):
                break

    return at1rec_all, s1rec_all

# vim:set ts=4 sw=4 tw=78:
