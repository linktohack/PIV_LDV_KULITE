#!/usr/bin/env python2
import os
import numpy as np

def datarate(fn, rn=None, debug=False):
    """Calculate datarate of a file fn or a set of file fn with range rn"""
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
    
    if rn is None:
        if debug:
            print fn
        return _dr(fn)
    else:
        pos = []
        dr = []

        for i in xrange(rn[0], rn[1]):
            fn_ = fn % i
            if os.path.exists(fn_):
                if debug:
                    print fn_
                pos_, dr_ = _dr(fn_)
                pos.append(pos_)
                dr.append(dr_)
            else:
                if debug:
                    print '-', fn_
                
        return pos, dr


def quantities(fi, rn=None, fo=None, rot=0, off=(), debug=True):
    """ Calculate statistical quantities of a file name fi or a set of file with range rn
    In that case, an optional ouput file can be set to export the data"""
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
                
                U = (u**2 + v**2)**0.5
                
                k = (2*(u*u).mean()**0.5 + (v*v).mean()**0.5)/2
                uvbar = np.abs((u*v).mean())**0.5
    
                return pos, k, uvbar
    
    if rn is None:
        return qt(fi, rot, off, debug)
    else:
        pos = []
        k = []
        uvbar = []

        for i in xrange(rn[0], rn[1]):
            fn_ = fi % i
            if os.path.exists(fn_):
                if debug:
                    print fn_
                pos_, k_, uvbar_ = _qt(fn_, rot, off, debug)
                pos.append(pos_)
                k.append(k_)
                uvbar.append(uvbar_)
            else:
                if debug:
                    print '-', fn_
        
        if fo is not None:
            po = fo[fo.rfind('/')]
            if not os.path.exists(po):
                os.makedirs(po)
            np.save(fo, [pos, k, uvbar])

        return pos, k, uvbar
                

