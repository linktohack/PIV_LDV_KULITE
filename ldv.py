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
            return dr, pos
    
    if rn is None:
        if debug:
            print fn
        return _dr(fn)
    else:
        dr = []
        pos = []
        for i in xrange(rn[0], rn[1]):
            fn_ = fn % i
            if os.path.exists(fn_):
                if debug:
                    print fn1
                dr_, pos_ = _dr(fn_)
                dr.append(dr_)
                pos.append(pos_)
            else:
                if debug:
                    print '-', fn1
                
        return dr, pos

