#!/usr/bin/env python2
import os
import numpy as np

def datarate(fn, rn=None, debug=False):
    """Calculate datarate of a file fn or a set of file fn with range rn"""
    def dr1(fn):
        with open(fn, 'r') as f:
            [f.next() for _ in xrange(6)]    
            firstline =  f.next()
        
            for lastline in f:
                pass

            first = np.array(firstline.strip().split('\t')).astype(np.float64)        
            last = np.array(lastline.strip().split('\t')).astype(np.float64)

            dr = (last[0] - first[0])/(last[1] - first[1])*1000   
            return dr
    
    if rn is None:
        if debug:
            print fn
        return dr1(fn)
    else:
        dr = []
        for i in xrange(rn[0], rn[1]):
            fn1 = fn % i
            if os.path.exists(fn1):
                if debug:
                    print fn1
                dr.append(dr1(fn1))
            else:
                if debug:
                    print '-', fn1
                
        return dr
