#!/usr/bin/env python2

def correlation(a, v):
    """Correlation between a & v"""
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) / np.std(v)

    return np.correlate(a, v, 'same')

# vim:set sw=4 ts=4 tw=78:
