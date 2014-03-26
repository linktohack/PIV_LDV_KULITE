#!/usr/bin/env python2

def correlation(a, v):
    """Correlation between a & v"""
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) / np.std(v)

    return np.correlate(a, v, 'same')

def argpeaks(s, deviation=3, lookahead=10, cb=None):
    """Find the peaks of the signal within `deviation`*sigma range,
    exclude the fake signal within `lookahead` range.

    Simple version.
    """
    ds = s[1:] - s[:-1]
    st = np.where(ds > deviation*s.std())[0]
    if cb:
        cb('First try, find the peaks', st)

    for i in xrange(st.shape[0]-1, 0, -1):
        if (st[i] - st[i-1]) < 10:
            if cb:
                cb('Fake peak in st: ', [st[i], st[i-1]])
            st = np.delete(st, i)
    return st

# vim:set sw=4 ts=4 tw=78:
