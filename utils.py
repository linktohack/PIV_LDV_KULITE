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

def norm_s(sig, Pxx):
    """Normalize PSD so its integral is equal to the energy"""
    def norm_s_1d(sig, Pxx):
        m = np.mean(sig)
        xcor = np.mean((sig-m)**2)
        A = np.sum(Pxx)

        A = 2*A/xcor
        return Pxx/A

    if len(sig.shape) == 1:
        return norm_s_1d(sig, Pxx)
    else:
        for i in sig.shape[1]:
            Pxx[:,i] = norm_s_1d(sig[:,i], Pxx[i])
        return Pxx

# vim:set sw=4 ts=4 tw=78:
