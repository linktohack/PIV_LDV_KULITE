#!/usr/bin/env python2

from lib import waveletrec

def correlation(a, v):
    """Correlation between a & v"""
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) / np.std(v)

    return np.correlate(a, v, 'same')

def argpeaks(s, type='+', deviation=3, lookahead=10, cb=None):
    """Find the peaks of the signal within `deviation`*sigma range,
    exclude the fake signal within `lookahead` range.

    Parameter `type` can be `+` (default, ) `-`, `+-` for positive, negative
    and both peaks.

    Simple version.
    """
    ds = s[1:] - s[:-1]
    sts = []

    for t in type:
        if t == '+':
            st = np.where(ds > deviation*s.std())[0]
            if cb is not None:
                cb('First try, find positive peaks', st)

            for i in xrange(st.shape[0]-1, 0, -1):
                if (st[i] - st[i-1]) < lookahead:
                    if cb is not None:
                        cb('Fake peak detected', [st[i], st[i-1]])
                    st = np.delete(st, i)
            sts.append(st)

        if t == '-':
            st = np.where(ds < -deviation*s.std())[0]
            if cb is not None:
                cb('First try, find negative peaks', st)

            for i in xrange(st.shape[0]-1, 0, -1):
                if (st[i] - st[i-1]) < lookahead:
                    if cb is not None:
                        cb('Fake peak detected', [st[i], st[i-1]])
                    st = np.delete(st, i)
            sts.append(st)
    return sts

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

def wavelet_rec(t, sig, fs=1e5, NFFT=4096, stw=0.3, err=1e-3):
    """Reconstruct signal by Wavelet Technique
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

def linear_rec(Fs, t, sig):
    """Reconstruct signal by Linear Interpolation"""
    t_rec = np.arange(np.floor(t[-1]*Fs))/Fs
    sig_rec = np.zeros(t_rec.shape)
    
    first = np.where(t_rec<t[1])[0][-1]
    end = np.where(t_rec>=t[-2])[0][0]

    sig_rec[:first+1] = interpolate((t[0], sig[0]), (t[1], sig[1]), t_rec[:first+1])    
    sig_rec[end:] = interpolate((t[-2], sig[-2]), (t[-1], sig[-1]), t_rec[end:])
    
    j = 1
    last = first
    X = []
    for i in xrange(first, end):
        if t_rec[i] >= t[j]:
            sig_rec[last:i] = interpolate((t[j-1], sig[j-1]), (t[j], sig[j]), X)
            # If the samples too close
            while t_rec[i] >= t[j]:
                j += 1
            last = i
            X = []
        X.append(t_rec[i])
        
    return t_rec, sig_rec


# vim:set sw=4 ts=4 tw=78:
