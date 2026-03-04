"""This file contains the analytic feedback model we developed for the paper (Wagg+2025).

Example usage:

>>> from analytic_feedback_model import sample_ccsn
>>> times, vels = sample_ccsn(n=1000, FeH=-0.5)"""

import numpy as np
from scipy.stats import beta, maxwell, norm


def get_ccsn_rate_parameters(FeH=0, a1=0.38, a2=0.47, a3=0.22,
                             a4=0.13, a5=0.1, a6=0.175,
                             t1=3.5, t2=6, t3=23, t4=28, t5=45.5, t6=200):
    """Get the metallicity-dependent parameters for the CCSN rate.

    Parameters
    ----------
    FeH : `float`, optional
        Metallicity, by default 0
    a_i : `float`, optional
        Solar metallicity values for a parameters
    t_i : `float`, optional
        Solar metallicity values for transition times
    """
    # add metallicity dependence
    a1 += 0.13 * FeH
    a2 += 0.05 * FeH
    a3 += 0.02 * FeH
    a6 += 0.05 * FeH
    t3 -= 6.5 * FeH
    t4 -= 6.5 * FeH
    t5 -= 16.5 * FeH

    a = [a1, a2, a3, a4, a5, a6]
    ts = [t1, t2, t3, t4, t5, t6]

    # calculate the psi values
    psi = [np.log(a[i + 1] / a[i]) / np.log(ts[i + 1] / ts[i]) for i in range(len(a) - 1)]

    return a, ts, psi

def cdf_piece(t, i, a, ts, psi):
    return a[i] * ts[i]**-psi[i] * (t**(psi[i] + 1) - ts[i]**(psi[i] + 1)) / (psi[i] + 1)

def cumulative_rate(t, FeH=0):
    a, ts, psi = get_ccsn_rate_parameters(FeH=FeH)

    cumul_below_t = np.array([0.0 for t in ts])

    for i in range(1, len(ts) - 1):
        cumul_below_t[i:] += cdf_piece(ts[i], i - 1, a, ts, psi)

    @np.vectorize
    def get_cdf(t):
        cdf = None
        if t < ts[0]:
            cdf = 0.0
        elif t < ts[1]:
            cdf = cdf_piece(t, 0, a, ts, psi) + cumul_below_t[0]
        elif t < ts[2]:
            cdf = cdf_piece(t, 1, a, ts, psi) + cumul_below_t[1]
        elif t < ts[3]:
            cdf = cdf_piece(t, 2, a, ts, psi) + cumul_below_t[2]
        elif t < ts[4]:
            cdf = cdf_piece(t, 3, a, ts, psi) + cumul_below_t[3]
        elif t <= ts[5]:
            cdf = ts[4] * a[5] * (np.exp(-1) - np.exp(-t / ts[4])) + cumul_below_t[4]
    
        return cdf

    return get_cdf(t)

def inverse_cdf_rate(u, FeH=0, norm=None, transitions=None):
    a, ts, psi = get_ccsn_rate_parameters(FeH)
    
    norm = 1 / cumulative_rate(200, FeH=FeH)
    transitions = [cumulative_rate(t, FeH=FeH) * norm for t in ts]
    
    @np.vectorize
    def get_val(u):    
        if u < transitions[1]: # 
            return (ts[0]**psi[0] * ((u / norm - transitions[0] / norm) * (psi[0] + 1) / a[0] + ts[0]))**(1 / (psi[0] + 1))
        elif u < transitions[2]:
            return (ts[1]**psi[1] * ((u / norm - transitions[1] / norm) * (psi[1] + 1) / a[1] + ts[1]))**(1 / (psi[1] + 1))
        elif u < transitions[3]:
            return (ts[2]**psi[2] * ((u / norm - transitions[2] / norm) * (psi[2] + 1) / a[2] + ts[2]))**(1 / (psi[2] + 1))
        elif u < transitions[4]:
            return (ts[3]**psi[3] * ((u / norm - transitions[3] / norm) * (psi[3] + 1) / a[3] + ts[3]))**(1 / (psi[3] + 1))
        else:
            return -ts[4] * np.log(1 / np.e - (u / norm - transitions[4] / norm) / (a[5] * ts[4]))

    return get_val(u)

def sample_ccsn_times(n, FeH=0):
    return inverse_cdf_rate(np.random.rand(n), FeH=FeH)

def sample_ccsn(n, FeH=0, v_disp=1.7):
    """Sample CCSN times and velocities.
    This function samples the CCSN times and then the progenitor velocities conditioned on the times.
    Full model is described in Wagg+2025

    Parameters
    ----------
    n : `int`
        Number of CCSN to sample
    FeH : `float`, optional
        Metallicity, by default 0
    v_disp : `float`, optional
        Velocity dispersion in km/s, by default 1.7 km/s

    Returns
    -------
    times : `ndarray`
        CCSN times
    vels : `ndarray`
        CCSN velocities
    """
    times = sample_ccsn_times(n=n, FeH=FeH)
    vels = maxwell(scale=v_disp / np.sqrt(3)).rvs(n)

    inds = np.arange(n)
    
    low_t_inds = inds[(times > 5) & (times < 45.5 - 16.5 * FeH)]
    eject_inds = np.random.choice(low_t_inds, int(len(low_t_inds) * 0.24), replace=False)
    
    low_t_inds = inds[(times > 45.5 - 16.5 * FeH) & (times < 60)]
    eject_inds = np.concatenate((eject_inds, np.random.choice(low_t_inds, int(len(low_t_inds) * 0.1), replace=False)))

    n_eject = len(eject_inds)
    f_no_MT = 0.14 - 0.12 * FeH
    n_no_MT = int(f_no_MT * n_eject)
    
    factor = -1.8 + np.sqrt(abs(FeH)) * 0.5
    f1 = factor + 1
    A = ((100**f1 - 5**f1) / f1)**(-1)
    u = np.random.rand(n_no_MT)
    eject_nomt = (u * f1 / A + 5**f1)**(1/f1)
    
    f_caseA = 0.12 + 0.035 * FeH
    n_caseA = int(f_caseA * n_eject)
    eject_caseA_mt = norm(loc=22 - 8 * FeH, scale=6 - 3 * FeH).rvs(n_caseA)

    f_caseB = 0.67 + 0.12 * FeH
    n_caseB = int(f_caseB * n_eject)
    eject_caseB_later_mt = beta(1.5 - 1.5 * FeH, 18 + 5 * FeH, loc=5, scale=100 - 5).rvs(n_caseB)

    n_CE = n_eject - n_no_MT - n_caseA - n_caseB
    eject_CE_mt = beta(5 - FeH * 4, 10, loc=5,scale=100 - 5).rvs(n_CE)

    eject_vels = np.concatenate((eject_caseA_mt, eject_CE_mt, eject_caseB_later_mt, eject_nomt))
    vels[eject_inds] = eject_vels

    return times, vels
