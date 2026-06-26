#!/usr/bin/env python3
"""
Strategy 1 -- Pc5 -> seismic-hum coherence-bridge detector  (target 3.70 mHz)

Tests the falsifiable prediction:  after a Pc5 ULF burst, does the seismic hum
amplitude at 3.70 mHz rise and PERSIST (>5x burst duration) with its PHASE
LOCKED to the driver?  Foam hypothesis -> yes (high ratio AND high PLV).
Standard physics -> no (ratio ~1, PLV ~0).

USAGE
  Real data:   python foam_coherence_detector.py pc5.csv hum.csv
               each CSV: two columns  ->  time_seconds,value   (uniform 60 s cadence)
               pc5.csv  = ground magnetometer, Pc5 band (nT)
               hum.csv  = vertical seismometer / superconducting gravimeter
  No data:     python foam_coherence_detector.py        # runs synthetic validation
"""
import sys, numpy as np
from scipy.signal import butter, filtfilt, hilbert

DT, F_HUM, F_PC5 = 60.0, 3.70e-3, 3.0e-3
BURST_DUR = 1.5*3600

def narrow(x, f0, bw=0.6e-3):
    b, a = butter(4, [(f0-bw/2)/(0.5/DT), (f0+bw/2)/(0.5/DT)], btype='band')
    z = hilbert(filtfilt(b, a, x)); return np.abs(z), np.angle(z)

def find_bursts(B):
    amp, _ = narrow(B, F_PC5, bw=5.0e-3)
    thr = amp.mean() + 2*amp.std()
    hot = amp > thr; onset_idx = np.where((~hot[:-1]) & (hot[1:]))[0]
    return onset_idx * DT

def detect(t, B, S):
    amp, phase = narrow(S, F_HUM)
    _, drv_phase = narrow(B, F_PC5, bw=5.0e-3)
    onsets = find_bursts(B)
    pre, post, plv = [], [], []
    for o in onsets:
        m_pre  = (t >= o-2*3600) & (t < o)
        m_post = (t >= o+BURST_DUR+0.5*3600) & (t < o+BURST_DUR+12*3600)
        if m_pre.sum() > 5 and m_post.sum() > 5:
            pre.append(amp[m_pre].mean()); post.append(amp[m_post].mean())
            k = np.argmin(np.abs(t-(o+5*3600)))
            plv.append(np.exp(1j*(phase[k]-drv_phase[min(k, len(t)-1)])))
    if not pre: return None
    ratio = np.mean(post)/np.mean(pre); PLV = abs(np.mean(plv))
    return ratio, PLV, len(pre)

def verdict(ratio, PLV):
    if ratio > 1.3 and PLV > 0.4:  return "FOAM-LIKE SIGNAL  (persistent, phase-locked amplification)"
    if ratio < 1.1 and PLV < 0.2:  return "NULL  (consistent with standard physics: no coupling)"
    return "AMBIGUOUS  (needs more events / cleaner data)"

def load(p):
    d = np.loadtxt(p, delimiter=',', skiprows=1)
    return d[:,0], d[:,1]

def synthetic():
    days = 45; N = days*1440; t = np.arange(N)*DT; rng = np.random.default_rng(7)
    onsets = np.sort(rng.uniform(0, t[-1]*0.9, days))
    B = 0.5*rng.standard_normal(N)
    for o in onsets:
        m = (t>=o)&(t<o+BURST_DUR)
        B[m] += 8*np.sin(np.pi*(t[m]-o)/BURST_DUR)*np.sin(2*np.pi*F_PC5*(t[m]-o))
    ocean = rng.standard_normal(N)
    S_foam = ocean.copy(); S_null = ocean.copy()
    for o in onsets:
        s=o+BURST_DUR; m=(t>=s)
        S_foam[m] += 1.2*np.exp(-(t[m]-s)/(10*3600))*np.sin(2*np.pi*F_HUM*(t[m]-s))
        mb=(t>=o)&(t<o+BURST_DUR); S_null[mb]+=0.05*np.sin(2*np.pi*F_HUM*(t[mb]-o))
    for name,S in [("FOAM world",S_foam),("NULL world",S_null)]:
        r,p,n = detect(t,B,S); print(f"  {name:12s}: ratio={r:.2f}  PLV={p:.2f}  ({n} bursts) -> {verdict(r,p)}")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        (tb,B),(ts,S) = load(sys.argv[1]), load(sys.argv[2])
        n=min(len(tb),len(ts)); t,B,S = tb[:n],B[:n],S[:n]
        r,p,k = detect(t,B,S)
        print(f"events={k}  post/pre amp={r:.2f}  PLV={p:.2f}\nVERDICT: {verdict(r,p)}")
    else:
        print("No data given -- synthetic validation:\n"); synthetic()
