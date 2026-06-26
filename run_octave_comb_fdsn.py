# -*- coding: utf-8 -*-
"""
Fetch REAL open broadband seismic data from any FDSN data center (no account
needed) and run the locked octave-comb detector on it.

This is the real-data driver for Arm 3. It needs internet + obspy, which this
sandbox does not have, so run it on your own machine:

    pip install obspy scipy numpy matplotlib
    python run_octave_comb_fdsn.py

It defaults to a very quiet Global Seismographic Network station (II.BFO, the
Black Forest Observatory) and one month of 1 Hz vertical data (LHZ). Every data
center below is free and open; change CLIENT/NET/STA/CHA to point anywhere.

  FDSN center      obspy short name        region
  ------------     -----------------       ------------------------
  EarthScope/IRIS  "IRIS"                  global (GSN: IU, II, IC) + USA
  GEOFON (GFZ)     "GFZ"                   global + Germany/Europe
  ORFEUS / EIDA    "ORFEUS"                Europe (federated)
  RESIF            "RESIF"                 France
  ETH / SED        "ETH"                   Switzerland
  INGV             "INGV"                  Italy
  GeoNet           "GEONET"                New Zealand
  AusPass          "AUSPASS"               Australia
  BGR              "BGR"                   Germany

Channel codes by sampling rate (pick by which comb teeth you need):
  VHZ = 0.1 Hz  (Nyquist 50 mHz)  -> teeth 0-4 (2..32 mHz); smallest files
  LHZ = 1   Hz  (Nyquist 0.5 Hz)  -> all 8 demo teeth (2..257 mHz)  [default]
  BHZ = 20-40Hz                   -> up to ~1 Hz octaves; large files
"""
import numpy as np
from scipy import signal

# ----------------------------- comb detector (general fs) -------------------
F0 = 2.005e-3
ALPHA = 0.05
TEMPLATES = {
    "flat (incoherent)": lambda m: np.ones(m),
    "dying 1/2^k":       lambda m: 0.5 ** np.arange(m),
    "dying 1/k":         lambda m: 1.0 / np.arange(1, m + 1),
    "climbing 2^k":      lambda m: 2.0 ** np.arange(m),
}
def _normw(w): return w / w.sum() * len(w)

def comb_bins(f0, fs, nper):
    df = fs / nper; bins, k = [], 0
    while True:
        b = int(round(f0 * 2 ** k / df))
        if b >= nper // 2: break
        bins.append(b); k += 1
    return np.array(bins)

def integer_ctrl(f0, fs, nper, mults=(3, 5, 6, 7)):
    df = fs / nper
    return np.array([int(round(f0 * m / df)) for m in mults
                     if int(round(f0 * m / df)) < nper // 2])

def whiten(P, protect, halfwin=60, guard=3):
    n = len(P); keep = np.ones(n, bool)
    for b in protect: keep[max(0, b-guard):min(n, b+guard+1)] = False
    B = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i-halfwin), min(n, i+halfwin+1)
        seg = P[lo:hi][keep[lo:hi]]
        B[i] = np.median(seg) if seg.size >= 5 else np.median(P[lo:hi])
    return P / B

def detect(P, fs, nper, f0=F0, nboot=200000):
    octv = comb_bins(f0, fs, nper); integ = integer_ctrl(f0, fs, nper)
    m = len(octv)
    if m < 2:
        raise SystemExit(f"Only {m} comb tooth below Nyquist ({fs/2*1e3:.2f} mHz). "
                         "Use a faster channel (LHZ or BHZ).")
    rw = whiten(P, np.concatenate([octv, integ]))
    lo = 60; pool = rw[lo:len(rw)-60]
    tau = np.quantile(pool, 1 - ALPHA)
    r = rw[octv]; S = float(r.sum()); C = int((r > tau).sum())
    g = np.random.default_rng(7); draws = g.choice(pool, size=(nboot, m), replace=True)
    p_sum = float(np.mean(draws.sum(1) >= S))
    p_cnt = float(np.mean((draws > tau).sum(1) >= C))
    k = np.arange(m); up = r > 1.0
    slope = np.polyfit(k[up], np.log(r[up]), 1)[0] if up.sum() >= 2 else np.nan
    print(f"\nlocked octave teeth ({m}): {[round(b*fs/nper*1e3,3) for b in octv]} mHz")
    for b, v in zip(octv, r):
        print(f"   {b*fs/nper*1e3:9.3f} mHz : {v:7.2f}{'  *' if v>tau else ''}")
    print(f"\nincoherent SUM  S={S:.2f}  p_sum={p_sum:.4f}")
    print(f"coincidence CNT C={C}      p_cnt={p_cnt:.4f}")
    det = (p_sum < ALPHA) and (p_cnt < ALPHA)
    print(f"octave teeth(2f..) mean {float(r[1:].mean()):.2f} vs "
          f"non-octave integer mean {float(rw[integ].mean()):.2f}")
    if det:
        print(f"\n>>> COMB DETECTED (both gates) | data octave slope {slope:+.2f}/octave")
        print(">>> NOTE: a detection is an UNEXPLAINED line, not proof of foam. Exclude "
              "instrumental, tidal-harmonic, and atmospheric-resonance origins first.")
    else:
        amin = 2*np.sqrt(np.log(1/ALPHA)/ (P.size))  # rough sensitivity scale
        print(f"\n>>> NO COMB (consistent with null) -> this is the UPPER BOUND result.")
    return det

# ----------------------------- fetch + run ----------------------------------
def main():
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    # ---- archive menu (all FDSN; pick one; the detector is identical for every one) ----
    # ObsPy built-in shortcuts -- no base URL needed, just set CLIENT to the name:
    #   "IRIS"   global GSN (IU/II/IC) + USA      "GFZ"     GEOFON: global + Germany
    #   "ORFEUS" Europe (EIDA federation)         "RESIF"   France
    #   "ETH"    Switzerland (SED)                "INGV"    Italy
    #   "NOA"    Greece                           "BGR"     Germany (federal)
    #   "KOERI"  Turkey (Bogazici)                "GEONET"  New Zealand
    #   "AUSPASS" Australia                       "SCEDC"   S. California (dense; urban=noisy)
    #   "NCEDC"  N. California                    "USGS"    USGS
    # Nodes without a shortcut (AFAD-Turkey, CSN-Chile, SSN-Mexico): pass the FDSN base
    # URL instead -> Client(base_url="https://<node>")  (endpoints: https://www.fdsn.org/datacenters/)
    CLIENT = "IRIS"
    # Quiet, low-noise, long-period stations (best for a 2 mHz comb), by archive:
    #   IRIS : ("II","BFO","00")  ("IU","TSUM","00")  ("IU","ANMO","00")  ("II","NNA","00")
    #   ORFEUS/GFZ : ("II","BFO","00") / ("GE","STU","")    GEONET : ("NZ","SNZO","10")
    NET, STA, LOC, CHA = "II", "BFO", "00", "VHZ"   # VHZ=0.1Hz -> ~1 MB/month, 5 octave teeth
    #                                                (use "LHZ"=1Hz for all 8 teeth, still tiny)
    T1=UTCDateTime("2023-06-01T00:00:00"); T2=T1 + 30*86400   # 30 days
    print(f"fetching {NET}.{STA}.{LOC}.{CHA} {T1.date}..{T2.date} from {CLIENT} ...")
    c=Client(CLIENT)
    st=c.get_waveforms(NET,STA,LOC,CHA,T1,T2,attach_response=True)
    st.merge(method=1, fill_value="interpolate")          # bridge small gaps
    tr=st[0]; fs=tr.stats.sampling_rate
    print(f"got {tr.stats.npts} samples at {fs} Hz "
          f"(Nyquist {fs/2*1e3:.1f} mHz, span {(T2-T1)/86400:.0f} d)")
    try:
        tr.remove_response(output="VEL")                  # optional; whitening also handles it
    except Exception as e:
        print("response removal skipped:", e)
    tr.detrend("linear"); tr.detrend("demean"); tr.taper(0.01)
    x=tr.data.astype(float)
    nper=1<<int(np.log2(fs*60000))                        # ~16-18 h segments
    nper=min(nper, len(x)//4)
    f,P=signal.welch(x, fs=fs, nperseg=nper, noverlap=nper//2, window="hann")
    print(f"PSD: nperseg={nper} -> df={fs/nper*1e6:.1f} uHz, "
          f"{len(x)//(nper//2)-1} segments")
    detect(P, fs, nper)

if __name__ == "__main__":
    main()
