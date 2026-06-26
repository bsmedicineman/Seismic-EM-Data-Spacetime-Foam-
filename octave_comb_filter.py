# -*- coding: utf-8 -*-
"""
Octave-comb matched filter for a locked-frequency resonance search (revised).

Design corrections after demonstration exposed two failure modes:
  (1) A concentrated amplitude template degenerates into a single-line detector.
      FIX: detection requires BOTH an incoherent SUM and a multi-tooth COINCIDENCE
      count, so one strong lone line at f0 cannot fake a comb. Amplitude
      "possibilities" are characterized AFTER detection by fitting the octave
      slope from the data, not used as a detection knob.
  (2) The analytic chi-square null is optimistic (~2x too many false alarms).
      FIX: a DATA-DRIVEN bootstrap null drawn from the off-comb whitened spectrum,
      which captures the true (heavier) tails and calibrates correctly.

Locked comb: teeth at f0 * 2**k. f0 is frozen, so there is no frequency scan.
Non-octave integer harmonics (3,5,6,7 x f0) are evaluated as a control to tell an
OCTAVE comb apart from a generic INTEGER-harmonic comb.

Plug a real Welch PSD into run_detector() to use this for the experiment.
"""
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

F0   = 2.005e-3          # Hz, frozen fundamental (c / 1 AU light-crossing)
FS   = 1.0               # Hz (superconducting gravimeter)
L    = 8192              # Welch segment length -> df = FS/L
NSEG = 48                # independent segments -> K averages
K    = NSEG
ALPHA = 0.05            # per-statistic false-alarm target

# ----------------------------- synthetic data -------------------------------
def synth_series(seed, inject=None):
    g = np.random.default_rng(seed)
    N = L * NSEG
    t = np.arange(N) / FS
    w = g.standard_normal(N)
    fr = np.fft.rfftfreq(N, 1 / FS); f = fr.copy(); f[0] = fr[1]
    shape  = 1.0 / np.sqrt(f / fr[1])                                  # red (1/f power)
    shape += 2.5 * np.exp(-0.5 * ((np.log(f) - np.log(0.18)) / 0.25) ** 2)  # microseism
    shape += 0.15                                                     # white floor
    x = np.fft.irfft(np.fft.rfft(w) * shape, n=N); x /= x.std()
    if inject:
        for fi, ai in inject.items():
            x += ai * np.sin(2 * np.pi * fi * t + g.uniform(0, 2 * np.pi))
    return x

def welch_psd(x):
    return signal.welch(x, fs=FS, nperseg=L, noverlap=0, window="hann",
                        detrend="constant")

# ----------------------------- whitening ------------------------------------
def whiten(P, protect_bins, halfwin=60, guard=3):
    n = len(P); keep = np.ones(n, bool)
    for b in protect_bins:
        keep[max(0, b - guard):min(n, b + guard + 1)] = False
    B = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - halfwin), min(n, i + halfwin + 1)
        seg = P[lo:hi][keep[lo:hi]]
        B[i] = np.median(seg) if seg.size >= 5 else np.median(P[lo:hi])
    return P / B, keep

# ----------------------------- comb machinery -------------------------------
def comb_bins(f0):
    df = FS / L; bins, k = [], 0
    while True:
        b = int(round(f0 * 2 ** k / df))
        if b >= L // 2: break
        bins.append(b); k += 1
    return np.array(bins)

def integer_control_bins(f0, mults=(3, 5, 6, 7)):
    df = FS / L
    return np.array([int(round(f0 * m / df)) for m in mults
                     if int(round(f0 * m / df)) < L // 2])

# ----------------------------- detector -------------------------------------
def run_detector(P, f0=F0, nboot=200000, label=""):
    octave  = comb_bins(f0)
    integer = integer_control_bins(f0)
    m = len(octave)
    rw, keep = whiten(P, np.concatenate([octave, integer]))

    # off-comb whitened powers used as the empirical null pool
    lo = 60; pool = rw[lo:len(rw) - 60][keep[lo:len(rw) - 60]]
    tau = np.quantile(pool, 1 - ALPHA)            # per-tooth coincidence threshold

    r_teeth = rw[octave]
    S_sum   = float(r_teeth.sum())                # incoherent sum (sensitivity)
    C_cnt   = int(np.sum(r_teeth > tau))          # coincidence count (specificity)

    g = np.random.default_rng(7)
    draws = g.choice(pool, size=(nboot, m), replace=True)
    S_null = draws.sum(1)
    C_null = (draws > tau).sum(1)
    p_sum  = float(np.mean(S_null >= S_sum))
    p_cnt  = float(np.mean(C_null >= C_cnt))
    detect = (p_sum < ALPHA) and (p_cnt < ALPHA)  # BOTH gates required

    # amplitude "possibility" characterized from data: slope of log power vs octave
    k_idx = np.arange(m)
    above = r_teeth > 1.0
    slope = np.polyfit(k_idx[above], np.log(r_teeth[above]), 1)[0] if above.sum() >= 2 else np.nan

    oc = float(r_teeth[1:].mean()); it = float(rw[integer].mean())
    return dict(label=label, octave=octave, integer=integer, m=m, tau=tau,
                r_teeth=r_teeth, r_integer=rw[integer], S_sum=S_sum, C_cnt=C_cnt,
                p_sum=p_sum, p_cnt=p_cnt, detect=detect, slope=slope,
                S_thr=float(np.quantile(S_null, 1 - ALPHA)),
                C_thr=float(np.quantile(C_null, 1 - ALPHA)),
                S_null=S_null, oc=oc, it=it, df=FS / L)

def show(r):
    df = r["df"]
    print(f"\n=== {r['label']} ===")
    print("whitened power at octave teeth (1.0 = local noise; * = above tooth thr "
          f"{r['tau']:.2f}):")
    for b, v in zip(r["octave"], r["r_teeth"]):
        star = " *" if v > r["tau"] else ""
        print(f"   {b*df*1e3:8.2f} mHz : {v:7.2f}{star}")
    print(f"incoherent SUM     S = {r['S_sum']:7.2f}  (95% null thr {r['S_thr']:.2f})"
          f"  p_sum = {r['p_sum']:.4f}")
    print(f"coincidence COUNT  C = {r['C_cnt']:7d}  (95% null thr {r['C_thr']:.0f})"
          f"  p_cnt = {r['p_cnt']:.4f}")
    print(f"octave-vs-integer    : octave(2f..) mean {r['oc']:.2f} vs "
          f"non-octave integer mean {r['it']:.2f}")
    verdict = "COMB DETECTED" if r["detect"] else "no comb (consistent with null)"
    extra = f" | data octave slope = {r['slope']:+.2f}/octave" if r["detect"] else ""
    print(f"VERDICT: {verdict}  (needs BOTH gates){extra}")

if __name__ == "__main__":
    octave = comb_bins(F0); df = FS / L
    print(f"f0 = {F0*1e3:.3f} mHz | df = {df*1e6:.1f} uHz | Nyquist = {FS/2} Hz")
    print(f"locked octave teeth ({len(octave)}): "
          f"{[round(b*df*1e3,2) for b in octave]} mHz")

    rA = run_detector(welch_psd(synth_series(11))[1], label="A. Pure colored noise (null)")
    show(rA)

    inj = {F0 * 2 ** k: 0.060 * 0.62 ** k for k in range(len(octave))}
    rB = run_detector(welch_psd(synth_series(12, inject=inj))[1],
                      label="B. Injected octave comb (decreasing, amplitude unknown)")
    show(rB)

    rC = run_detector(welch_psd(synth_series(13, inject={F0: 0.060}))[1],
                      label="C. Single f0 line only (specificity test)")
    show(rC)

    # calibration of the JOINT detector (both gates) over many nulls
    Nrep = 200; fa = 0
    for s in range(Nrep):
        fa += run_detector(welch_psd(synth_series(2000 + s))[1], nboot=40000)["detect"]
    print(f"\n=== calibration ===")
    print(f"JOINT false-alarm rate over {Nrep} null trials: {fa/Nrep:.3f}  (target <= {ALPHA})")

    # figure
    fB, PB = welch_psd(synth_series(12, inject=inj))
    fig, ax = plt.subplots(2, 1, figsize=(9, 7.2))
    ax[0].loglog(fB[1:], PB[1:], lw=0.7, color="#555", label="response PSD (with comb)")
    for b in rB["octave"]:
        ax[0].axvline(b * df, color="#1F3864", lw=1.0, ls="--", alpha=0.85)
    ax[0].axvspan(0.13, 0.30, color="#f4d4d4", alpha=0.5, label="microseism band")
    ax[0].set_xlabel("Frequency (Hz)"); ax[0].set_ylabel("PSD")
    ax[0].set_title("Locked octave comb (dashed = teeth at f$_0\\cdot2^k$)", color="#1F3864")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[1].hist(rA["S_null"], bins=120, color="#cfd8e8", density=True, label="null sum S")
    ax[1].axvline(rA["S_thr"], color="#2e7d32", lw=1.5, ls="--", label=f"95% thr {rA['S_thr']:.1f}")
    ax[1].axvline(rB["S_sum"], color="#b00020", lw=2.0, label=f"comb S={rB['S_sum']:.0f} (p={rB['p_sum']:.3f})")
    ax[1].axvline(rA["S_sum"], color="#888", lw=2.0, label=f"noise S={rA['S_sum']:.1f} (p={rA['p_sum']:.2f})")
    ax[1].set_xlabel("incoherent comb sum"); ax[1].set_ylabel("density")
    ax[1].set_title("Data-driven null vs observed", color="#1F3864"); ax[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig("octave_comb_demo.png", dpi=160, facecolor="white")
    print("\nsaved figure: octave_comb_demo.png")
