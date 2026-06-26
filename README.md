# Seismic-EM-Data-Spacetime-Foam

**Detection of spacetime foam signatures through octave-comb resonance and coherence analysis in seismic and electromagnetic data.**

This repository contains Python tools, test data references, and analysis scripts for searching ultra-low-frequency (ULF) geophysical data for structured resonant signatures predicted by certain spacetime foam / resonant scalar medium models.

## Overview

The project implements two complementary detection strategies:

1. **Locked Octave-Comb Detector** — Searches for statistically significant power excess at frequencies \(f_0 \cdot 2^k\) (octave harmonics), with base frequency \(f_0 \approx 2.005\) mHz (corresponding to the light-travel time across 1 AU).

2. **Pc5 → Seismic Hum Coherence Bridge** — Tests whether Pc5 ULF geomagnetic bursts are followed by persistent, phase-locked amplification in the seismic hum band near 3.7 mHz.

Both methods use rigorous statistical controls (bootstrap null distributions, control frequency bins, and phase-locking metrics) to distinguish potential foam-related signals from conventional geophysical or instrumental effects.

## Repository Contents

| File | Description |
|------|-------------|
| `octave_comb_filter.py` | Core implementation of the octave-comb matched filter and statistical detector. Includes synthetic data generation, Welch PSD computation, whitening, bootstrap-based significance testing, and dual-criterion detection (incoherent sum + multi-tooth coincidence). |
| `run_octave_comb_fdsn.py` | Driver script that fetches real broadband seismic data directly from any FDSN web service (IRIS, GEOFON, ORFEUS, etc.) using ObsPy and runs the octave-comb detector. |
| `foam_coherence_detector.py` | Detector for phase-locked coherence between Pc5 ULF bursts (ground magnetometer) and seismic hum. Includes synthetic validation for "foam-like" vs "null" worlds. |
| `Seismic_Test_Data_Reference.md` | Detailed documentation of test data windows, FDSN download links, station parameters (PFO, BFO), and preliminary results across multiple solar and earthquake events. |
| `Solar_Seismic_Foam_Findings_Writeup.docx` | Detailed writeup of solar flare / geomagnetic storm analysis. |
| `Spacetime_Foam_Report.docx` | Broader theoretical and experimental context report. |
| `octave_comb_demo.png` | Example visualization of octave-comb filter output on synthetic or real data. |

## Key Theoretical Anchor

The base frequency is deliberately set to:
```
f₀ ≈ 2.005 mHz  (= c / 1 AU light-crossing time)
```

This choice anchors the search to a macroscopic solar-system scale, consistent with certain spacetime foam and resonant navigation frameworks.

## Installation

```bash
pip install numpy scipy obspy matplotlib
```

- `obspy` is only required for `run_octave_comb_fdsn.py` (FDSN data access).
- The core filter and coherence detector can run with just NumPy + SciPy.

## Usage

### 1. Run Octave-Comb Detector on Real FDSN Data

```bash
python run_octave_comb_fdsn.py
```

Edit the parameters at the top of the script (`CLIENT`, `NET`, `STA`, `CHA`, `T1`, `T2`, `F0`, etc.) to test different stations, time windows, or channels (VHZ, LHZ, BHZ, etc.).

### 2. Run Pc5–Seismic Hum Coherence Detector

**With real data (CSV format):**
```bash
python foam_coherence_detector.py pc5.csv hum.csv
```

**Synthetic validation only:**
```bash
python foam_coherence_detector.py
```

CSV files should have two columns: `time_seconds, value` at uniform 60-second cadence.

## Current Status & Preliminary Results

Multiple clean solar-driven test windows (including the May 2024 Gannon superstorm and the October 2024 X9.0 flare period) have been analyzed. **No statistically significant locked octave comb** has been detected in these intervals. 

A positive control using the 2006 Tonga Mw 8.0 earthquake correctly identifies known free oscillation modes (e.g., ₀S₁₂) but does not produce an octave-comb structure. These results are documented in `Seismic_Test_Data_Reference.md`.

## Requirements for Real Data Analysis

- Broadband vertical seismic data (VHZ / LHZ channels preferred for low-frequency work)
- Access to any FDSN-compatible data center (no account required)
- Sufficient data length (typically weeks to months) for statistical power

## Related Research

This work forms part of ongoing experimental efforts to detect macroscopic signatures of spacetime foam through resonant geophysical observables. It complements theoretical and hardware development in resonant scalar systems, 4D navigation, and CryoGrid technologies.

## Contributing

This is active research code. Issues, improvements to the statistical methods, additional test windows, and visualization enhancements are welcome.



Would you like a slightly shorter version, a more technical/math-heavy version, or one that includes specific links/badges? I can also adjust the tone if you want it more formal or more exploratory.
