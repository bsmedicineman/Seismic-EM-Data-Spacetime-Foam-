# Seismic Octave-Comb Tests — Data Links and Event Parameters

**Detector target:** locked octave comb at f₀ = 2.005 mHz (= c / 1 AU light-crossing) and octaves f₀·2ᵏ.
**Primary station:** II.PFO — Pinyon Flat Observatory, California (GSN/IDA, very-broadband).
**Cross-check station:** II.BFO — Black Forest Observatory, Germany.
**Data service (all links):** IRIS/EarthScope FDSN — `https://service.iris.edu/fdsnws/dataselect/1/query`
(The same query string works at any FDSN node, e.g. GEOFON, ORFEUS/EIDA, RESIF, by swapping the host.)

### Channels and testable teeth
| Channel | Sampling | Nyquist | Comb teeth tested | Control teeth (non-octave) |
|---|---|---|---|---|
| VHZ | 0.1 Hz | 50 mHz | 5: 2.0, 4.0, 8.0, 16.0, 32.1 mHz | 6.0, 10.0, 12.0, 14.0 mHz (3,5,6,7×f₀) |
| LHZ | 1.0 Hz | 500 mHz | 8: 2, 4, 8, 16, 32, 64, 128, 257 mHz | same ratios, extended |

---

## Window 1 — PFO 2006 (Tonga positive control)
**Download:**
```
https://service.iris.edu/fdsnws/dataselect/1/query?net=II&sta=PFO&loc=00&cha=VHZ&start=2006-04-25T00:00:00&end=2006-05-14T00:00:00
```
**Request parameters:** net=II, sta=PFO, loc=00, cha=VHZ, 2006-04-25 → 2006-05-14 (~20 d), ~2 MB.

**Event in window — earthquake:**
- **2006 Tonga, Mw 8.0** — origin 2006-05-03 15:26:40 UTC; epicenter 20.187°S, 174.123°W (~47 km SSE of Pangai, Ha'apai, Tonga); depth 55 km. USGS-NEIC body/surface magnitudes mb 7.2 / MS 7.9; intraslab compressional thrust.

**Result:** post-quake the 2 mHz tooth reached **6× background** — the real fundamental free oscillation ₀S₁₂ — but as a single line, not a comb. Pre-quake baseline quiet. Positive control: detector sees real physics, coincidence gate refuses to over-read it. **No comb.**

---

## Window 2 — BFO 2023 cross-check (VHZ)
**Download:**
```
https://service.iris.edu/fdsnws/dataselect/1/query?net=II&sta=BFO&loc=00&cha=VHZ&start=2023-06-15T00:00:00&end=2023-06-16T00:00:00
```
**Request parameters:** net=II, sta=BFO, loc=00, cha=VHZ, 2023-06-15 (1 day).
**Context:** independent quiet baseline. Nearest sizeable event was a Mw 7.2 in the Tonga region (~2023-06-15 18:06 UTC) — below the threshold for strong low-frequency mode excitation. **No comb.**

## Window 3 — BFO 2023 cross-check (LHZ, all 8 teeth)
**Download:**
```
https://service.iris.edu/fdsnws/dataselect/1/query?net=II&sta=BFO&loc=00&cha=LHZ&start=2023-06-15T00:00:00&end=2023-06-16T00:00:00
```
**Request parameters:** net=II, sta=BFO, loc=00, cha=LHZ, 2023-06-15 (1 day). Reaches all 8 teeth to 257 mHz. **No comb.**

---

## Window 4 — PFO May 2024 (Gannon superstorm — flare test)
**Download:**
```
https://service.iris.edu/fdsnws/dataselect/1/query?net=II&sta=PFO&loc=00&cha=VHZ&start=2024-03-12T00:00:00&end=2024-07-09T00:00:00
```
**Request parameters:** net=II, sta=PFO, loc=00, cha=VHZ, 2024-03-12 → 2024-07-09 (~120 d). Split at 2024-05-10.

**Events in window — solar:**
- **Active region AR3664 (NOAA 13664)** produced **12 X-class flares, 1–15 May 2024**, plus ~7–10 halo CMEs.
- **X8.7 flare — 2024-05-14 (~16:51 UTC)** — the largest of Solar Cycle 25 at the time, as AR3664 rotated past the west limb; R3 radio blackout.
- **G5 "Gannon" geomagnetic superstorm, 10–11 May 2024** — shock arrival (SSC) 2024-05-10 17:05 UTC; **peak Dst ≈ −412 nT at ~02 UTC on 11 May**; Kp = 9 (twice); strongest storm since the 2003 Halloween storms.

**Earthquake screening:** **no M8 anywhere in 2024** — making this one of the two cleanest possible solar-driving tests.
**Result (after de-glitching telemetry spikes to 10⁸ counts):** comb band statistically identical before vs after the flare. **No comb; flares produce nothing.**

---

## Window 5 — PFO October 2024 (X9.0 flare + severe storm — flare test)
**Download:**
```
https://service.iris.edu/fdsnws/dataselect/1/query?net=II&sta=PFO&loc=00&cha=VHZ&start=2024-08-09T00:00:00&end=2024-12-06T00:00:00
```
**Request parameters:** net=II, sta=PFO, loc=00, cha=VHZ, 2024-08-09 → 2024-12-06 (~120 d). Split at 2024-10-03.

**Events in window — solar:**
- **X9.0 flare — 2024-10-03 (~12:18 UTC)** — from NOAA AR13842; the largest flare of Solar Cycle 25.
- **Severe geomagnetic storm, 10–11 October 2024** — **Dst ≈ −340 nT (G4)**, driven by an X1.8 flare/CME on 8–9 Oct; aurora seen to low latitudes.

**Earthquake screening:** no M8 in window. Clean solar-driving test.
**Result:** comb band flat before vs after the flare; every tooth ~1× background. **No comb.**

---

## Window 6 — PFO April–May 2026 (recent baseline; Sanriku M7.5)
**Download:**
```
https://service.iris.edu/fdsnws/dataselect/1/query?net=II&sta=PFO&loc=00&cha=VHZ&start=2026-04-01T00:00:00&end=2026-06-01T00:00:00
```
**Request parameters:** net=II, sta=PFO, loc=00, cha=VHZ, 2026-04-01 → 2026-05-31 (61 d). Split at 2026-04-20.

**Event in window — earthquake:**
- **2026 Sanriku, off Honshu, Japan (Japan Trench)** — origin 2026-04-20 07:53 UTC; ~100 km ENE of Miyako, Iwate; depth ~35 km. **Mw 7.4 (USGS) / 7.5 (GCMT) / MJMA 7.5→7.7**; thrust on the subducting Pacific plate interface. USGS event id us6000sri7.

**Result:** weak ~2× perturbation near 2 mHz after the quake — far below Tonga's 6×, proportionate to the smaller moment. **No comb.**

---

## Contamination screening (great earthquakes checked across the program)
A great earthquake rings the free-oscillation band for 1–2 weeks at the comb frequencies, so each window was screened before it was trusted.

| Event | Date (UTC) | Magnitude | Role |
|---|---|---|---|
| Sumatra–Andaman | 2004-12-26 | Mw 9.1 | historical — screened OUT of candidate windows |
| Tonga | 2006-05-03 | Mw 8.0 | **used as positive control (Window 1)** |
| Indian Ocean (off Sumatra) | 2012-04-11 | Mw 8.6 | historical — screened OUT |
| Chiapas, Mexico | 2017-09-08 | Mw 8.2 | historical — screened OUT |
| (all of 2024) | — | no M8 | why May & Oct 2024 are the cleanest flare windows |
| Kamchatka | 2025-07-30 | Mw 8.8 | sits BETWEEN the 2024 and 2026 windows; fully decayed before either window used |
| Sanriku, Japan | 2026-04-20 | Mw 7.4–7.5 | moderate; present in Window 6 |

**Net result across all windows, both stations, both channels, every before/after split:** no octave comb. The band responds to earthquakes in proportion to seismic moment and to solar flares not at all.
