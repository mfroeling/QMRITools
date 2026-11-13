# EPG T2: agent guide

This is a map of the extended phase graph (EPG) T2 code: the EPG section of
[QMRITools/Kernel/RelaxometryTools.wl](../QMRITools/Kernel/RelaxometryTools.wl) (`EPGSignal` to `EPGDictionary`,
about L440–L1340) and the slice profile functions in [SimulationTools.wl](../QMRITools/Kernel/SimulationTools.wl)
(`GetPulseProfile`, `SimulateSliceEPG`). It fits water T2, B1 and the water and fat signal of multi-echo spin-echo
(MESE) data with a two-component EPG model. Read this before changing that code.

- Code state as of 2026-11-11. Line numbers (`L123`) drift over time; grep for the function name if one no longer
  matches.
- Index of all agent docs: [AGENTS.md](../AGENTS.md). House style: [CodeStyle.md](CodeStyle.md). Units:
  [DataConventions.md](DataConventions.md). The main consumer is the BIDS pipeline: [MuscleBidsTools.md](MuscleBidsTools.md)
  (`MuscleBidsProcessI`, method `EPGT2`).
- Method background:
  - Weigel, JMRI 2015, 10.1002/jmri.24619 (EPG formalism, used by `EPGSignal`).
  - Marty et al., NMR Biomed 2016, 10.1002/nbm.3459 (dictionary fit of water T2, B1 and fat).
  - Keene et al., MRM 2020, 10.1002/mrm.28290 (slice profile with fat chemical shift, 90/10 adipose tissue model,
    calibration methods A, B, C). The current code follows this paper (§3, §6).

---

## 1. The pipeline (BIDS method `EPGT2`)

```
mese data {z, echo, y, x}, json (EchoTime s, SliceThickness, MagneticFieldStrength)
  → background mask, NormalizeData
  → GetPulseProfile[ex, ref, SliceRange -> 3 thk, SliceRangeSamples -> 6 thk, FatFieldStrength -> B0]
        half slice profiles of excitation and refocusing angles (0.25 mm steps) + fat shift in samples   (§3)
  → EPGT2Fit[data, 1000 echos, angle, DictT2IncludeWater -> True, EPGFitFat -> False, EPGCalibrate -> True, EPGFitPoints -> 100,
             EPGFatShift -> shift, DictB1Range -> {0.3, 1.4, 0.02}, DictT2Range -> {10, 60, 0.2}, ...]
        CalibrateEPGT2Fit ("2compF" = Keene method C) on subcutaneous fat → T2fat (dictionary search) (§6)
        CreateT2Dictionaryi (water part cached) → brute force dictionary fit → B1 median smoothing → refit (§4, §5)
  → t2w, t2f, b1, wat, fat, fatfr (lipid fraction), res; json gets EPGSettings and EPGCalibration {T2fat, B1}
```

`Settings` may also be two plain numbers (`[90, 180]`). Then no profile and no shift are used, and B1 is ill-posed
(§2, Q14). The BIDS settings are hard-coded in the `EPGT2` branch, only `Settings` and `Shift` come from the config.

## 2. Signal model and conventions

- **Units.** T1, T2 and echo spacing in ms, angles in degrees, B1 as a factor (1 = nominal). B1 scales both the
  excitation and the refocusing angle. `EPGT2Fit` takes `{nEcho, echoSpacing}` or the list of echo times; it then
  assumes the first TE equals the spacing (true for CPMG MESE).
- **State vector.** `{F+0, F-0, Z0, F+1, F-1, Z1, ...}`, real. `MixMatrix` shifts the states, `RotMatrixT` is the RF
  rotation, `MakeDiagMat` repeats a 3×3 block per state. One echo is relax(τ) → shift → RF → relax(τ) → shift
  (`MoveStates`), τ = spacing/2, with T1 recovery to 1 of Z0.
- **CPMG.** The real path is exact CPMG: it matches an independent complex Weigel EPG (excitation phase 90°,
  refocusing 0°) to ~1e-9 for refocusing angles 120–180° and B1 0.7–1. `EPGRefocussingPhase -> 90` gives CP. A
  non-zero phase or off-resonance switches to the complex path (`RotMatrixTI`, `MoveStatesI`).
- **Off-resonance** (`fp` frequency) has no effect on the magnitude. That is correct for a crushed CPMG with equal
  echo spacing.
- **State truncation.** `states = nEcho` below 20 echoes, else `Max[nEcho/2, 20]`. The error is <1e-4 for muscle T2 and
  0.008 for T2 = 200 ms with 50 echoes.
- **Slice profile.** `EPGSignali` with angle lists simulates each profile sample and averages
  `Join[sig, sig[[2;;]]]`, so it expects a **half** profile (centre first), as `GetPulseProfile` returns.
  `ShiftPulseProfile` returns **full** (mirrored) profiles; these are then re-mirrored, which gives the edge sample
  half weight (Q16, negligible).
- **B1 symmetry.** With scalar angles, B1 = b and B1 = 2 − b give identical signals (refocusing 180·b and 360 − 180·b
  are the same rotation in opposite directions). A slice profile breaks this: the difference is 0.06 at B1 0.9 and 0.17
  at 0.7, smaller close to 1. Always fit with slice profiles (Q14).
- **B1 above 1 is weakly determined, even with a profile.** Folding angles above 180° (α → 360 − α), B1 > 1 moves the
  ring where the refocusing angle passes 180° outward but keeps a similar spread of angles; below 1 no part of the slice
  reaches 180°. The signal shape changes 0.03 per 0.1 B1 at 1.0 but 0.005 at 1.3. At typical SNR, voxelwise B1 between
  about 1.05 and 1.4 is noise (simulation: true 1.16 gives 11% at the 1.4 bound). T2 is hardly affected. Hence the
  pooled B1 (§5).
- **B1 below 0.5 is real at 3 T** in the centre of the thighs (dielectric shading): echo2/echo1 of 1.7–1.8, above the
  model maximum at 0.5. The BIDS range starts at 0.3.
- **Output order.** `EPGT2Fit` returns `{{T2map, T2fatMap, B1Map}, {wat, fat, fatMap}, residual}`, plus the calibration
  when `OutputCalibration -> True` and `EPGCalibrate -> True`. `residual` is the RMS over echoes in signal units.

## 3. Fat model and chemical shift

- **Adipose tissue** is 90% lipid and 10% water (Keene §2.2, method C). The constants and the mix live in one place:
  `fatWatFr = 0.1`, `fatWatT2 = 20.` and `MixFatWater[fat, wat]` (L646). The water-in-fat component uses the muscle T1
  and the **unshifted** (water) profile; the lipid uses `T1fat` and the shifted profile. The dictionary uses this model
  when `DictT2IncludeWater -> True`, the calibration methods `"2comp"` and `"2compF"` always use it.
- **Lipid fraction.** With `DictT2IncludeWater -> True`, `EPGT2Fit` splits the fitted fat (adipose) amplitude:
  `wat = wat + 0.1 fat`, `fat = 0.9 fat`. So `fatMap` is the lipid fraction, the quantity Dixon reports, and
  subcutaneous fat reads about 0.9. Before 2026-11 it was the adipose tissue fraction (about 10% higher).
- **Chemical shift in the slice direction.** Fat sees each pulse's slice shifted by ΔL = Δf·B0/(γ·G)
  (Keene §2.1). With different gradients for excitation and refocusing, the fat profiles are shifted relative to each
  other by ΔL_ex − ΔL_ref. `GetPulseProfile` computes it with `FatFieldStrength -> B0`:
  `3.4 B0 42.5775 (thk_ex/BW_ex − thk_ref/BW_ref)` in mm (thk = slice thickness of the pulse, BW in Hz), divided by the
  profile step `SliceRange/SliceRangeSamples/2` (0.25 mm with the BIDS settings), absolute value. It is then the third
  output: `{ex, ref, shift, plots}`.
- **Example.** The BIDS example pulses (`sg_150_100_167` G 1.5991 mT/m, BW 538 Hz; `sinc_centre` G 1.3326 mT/m, BW 340 Hz)
  give 1.276 mm = 5.10 samples at 3 T.
- **Applying it.** `EPGFatShift -> shift` (samples, fractional allowed, sign irrelevant because the profiles are
  symmetric). `ShiftPulseProfile` mirrors both profiles and shifts the refocusing profile of fat by interpolation.
  Leaving the shift out biased muscle T2 at 78% lipid by +6 ms in the phantom (§8).

## 4. Dictionary (`CreateT2Dictionaryi`, `EPGDictionary`)

- Water signals: `T2Range` × `B1Range` (BIDS 251 × 56, T2 10–60 ms as in Keene, B1 down to 0.3 because thigh centres
  at 3 T reach B1 ≈ 0.35–0.5), each with the full profile. This is
  the expensive part (about 20–60 s for 151 × 46, proportionally longer for the BIDS grid) and is cached by `EPGDictionary` (a cached table of EPG signals for any T1 and profile, also used for the calibration lipid signals) on `{echo, T1, profile, t2 values, b1 values, phase}`, so it is
  built once per study.
- Fat signals are rebuilt each call (about 1 s for one T2fat, a few s for a T2fat range). A new calibrated T2fat per
  subject therefore costs nothing extra.
- Dimensions: fixed T2fat → `dict {t2m, b1, echo, 2}`, `vals {t2m, b1, {t2m, b1}}`; T2fat range →
  `dict {t2m, t2f, b1, echo, 2}`, `vals {t2m, t2f, b1, {t2m, t2f, b1}}`.

## 5. `EPGT2Fit` (L857)

1. Optional calibration (§6) sets the fixed T2fat (`t2fval`, the calibration median, not rounded).
2. `EPGFitFat -> False` (BIDS): 2D dictionary with the fixed T2fat. `EPGFitFat -> True` (default): 3D dictionary; the
   first pass uses the T2fat slice nearest to `DictT2fValue`.
3. Fit per voxel, `EPGMethod`:
   - `"dictionary"` (default): brute force over all entries, linear least squares for the water and fat amplitudes
     (`DictionaryMinSearchi` with the precomputed pseudo-inverses).
   - `"dictionaryM"`: integer `NMinimize` over the dictionary. Broken (Q10).
   - `"NLLS"`: `FindMinimum` on T2 and B1 (`NonLinearEPGFiti`), fixed T2fat only (Q9).
4. `EPGSmoothB1` (default True = radius 2, or a radius; volumes and dictionary methods only): instead of the voxelwise
   fit, `DictionaryB1Cost` gives each voxel's cost per B1 (minimum over T2), `BoxSum` (separable 1D sums, 40× faster than `MeanFilter`, same minimum) pools it over the 3D
   neighbourhood (B1 is smooth in-plane and through-slice), and the B1 with the lowest pooled cost is used. Then a refit
   at that fixed B1: water T2 only (fixed T2fat) or water T2 and T2fat (`EPGFitFat -> True`). Background voxels have zero
   cost and do not affect the minimum. Equal noise in all voxels makes the summed cost the joint likelihood.
   - Two grids: the B1 cost uses a coarse T2 grid (same range, step `Max[1, step]` ms), the refit the fine
     `DictT2Range` grid at the chosen B1 only. The minimum over T2 barely depends on the T2 step, so the expensive pass
     shrinks about 5× and a finer B1 or T2 grid stays affordable. Both water dictionaries are cached.
   - Fit quality weighting: each voxel's cost curve is divided by 1 + (r/r_median)², r its own best residual, so voxels
     the model fits badly (vessels, marrow edges, partial volume) count less in the pooling.
   - Why not a median of voxelwise B1 (the method before 2026-11): with a slice profile, B1 above about 1.05 barely changes
     the signal (the profile always keeps a ring near 180°, so the folded angle distribution hardly changes), so voxelwise
     B1 > 1 is mostly noise and a median of noise does not help. Stack T22 of the motion study, slice 15: neighbour
     |ΔB1| 90th percentile 0.18 voxelwise, 0.12 median filtered, 0.04 pooled 5×5 (2D test).
5. Lipid split (§3), then with `EPGRepetitionTime -> TR` (ms, BIDS passes the json TR) the water and fat amplitudes
   are divided by 1 − exp(−(TR − nEcho·ΔTE)/T1) with T1 from `EPGRelaxPars` (assumes Mz ≈ 0 after the echo train).
   Without it the fat fraction is biased up: at TR 2.5 s muscle recovers to about 0.82, fat fully. Then the fat fraction
   is clipped to {−0.1, 1.1}, the amplitudes to {−0.1, 1.5}·max.

## 6. Calibration (`CalibrateEPGT2Fit` L1119)

- **Fat mask.** 4D: last echo > 50 after `NormalizeData` (median first echo = 100), times the outside of the eroded
  first-echo mask, the two largest components. Then at most `EPGFitPoints` (default 50, BIDS 100) voxels, every
  step-th one. Measured on a real stack (2026-11-12, 3158 fat voxels): per-voxel T2fat sd 16 ms, median of 50 voxels
  130.2 vs 132.0 for 452 voxels; the median moves about 2–3 ms with 50 voxels and 1.3–2 ms with 100. The fit takes
  about 3.6 s per voxel on one kernel (parallel in the calibration). Stacks of one session differ by up to 15 ms
  (130–145), more than this noise, so the calibration stays per stack.
- **Threshold depends on the protocol.** It relies on fat being bright at the last echo. A synthetic phantom with
  full T1 recovery finds no fat; real data with a finite TR (muscle saturated) does. Keene uses the mean of the last
  echo instead (Q13).
- **Methods** (Keene §2.6.3). Each fits B1 too:

| `EPGMethodCal` | Keene | Model of the fat voxels | Fitted |
| --- | --- | --- | --- |
| `"1comp"` | A | pure lipid, free amplitude | T2fat, B1, S0 |
| `"2comp"` | B | 90% lipid + 10% water, one free amplitude | T2water-in-fat, T2fat, B1 |
| `"2compF"` (default, BIDS) | C | 90% lipid + 10% water with T2 20 ms, one free amplitude | T2fat, B1 |

- **Fit by dictionary search** (since 2026-11-13, replaces `FindMinimum`): lipid signals from `EPGDictionary` (cached per
  protocol, T2fat 100–300 ms or 20–300 for `"1comp"`, 1 ms step), mixed with the water-in-fat signal at the water signal
  fraction, one free amplitude per entry, best entry per fat voxel, median over voxels. Global, so no start-value
  dependence; 0.6 s per calibration once cached instead of about 35 s.
- **B1 limited to 0.5–1.1 in the calibration.** In subcutaneous fat the summed cost has a second, lower minimum at the
  upper B1 bound (residual about 20% lower at B1 ≈ 1.45 on stack T22). Most likely the model mimics something it lacks,
  such as fat J-coupling, with an unrealistic B1; the neighbouring muscle shows B1 ≈ 1. T2fat barely depends on it (142
  vs 143.5 ms), but the limited range keeps the logged calibration B1 physical (1.02 ± 0.07).
- The water fraction is fixed on purpose (user decision 2026-11-11, and Keene): a fitted fraction is too variable.
- Output `{median, sd}` per parameter. With fewer than 5 fat voxels: message `CalibrateEPGT2Fit::data`, no fit, and the
  defaults (`{150, 1}` for `"2compF"`). BIDS logs the calibrated values and writes them as `EPGCalibration` to the json,
  so a fallback shows as exactly `{150, 1}`.
- Only use `"2comp"`/`"2compF"` together with `DictT2IncludeWater -> True` (and `"1comp"` without), otherwise the
  calibration and the dictionary use different fat models.
- **T2fat precision matters.** Muscle T2 in fatty muscle is sensitive to T2fat: in the phantom a 5 ms error (160 instead
  of 155) moved muscle T2 at 78% lipid from 31.8 to 37 ms. The old `Ceiling[t2fval, 5]` did exactly that (Q7); the value is no longer rounded.

## 7. Known quirks and likely bugs

| # | Where | Status |
| --- | --- | --- |
| Q1 ✔ | `GetPulseProfile` `FatFieldStrength` | **Fixed.** Divided both slice thicknesses by the excitation bandwidth (1.53 mm instead of 1.28 mm for the example pulses). |
| Q2 ✔ | BIDS `EPGT2` shift | **Fixed.** Own formula with hard-coded 3 T, rounded to whole samples, negative shifts dropped. Now uses `GetPulseProfile` with the json field strength (3 T if missing). |
| Q3 ✔ | `shift > 0` guards | **Fixed.** `CalibrateEPGT2Fit` and `CreateT2Dictionaryi` now apply any non-zero shift. |
| Q4 ✔ | water in adipose tissue | **Fixed.** Used `T1fat`; now `T1m`, through `MixFatWater`. |
| Q5 ✔ | `"2comp"`, `"2compF"` | **Fixed.** Free water amplitude, so not Keene B/C and a different fat model than the dictionary. Now fixed 10% water. |
| Q6 ✔ | fat fraction | **Fixed.** Was the adipose tissue fraction; now the lipid fraction (`DictT2IncludeWater -> True`). Changes BIDS results by about −10% relative. |
| Q7 ✔ | `t2fval` rounding | **Fixed.** `Ceiling[t2fval, 5]` (to reuse dictionaries) turned 155.00000006 into 160 and always biased T2fat upward. Now not rounded at all: only the water dictionary is cached, the fat part is rebuilt each call. |
| Q8 ✔ | empty fat mask | **Fixed.** `step = 0` broke the subsampling and the fits ran on nothing. Now no fit, defaults returned. |
| Q9 | `EPGMethod -> "NLLS"` with `EPGFitFat -> True` | **Open.** `met` becomes `"dictionaryM"` but the dictionary is only built when the *option* is not `"NLLS"` (`OptionValue[EPGMethod]` instead of `met`, also in the smoothing check). Crashes. |
| Q10 | `"dictionaryM"` | **Open.** `NMinimize` over `Integers` returns unevaluated (`Set::shape`), every call. |
| Q11 | single voxel, `EPGFitFat -> True` | **Open.** The refit step (`dataf = RotateDimensionsLeft[{datal, B1Int}]`) assumes a volume. |
| Q12 | `CalibrateEPGT2Fit` 3D data | **Open.** `fitData = ... echo - 10^-10` multiplies by the `{nEcho, spacing}` list, and the `GetMaskData` call gets wrong dimensions. Only 4D works. |
| Q13 | calibration fat mask | **Open.** Absolute threshold 50 after normalisation, protocol dependent (§6). |
| Q14 | scalar angles | **Intentional.** B1 mirror b ↔ 2 − b is physics; always use slice profiles. Not a bug. |
| Q15 ✔ | B1 smoothing | **Fixed.** The median filter of voxelwise B1 set the background to 1 and could not remove noise on the B1 > 1 plateau. Replaced by pooled cost (§5). The plateau itself is physics: pooled patches can still sit at the upper bound. |
| Q16 | shifted profiles | **Open, negligible.** `EPGSignali` re-mirrors the full profiles from `ShiftPulseProfile`: edge sample half weight. |
| Q17 | dead code, usage text | **Open.** `wMat`/`dictMatW` computed and unused; `WaterFatShift`, `WaterFatShiftDirection` do nothing; usage of `EPGRelaxPars` (says `{1400,365,137}`), `EPGFitPoints` (says 200), `NonLinearEPGFit` (argument order) and `EPGSignal` ("b1 between 0 and 1") are stale. |

Ask the user before fixing an open one.

## 8. Validation (synthetic phantom, 2026-11-11)

Pulses from the BIDS example, thk 8 mm, 17 echoes × 10 ms, B1 0.85, muscle T2 30, T2fat 155, adipose 90/10, shift
applied in the data, muscle scaled 0.6 (TR saturation, §6). 4 slices with lipid fractions 0 / 0.265 / 0.562 / 0.783.

| Run | Calibration | Muscle T2 per slice | Lipid fraction per slice | Subcutaneous |
| --- | --- | --- | --- | --- |
| BIDS settings | 155.0 ms, B1 0.850 | 30.4, 30.4, 30.8, 31.8 | −0.001, 0.266, 0.566, 0.787 | 0.908 |
| T2fat 160 (old rounding) | – | 30.4, 30.8, 32.0, 37.0 | −0.001, 0.258, 0.550, 0.765 | 0.869 |
| no shift (T2fat 160) | – | 30.4, 31.0, 33.2, 43.2 | −0.001, 0.253, 0.548, 0.767 | 0.883 |

- The muscle T2 offset (0.4–1.8 ms) comes from the B1 grid: 0.85 is not on the 0.02 grid, the fit lands on 0.86.
- `"2comp"` recovers T2water-in-fat 20.0, T2fat 155.0, B1 0.850.
- A negative shift gives results identical to the positive one.
- Timing (phantom, `FindMinimum` calibration, before 2026-11-13): `"2compF"` ~100 s, `"2comp"` ~4.5 min. The dictionary calibration takes under a second once the lipid table is cached.

**Real data, motion study stack T22 (slices 10–15, 2026-11-13).** Old median-filtered B1 vs pooled weighted B1 (range
0.3–1.4): neighbour |ΔB1| 90/99% in-plane 0.10/0.40 → 0.02/0.04, through-slice 0.16/0.40 → 0.02/0.06, voxels at a bound
6.3% → 0.4%, muscle T2 median 28.2 → 27.6 ms. Coarse (1 ms) vs fine (0.2 ms) T2 grid for the B1 pass: identical B1 in
99.1% of voxels, T2 identical in ≥ 99%; the cost pass takes 8.2 s instead of 39.1 s for 3 slices (4 kernels).

## 9. Recipes

- **Phantom test.** Build signals with `EPGSignal` on the profiles from `ShiftPulseProfile[GetPulseProfile[...][[1;;2]],
  shift]`: muscle `0.6 EPGSignal[..., {1400, 30}, ang, b1]`, adipose
  `0.9 EPGSignal[..., {365, t2f}, angS, b1] + 0.1 EPGSignal[..., {1400, 20}, ang, b1]`. Use a muscle disk with a thin
  adipose ring (for example radius 18 and 21 in a 48 × 48 image) so the fat mask finds the ring. Expected lipid
  fraction for adipose amplitude a: `0.9 a/(0.6 (1 − a) + a)`.
- **Distribute private functions by full name in test scripts** (`` QMRITools`RelaxometryTools`Private`DictionaryB1Cost ``,
  plus `ErrorC`, `LeastSquares2C`). Distributing an alias (`f = ToExpression[...]`) does not send the definitions: the
  subkernels then return the unevaluated call with the whole dictionary inside, and the main kernel runs out of memory
  (30 GB for one slice) and does the work itself.
- **Profiling.** `` QMRITools`$Debug = {"EPG"} `` prints a timestamped line (with memory) at each `EPGT2Fit` step, when
  a water dictionary is built, and in the calibration. Full stack T22 (16 × 160 × 160, 13 echoes, 16 kernels,
  2026-11-13): calibration 35 s, fine water dictionary 40 s and coarse 8 s (first stack only), B1 cost 12 s, pooling
  4 s, refit 11 s; 111 s for the first stack, 70 s for later ones. Old code (10-2025): 45–85 s per stack.
- **Keep tests small.** Each slice-profile fit takes minutes; launch 4 kernels, one case per run. The water dictionary
  cache only lives within one kernel session.
- **Check a fit model change** against an independent EPG first (complex Weigel implementation, full states): the real
  path must match it to ~1e-9.

## 10. Function index

| Function | File, line | Role |
| --- | --- | --- |
| `EPGSignal`, `EPGSignali` | Relaxometry L458, L464, L469 | EPG echo train, single angles or profile |
| `MixMatrix`, `RotMatrixT(I)`, `MoveStates(I)`, `MakeDiagMat` | Relaxometry L531–L588 | EPG building blocks (compiled) |
| `ShiftPulseProfile` | Relaxometry L603 | mirror profiles, shift the fat refocusing profile |
| `ErrorFunc2M` | Relaxometry | cost function of the NLLS fit |
| `fatWatFr`, `fatWatT2`, `MixFatWater` | Relaxometry L646–L650 | adipose tissue model |
| `LeastSquaresC`, `PseudoInverseC`, `LeastSquaresErrorC`, ... | Relaxometry L658– | compiled linear amplitude fits |
| `NonLinearEPGFit(i)` | Relaxometry L689 | NLLS fit of T2 and B1 |
| `DictionaryMinSearch(i)` | Relaxometry L724–L812 | brute force and `NMinimize` dictionary search |
| `DictionaryB1Cost` | Relaxometry, after the brute force search | per-voxel cost per B1 (minimum over T2), for pooled B1 smoothing |
| `EPGT2Fit` | Relaxometry L857 | main fit |
| `CalibrateEPGT2Fit` | Relaxometry L1119 | T2fat calibration on subcutaneous fat |
| `CreateT2Dictionary(i)` | Relaxometry L1241, L1247 | dictionary (fat part uncached) |
| `EPGDictionary` | Relaxometry, after `CreateT2Dictionaryi` | cached table of EPG signals `{b1, t2, echo}` (water, lipid), distributes its own definitions |
| `GetPulseProfile`, `Pulses` | Simulation L548, L568, L390 | Bloch slice profiles, fat shift |
| `SimulateSliceEPG` | Simulation L653 | profile and EPG plot, optional fit report |
