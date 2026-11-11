# DixonTools: agent guide

This is a map of [QMRITools/Kernel/DixonTools.wl](../QMRITools/Kernel/DixonTools.wl) (~1570 lines). The file holds
chemical-shift-encoded water-fat separation: the field-map initialisation (`DixonPhase`), the iterative IDEAL fit with a
multi-peak or composition-based fat model (`DixonReconstruct`), fat fraction calculation, signal simulation, echo-time
optimisation and phase unwrapping. Read this before changing that file.

- Code state as of 2026-11-10. Line numbers (`L123`) drift over time; grep for the function name if one no longer
  matches.
- Index of all agent docs: [AGENTS.md](../QMRITools/AGENTS.md). House style: [CodeStyle.md](CodeStyle.md). Units and
  array layouts: [DataConventions.md](DataConventions.md) §5. The main consumer is the BIDS pipeline:
  [MuscleBidsTools.md](MuscleBidsTools.md) (`MuscleBidsProcessI`, Dixon branches).
- Method background: Froeling & Heskamp, NMR Biomed 2024, 10.1002/nbm.5217 (fat models, constrained vs unconstrained
  initial phase). The code has been updated since that paper (see Q19).
- `CholeskySolve` (used inside `DixonFitC`) lives in **GeneralTools.wl**. `T2Fit` is in RelaxometryTools,
  `PCADeNoise` in DenoiseTools.

---

## 1. The pipeline (BIDS method `Dixon-A`, raw complex echoes)

```
real, imag {z, echo, y, x}
  → background mask (Mask 15, MaskComponents -> 2, dilation 1) → PCADeNoise (complex, patch, kernel 5)
  → DixonPhase       : start maps {b0, t2star, phi0, phib} + echoes used           (§3)
  → DixonReconstruct : IDEAL per voxel, "CallDB" (>6 echoes) or "Fixed", all phases (§4)
       → DixonFit/DixonFitC (complex amplitudes; constrained mode: two stages)
       → R2 clip, φ0 = Arg(W+F), median filter field maps, DixonFitFC refit
       → in/out-phase images, proton-density conversion, DixonToPercent
  → DixonPrecessionCheck on the fitted FF (log + json)                          (Q21)
  → export + json (source json + config Process block + version + Dixon settings, MuscleBidsTools.md §6)
```

On real data (320 × 320 × 31, 10 echoes, 3 T) one stack takes about 2 min: denoise 25–55 s, `DixonPhase` 20–40 s,
`DixonReconstruct` 20–30 s with `Parallelize -> True`, peak memory about 4.4 GB.

## 2. Signal model, conventions and units

**Signal.** `S(t_n) = exp(ψ(t_n)) · Σ_species ρ_s A_s(t_n)` with (Froeling 2024 eq. 7)

`ψ(t) = −R2'·t + 2πi(B0·t + θ0·(−1)^n + φ0 + θ1·(−1)^n·t)`

| Parameter (index in `DixonPhases`) | Symbol | Unit | Output name (BIDS) |
| --- | --- | --- | --- |
| 1 R2 (R2* or R2' with `DixonFixT2`) | R2' | 1/s | `r2star`, `t2star` = 1/R2 in s |
| 2 B0 | B0 | Hz | `b0` |
| 3 bipolar phase | θ0 | cycles | `phbp` |
| 4 initial phase | φ0 | cycles | `phi` |
| 5 bipolar phase drift | θ1 | Hz | `phbpt` |

- **Echo times are in seconds.** Frequencies are `γ · B · ppm` in Hz with water at 0 (`dixFreq` holds ppm − 4.7).
  `DixonPrecessions` (default +1) sets the rotation sense of the fat peaks relative to water. It is vendor/convention
  dependent and is not detected yet (Q21).
- `n` runs from 1, so echo 1 has bipolar sign −1 (`vp = (-1)^Range[ne]` in both `DixonPhase` and `DixonReconstruct`).
- **Amplitude matrix.** `matA` columns: water, then the fat basis. Complex amplitudes are fitted in real block form
  `[[Re A, −Im A], [Im A, Re A]]` (`matAc`); constrained (real) amplitudes as `Join[Re A, Im A]`.
- **Phase matrices.** `matC` (complex, `{1, I, I, I, I} vec`) demodulates; `mat` (real, `{1,-1,-1,-1,-1} vec` stacked
  with `vec`) times the fitted signal (`sol`) or its real/imag swap (`solR`) gives the Jacobian columns of ψ. Together
  with `matA` that is the IDEAL matrix B (Yu 2007 eq. A.4/A.5). Both are built for the selected phases only.

**Fat models (`GenerateAmps`).** The 10 peaks of `dixFreq` with proton counts from Hamilton (Froeling 2024 Fig. 1D):
`{9, 6(cl−4)−8db+2idb, 6, 4(db−idb), 6, 2idb, 2, 2, 1, 2db}`, total `6cl − 2db + 2`. A string model gives a linear basis:
water column (amplitude `Round[100/(nVar+1)]`, conditioning only), the constant fat column C0 and one column per free
variable. "Fixed" = {17.31, 2.62, 0.46}; "CallDB" = cl = 16.29 + 0.39 db, idb = −0.73 + 0.45 db (newer than the paper,
Q19). A matrix of amplitudes instead of a string gives the classic fixed multi-peak model (`dixAmp`). After the fit
water and fat are converted to proton density with the model's total proton count; extra signals (e.g. `dbond`) are
ratios to the fat amplitude.

## 3. `DixonPhase` (L242): start maps

Returns `{{b0 (Hz), t2star (s), phi0 (cycles), phib (cycles)}, {e1, e2, n}}`. It iterates (`MaxIterations` 15) over three
estimates, each removed from the data before the next round:

1. **B0** from the phase difference of echoes `e1, e2` (`FindInPhaseEchos`, both near in-phase for the multi-peak fat).
   For the motion protocol (TE 1.24 ms + 1 ms) that is echoes {2, 4}, dt = 2 ms. The phase is unwrapped with the
   weighted least-squares DCT solver (`UnwrapDCTi`, mask weights).
2. **Bipolar phase θ0** from echo triplets, `S_n² / (S_{n−1} S_{n+1})`, in which B0 cancels.
3. **Initial phase φ0** from Bydder (2011) eq. 10, `½ arg{(AᴴS)ᵀ Re(AᴴA)⁻¹ (AᴴS)}`; `ph1` holds 2φ0.

Then, after the loop:

- **`RestoreDC`** (inside the unwrap function, every iteration): the DCT solution has zero mean, which loses the global
  phase. The constant is restored as the circular mean of (wrapped − unwrapped) inside the mask.
- **`SelectAlias`** (L370): B0 from a pair 2 ms apart is only known modulo 1/dt = 500 Hz. Per connected object of the
  mask (the two legs; the BIDS mask has at most 2 components) it compares B0 + k/dt, k = −1, 0, 1, by the residual of a
  water-fat model with free complex amplitudes and keeps the lowest. A 0.1% penalty per unit of |mean field| prefers the
  smallest field and breaks exact ties (with even echo spacing ±1/ΔTE aliases are identical).
- φ0 gets one global sign correction (phase of the summed W+F), not a per-voxel one: per voxel it soaks up the B0 fat
  bias and noise.
- T2* from a log fit of the magnitude from the first echo after 0.75 · in-phase period.

Known limitation: echoes 2 and 4 are 0.15 cycles apart in fat phase, so B0 in fat starts ≈ 70 Hz off (Q16).

## 4. `DixonReconstruct` (L480) and the IDEAL fit

- Overloads fill missing start maps with 0. Options: §9 function index, usage strings in the file.
- `sel` = requested phases (`DixonPhases`, default R2 and B0 only). **With complex amplitudes
  (`DixonConstrainPhase -> False`, default) φ0 is not fitted** (`fitSel` drops index 4): a common phase can sit in φ0 or
  in the amplitude phases, which makes B singular (condition ≈ 10¹⁶). The φ0 output is then `TotalPhase` = Arg(W+F) at
  TE = 0, so the output structure is unchanged.
- **Constrained mode** (real amplitudes) first runs the complex-amplitude fit and uses its R2, B0, θ and Arg(W+F) as start
  values. Real amplitudes cannot absorb a phase error, and from the `DixonPhase` start (fat-biased B0) the fit diverged
  in fat (18% swaps in the phantom before, 0–2 after).
- `DixonFit` (L677) builds `matC`/`mat` for a phase selection and runs `DixonFitC` either directly (Listable) or with
  `ParallelMap` over slices (`Parallelize`, default True).
- **`DixonFitC`** (L707) per voxel: demodulate, ρ = A⁺ s, residual, B = [Jacobian of ψ · fitted signal, A], solve the
  normal equations with `CholeskySolve` and keep the first `Length[sel]` entries, update. Stop when
  `Total[{1,1,100,100,100}[[sel]] · |Δφ|] < DixonTolerance` (R2 in 1/s, B0 in Hz, phases ×100 in cycles) or after
  `DixonIterations`. The convergence check comes before the update, so the returned ρ/rms belong to the returned φ.
  Compiled with `"InlineExternalDefinitions" -> True`: no `MainEvaluate`, no `CompiledFunctionCall` (CodeStyle.md).
- Output filter (`DixonFilterOutput`, default True): R2 clipped to [1, 200] s⁻¹, field maps median filtered
  (`DixonFilterType`, `DixonFilterSize`), then `DixonFitFC` refits the amplitudes with the final maps. `res` is always
  scaled like the signals.
- Output: `{{watFr, fatFr}, {watSig, fatSig, model params…}, {inph, outph}, {{b0, phb, ph0, phbt}[[sel]], {t2, r2}}, itt,
  res}` (the phase list holds only the selected phases).

## 5. Other functions

- **`DixonToPercent`** (L842): magnitude discrimination (Liu 2007) in both modes: the minor component is
  1 − major/|W+F|, which can be slightly negative from noise. `clip` only sets the range ({0, 1} or {−0.1, 1.1}). For ROI
  means use `DixonClipFraction -> False`; clipping at 0 biases low FF upward.
- **`SimulateDixonSignal`** (L883): `fr` is the **water** fraction, echo times in s, T2 in s (Q18).
- **`FindInPhaseEchos`** (L906): ranks echo pairs by in-phase-ness, phase difference, echo distance and first echo;
  falls back to `{1, n}` without a phase check.
- **`OptimizeDixonEcho`**: Manipulate with the condition number of the amplitude matrix A only (not the full Jacobian).
- **Unwrapping**: `UnwrapDCT`/`UnwrapDCTi` (least squares, Ghiglia-Romero PCG with DCT preconditioner, zero-mean
  output) are used by `DixonPhase`. `Unwrap`/`UnwrapSplit` (Herráez/Abdul-Rahman quality-guided path) are only used by
  BIDS `Dixon-B`.

## 6. Why water-fat swaps happen (verified on a synthetic phantom)

- **IDEAL is a local Gauss-Newton fit per voxel.** Its capture range is about half the fat shift (≈ ±200 Hz at 3 T).
  Everything depends on the start B0.
- **With ΔTE = 1 ms a 500 Hz B0 error is exactly a sign flip on every other echo**, which the bipolar phase models
  exactly (θ0 = 0.25 cycles plus a φ0 shift). With bipolar terms fitted, IDEAL cannot correct a wrong alias; only
  `DixonPhase` can (θ0 from triplets is B0-free there).
- Before the 2026-11 fixes, `DixonPhase` lost the global B0 offset (zero-mean unwrap): global offsets of 200 Hz gave 35%
  swapped voxels, ≥ 260 Hz about 98%. Now offsets up to 300 Hz, −434 Hz (scanner locked on fat) and a 300 Hz jump
  between the legs all give the correct B0.
- Remaining: isolated voxels (0.2–0.4% in the phantom) at fat/muscle edges from the B0 fat bias (Q16), and local swaps
  inside one object where the field is steep (Q15).

## 7. Known quirks and likely bugs

✔ marks items confirmed with wolframscript on 2026-11-10. Ask the user before fixing any open one.

| # | Where | Issue |
| --- | --- | --- |
| Q1 ✔ | `DixonToPercent` | **Fixed.** The `clip=True` branch used \|F\|/\|W+F\| in water-dominant voxels (the noise-biased estimator, the reverse of Liu): FF 2% read 4%. Both branches now use Liu; `clip` only sets the range. |
| Q2 ✔ | MuscleBids `Dixon-B` start B0 | **Fixed.** Divided by n·ΔTE instead of the first-to-last span (17% low with 6 echoes). Now `Last − First` echo time. |
| Q3 ✔ | usage strings | **Fixed.** Units (s, cycles), defaults (`DixonPrecessions` 1, `DixonTolerance` 1, absolute `DixonMaskThreshold` 0.1), `DixonFixT2` values, output structure, `UnwrapThresh` range, `OptimizeDixonEcho` forms. |
| Q4 ✔ | `DixonPhase` unwrap | **Fixed.** The zero-mean DCT solution lost the global phase; `RestoreDC` restores it. |
| Q5 ✔ | `DixonPhase` B0 alias | **Fixed.** B0 is known modulo 1/dt; `SelectAlias` picks the alias per object. This was the main cause of whole-leg swaps (motion study, ses-006 calf). |
| Q6 ✔ | `DixonPhase` φ0 | **Fixed.** The final sign step added Arg to the variable that holds 2φ0 (half correction), and the DC loss removed a constant φ0. Now a global sign on the restored map. |
| Q7 ✔ | `DixonReconstruct` φ0 with complex amplitudes | **Fixed.** Fitting φ0 next to complex amplitudes was singular (PseudoInverse hid it, the φ0 map was arbitrary). φ0 is no longer fitted there; output = Arg(W+F). |
| Q8 ✔ | constrained mode | **Fixed.** Diverged in fat voxels (18% swaps) from the fat-biased start; now starts from the complex-amplitude fit. Relevant for the constrained rows of Froeling 2024 Fig. 7. |
| Q9 ✔ | `DixonFitC` speed | **Fixed.** `PseudoInverse` ran through `MainEvaluate` every iteration. Now an inlined `CholeskySolve`: ~7× faster, same results. A compiled-function call inside it crashed the parallel subkernels; that is why `CholeskySolve` is a plain `Function` (CodeStyle.md). |
| Q10 ✔ | `DixonFitPhase -> True` | **Fixed.** The start-map loop ran over `Length[sel]` instead of all 5 maps, breaking partial phase selections. |
| Q11 ✔ | `DixonFitC` output | **Fixed.** ρ/rms belonged to the φ before the last update. |
| Q12 ✔ | `res` | **Fixed.** Scaled only when the output filter was on. |
| Q13 | output median filter | **Intentional, keep.** The median pulls B0 toward 0 at the mask edge. A background fill (`FillBackground`) was tried and **removed on request**: with the extended BIDS masks the effect is negligible and the fill complicates consistency. Do not re-add. |
| Q14 ✔ | MuscleBids `Dixon-B` | **Fixed.** Did not pass `DixonFieldStrength`, so it always assumed 3 T. Now in `dixOpts`. |
| Q15 | local swaps | **Open (next research step).** Inside one leg where B0 is steep (calf end slices, B0 +200…+300 Hz) swaps remain, in old and new code. Needs a spatial method (region growing / graph cut). Motion study, full rerun: ses-005 left calf still has fibularis brevis/tertius and FHL at ~96% FF and EHL at 50% (old and new alike), while the whole-leg swap of the ses-006 right calf (same four muscles, 57–85%) is gone. |
| Q16 | `DixonPhase` B0 in fat | **Open.** Echoes {2, 4} are not exactly in-phase for fat: B0 starts ≈ 70 Hz off in fat. Causes the 0.2–0.4% isolated swaps. |
| Q17 | input shapes | **Open, existing.** `DixonReconstruct` needs 2D/3D images (mask uses `Closing`); single voxels or `{echo, n}` fail. `DixonPhase` hangs on a single 2D image. |
| Q18 | `SimulateDixonSignal` | **Kept.** `fr` is the water fraction; usage documents it, behaviour kept for existing notebooks. |
| Q19 | "CallDB" coefficients | **Intentional.** Newer than the paper (16.29 + 0.39 db vs 16.32 + 0.38 db). |
| Q20 | BIDS settings | **Intentional.** `DixonFixT2 -> False` and no T1 correction, for stability while T1/T2 values are uncertain. To be evaluated last. |
| Q21 | `DixonPrecessions` | **Checked, not auto-set.** `DixonPrecessionCheck` (after the fit) compares the current and the opposite precession sense on a fixed random sample of 20000 voxels above the 0.9 FF quantile: B0 × R2 grid, water-fat model with complex amplitudes, bipolar phase removed. Ratio < 0.9 "ok", > 1.1 warning (`Message`), else "undecided". The BIDS pipeline logs it and writes it to the json (`DixonPrecessionCheck`). Voxel selection was tested on 4 real stacks (original and conjugated): high-FF voxels separate best (worst 0.59 vs 1.70), water voxels carry nothing (≈ 1.0), "brightest" was marginal in the calf (0.93 vs 1.07). Phantom: 0.08 vs 12.3. Final function on the 4 stacks: hip 0.28/0.33 ok, calf 0.52/0.63 ok; conjugated 1.59–3.59, all warned; ~6 s per stack. |
| Q22 | left-right FF gradient | **Open (research).** The residual left-right FF differences (paper: 0.5–0.64% unconstrained) form a smooth gradient along the readout (left-right). Candidates: odd/even **amplitude** alternation of the bipolar readout (not modelled; Yu 2010 corrects it), coupling of the fitted θ0 with FF. B1 in the T1 correction gives only ±0.16% FF. Chemical-shift displacement (0.31 px per echo, alternating) is an edge effect and does not explain a smooth gradient. Not affected by the Q1–Q12 fixes: motion study ses-005/006, 32 swap-free muscles each, mean Dixon R−L FF −0.54% / −1.00% new vs −0.53% / −0.99% old (left higher in 20/32 and 26/32), while T2-mapping FF of the same muscles gives R−L +0.03% / +0.02% (15/32) and correlates weakly with the Dixon difference. So it is a Dixon bias, not anatomy. Next test: odd/even amplitude map. |
| Q23 ✔ | `Unwrap` family | **Fixed.** Unused Block variables removed. `UnWrapC` merges groups with full-array operations; a union-find version gave identical output but only 7.2 → 6.8 s on two calf stacks (`UnwrapSplit`, 3D), so it was not kept. Used only by `Dixon-B`. |

## 8. Real-data validation

Motion study (`7_motion_long`, author's data), subject 1, sessions 1 and 6, stacks DIXON1 (hip) and DIXON6 (calf),
Dixon-A settings, new code against the existing derivatives:

- ses-006 calf: old output had the **whole left leg swapped**; new output correct (both legs B0 ≈ −150 Hz), except the
  end slice (Q15).
- Hips and ses-001 calf: practically identical to old (median B0 difference 0 Hz, 85–96% of voxels within 10 Hz); the
  rest is the noise rim of the extended mask.
- B0 spans about −100 … +300 Hz across each stack (shim), which is the range where the old code failed.

## 9. Recipes

**Synthetic phantom.** Two discs (legs), subcutaneous fat ring (FF 0.92), muscle (0.03), marrow (0.95); smooth B0 plus
a global offset; φ0 = 0.15 cycles; θ0 left-right gradient; R2 40 s⁻¹; noise 1%; echoes `(1.24 + Range[0, 9])/1000.`.
Signal per echo `msk Exp[-40 t + 2πi(φ0 + B0 t + (−1)^e θ0)] ((1 − FF) + FF fatBasis[t])` with `fatBasis` from
`GenerateAmps["Fixed"]` normalised. Layout for `DixonPhase`/`DixonReconstruct`: `{z, echo, y, x}`. Score: median B0
error and voxels with |FF − true| > 0.4 (swaps).

**Light tests** (64 GB machine, the user often runs parallel kernels; see the memory notes):

- One kernel, `Parallelize -> False`, 2–4 slices. A phantom check takes seconds.
- Anything that goes through `ParallelMap`/`DistributeDefinitions`: also run once with `LaunchKernels[2]` and
  `CloseKernels[]` (Q9 crashed only there).
- Old against new in one kernel: run, then `Unprotect["QMRITools`DixonTools`*"]` and `Get` the HEAD version of the
  file (from `git show HEAD:...`), run again. Without `Unprotect` the reload silently fails.
- Check compiled code with `CompilePrint`: no `MainEvaluate`, no `CompiledFunctionCall`.
- Wrong precession sign = conjugated data (imag → −imag), a free test case on real data.

## 10. Function index

| Function | Line | Role |
| --- | --- | --- |
| `DixonPhase` | L242 | start maps (B0, T2*, φ0, θ0) |
| `RestoreDC`, `SelectAlias` | L366, L370 | constant after unwrap; B0 alias per object |
| `ApplyPhase`, `DotAc` | L396, L405 | compiled helpers of `DixonPhase` |
| `FixDixonFlips`, `FindComplexFlips` | L418 | 180° volume flips (commented out in the BIDS pipeline) |
| `DixonReconstruct` | L480 | IDEAL reconstruction and outputs |
| `DixonFit`, `TotalPhase` | L677, L691 | serial/parallel fit call; Arg(W+F) in cycles |
| `InOutPhase`, `DixonFitC`, `DixonFitFC` | L698, L707, L760 | compiled: in/out phase, IDEAL loop, amplitude refit |
| `GenerateAmps` | L792 | fat model basis from strings or amplitudes |
| `DixonToPercent` | L842 | fractions with magnitude discrimination |
| `DixonPrecessionCheck` | L893 | precession-sense check after the fit (Q21) |
| `SimulateDixonSignal`, `FindInPhaseEchos` | L883, L906 | simulation; echo pair for B0 |
| `Wrap`, `UnwrapList`, `Unwrap`, `UnwrapSplit` | L946– | quality-guided unwrapping (`Dixon-B`) |
| `UnwrapDCT`, `UnwrapDCTi`, `MakeWeights`, `MaskWeights`, `SolvePoisson` | L1283– | least-squares DCT unwrapping |
| `OptimizeDixonEcho` | L1442 | echo-time condition number explorer |
