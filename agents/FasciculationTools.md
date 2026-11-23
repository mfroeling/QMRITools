# FasciculationTools: fasciculation detection (for agents)

Code map of `QMRITools/Kernel/FasciculationTools.wl`: detection of fasciculations as short signal voids in dynamic
MUMRI or DTI data, selection by size, analysis per mask, and the per-event table for curation. Index of all agent
docs: [AGENTS.md](../QMRITools/AGENTS.md). Related: [CodeStyle.md](CodeStyle.md) §11 (the intentional designs below),
[DataConventions.md](DataConventions.md) (array layouts), [MuscleBidsTools.md](MuscleBidsTools.md) §6.3 (the DTI step
that calls it). Line numbers are approximate. Grep for the function name.

The published method is Heskamp et al., *Fasciculation distribution in a healthy population assessed with diffusion
tensor imaging*, Physiol Rep 2025;13:e70247 (doi 10.14814/phy2.70247). Its thresholds 0.6/3.0 (`fr`/`sc`) were
validated with Monte Carlo simulations and expert scoring.

---

## 1. Pipeline

```text
{act, norm, mn, thr} = FindActivations[dat, mask, ActivationThreshold -> {sc, fr}, ..., ActivationOutput -> All]
{sel, heat}          = SelectActivations[act, {regionMask, backMask}, vox, SizeMethod -> "Area"]   (3D mask)
row                  = AnalyzeActivations[sel, regionMask, vox]                                     (summary)
events               = ListActivations[sel, vox, seg, {norm, mn}]                                  (one row per event)
selC                 = RemoveActivations[sel, keys]                                                 (curation)
```

- Arrays: data, `act`, `norm` are `{z, vol, y, x}`, masks and `mn` `{z, y, x}`, `thr` `{z, 3, y, x}`
  (used / SD-based / fraction threshold × mean), `vox = {z, y, x}` mm. `act` and `sel` are 0/1 SparseArrays.
- `FindActivations[data, mask]` normalises each volume to the same mean in the (dilated) mask
  (`NormalizeData[..., NormalizeMethod -> "Volumes"]`, signal ~100), then runs `FindActC` on every voxel time series
  (`RotateDimensionsLeft[Transpose[...]]`, Listable compiled). `IgnoreSlices` pads the skipped slices back with 0.
- **Event**: one 8-connected component (`MorphologicalComponents`) in one slice of one volume. Events are not linked
  across slices or volumes. **Observation**: one slice in one volume.
- `SelectActivations` keeps components of at least `ActivationSize` voxels (`"Voxels"`, default 4) or mm²
  (`"Area"`, default 30, rounded up to whole voxels). With a 4D mask it works per mask.
- Muscle-BIDS (DTI step, MuscleBidsTools.wl ~l.4177-4192): `FindActivations[norm, IgnoreSlices -> {0, 0},
  ActivationThreshold -> {3., .6}, ActivationIterations -> 10, ActivationBackground -> 10]` (the paper values), then
  `SelectActivations[..., ActivationSize -> 4, SizeMethod -> "Voxels"]`. It writes `fasc`, `fascm`, `norm`; there is
  no per-muscle analysis and no event table there.

## 2. The detection algorithm (`FindActC`, ~l.219)

Per voxel time series `t`, with `sc` (number of SDs) and `fr` (fraction of the mean):

1. `Mean[t] <= back` (`ActivationBackground`): background, nothing detected.
2. Zero samples (no data, e.g. out of FOV after registration) are excluded from all statistics and never flagged.
3. Loop (`i <= it`, `ActivationIterations`), `ti` = the non-active samples, at the start all non-zero samples:
   - `mn = Mean[ti]`
   - **SD from the samples above `fr`·mn only**: `ts = samples > fr mn`
   - `tr = Max[0.1, Min[1 - sc SD(ts)/mn, fr]]`
   - new `ti` = samples > `tr`·mn. Stop when its length no longer changes.
   - ≤ 5 samples left in either selection: error state, nothing detected in that voxel.
4. Final: `tr = Min[1 - sc SD(ti)/mn, fr]`, with SD over the non-active samples and **no floor**. A negative `tr` means
   nothing is detected (very noisy voxel). Active = `t <= tr mn`, zeros masked out.

`MeanThresh` (~l.275) recomputes `mn` and `thr` for `ActivationOutput -> All` with the same final rule (non-active,
non-zero samples, `Ramp` only on the shown SD threshold).

Design points from the paper (keep them):

- The two thresholds combine as `Min` on purpose: in variable voxels (blood vessels, low SNR, changing DTI contrast)
  T_SD drops below T_FRAC to prevent false positives. `ThresholdMethod` can switch to one of the two.
- The `fr` pre-selection for the SD keeps the drops out of the spread, so detection holds up to single-voxel rates of
  ~25% at SNR ≥ 10. Above that the method fails (paper Fig. S9).
- The method is conservative: some true events are missed in low-SNR regions (sensitivity ≥ 60% at SNR > 5).

## 3. History: the 4.5.3-5.3 regression (fixed in 5.4)

The rewrite of 2026-01-15 (commit `01f0c423`, "changed fac algorithm", first released in 4.5.3) replaced `Select` with
`Pick` for speed and merged the two selections: the SD was computed on `ti` (all non-active samples) instead of on the
samples above `fr`·mean, and the loop floor went from 0.1 to 0. That looks equivalent but is a fixed-point error: drops
that are not yet detected stay in `ti`, inflate the SD, push the threshold down, and so stay undetected. Found from
Newcastle data sets with many fasciculations (2026-11-19).

Simulation (single voxels, 100 repetitions, 0.6/3.0, Rician noise; sensitivity at a 25% rate):

| SNR | 5.3 (regression) | 5.4 (= 4.3.6) |
| --- | --- | --- |
| 5 | 0.00 | 0.10 |
| 10 | 0.02 | 0.73 |
| 20 | 0.05 | 0.80 |

At low rates both are equal. Precision is ~1.0 at SNR ≥ 10 for both. 5.4 is bit-identical to the 4.3.6 reference on
all simulated cases (zeros excluded) and runs ~1 s per 6×192×192×100 data set (5.3: 0.7 s, it stopped after one
iteration). Fasciculation results made with 4.5.3-5.3 change on a rerun. On whole-leg DTI of a healthy cohort
(Muscle-BIDS settings, 3 subjects × 6 stacks) 5.4 finds 2.7× more events, 1.2-3.9× per stack, with the largest gain in
the stacks with the most fasciculations, so 5.3 maps are not only scaled down but biased between regions.

| Versions | SD in the loop | Zero samples (no data) | Use |
| --- | --- | --- | --- |
| ≤ 4.5.2 (paper: 4.3.6) | published | counted as activations | trustworthy where the data has no zeros |
| 4.5.3-5.2 | regression | counted as activations | not trustworthy |
| 5.3.0 | regression | excluded | not trustworthy |
| ≥ 5.4 | published | excluded | correct |

Rejected alternatives (simulated, see §6): a median/MAD baseline and spread (no better than the fix at SNR ≥ 10, worse
at SNR 5), the fix plus MAD (best at 50% rates, but up to 1.4% false positives in vessel voxels), and the `fr`
pre-selection also in the final step (more sensitive at SNR 5, ~4× more vessel false positives). The user chose the
exact published algorithm without an option.

## 4. Event table and curation

- `ListActivations[act, vox, seg, {norm, mn}]` (~l.484): one association per event, ordered by slice, volume,
  component. Columns: `event`, `slice`, `volume`, `component`, `size [vox]`, `area [mm2]`, `centroid y/x`,
  `bbox ymin/ymax/xmin/xmax` (separate columns, a list cell doesn't export to xlsx); with seg `label` (majority),
  `label fraction`, `obs active [%]` (active seg voxels / seg voxels of that slice); with `{norm, mn}` `signal`,
  `baseline`, `drop [%]`; always `events in obs`, `slices in vol`. Many events in a volume across many slices points
  to motion.
- Counts are the same as `AnalyzeActivations` (`Length` = `"Amount"`, distinct `{slice, volume}` = `"Amount Obs"`).
  Keep it that way when changing either.
- `RemoveActivations[act, keys]` (~l.570), `keys = {{slice, volume, component}, ..}`: relabels each listed
  observation with `MorphologicalComponents` and zeroes the listed components. It only works on the **same** map the
  table was made from, because component numbers follow the scan order of that map. Unknown component numbers are
  ignored, out-of-range slice/volume gives `RemoveActivations::key`. A whole volume: `act[[All, vol]] = 0`.
- Muscle names are not handled here. The caller maps `label` (e.g. `ImportITKLabels["MuscleLabels", "List"]`).

## 5. Quirks

| # | Where | Status | Description |
| --- | --- | --- | --- |
| Q1 | `FindActC` | fixed 5.4 | 4.5.3-5.3 computed the loop SD on all non-active samples (§3). |
| Q2 | `AnalyzeActivationsI` ~l.429 | fixed 5.4 | An observation with exactly one active voxel was skipped, while two separate single voxels counted as two. Only affected maps without `SelectActivations`. |
| Q3 | `FindActC`, `MeanThresh` | intentional | Zero samples are excluded and never flagged (CodeStyle §11). |
| Q4 | `FindActC` | intentional | `Min` of the two thresholds and no final floor: a voxel with SD/mean > (1 - 0)/sc gets no detection. Part of the published method. |
| Q5 | `FindActC` loop | open, minor | Convergence compares the number of non-active samples, not the sets. Same length but a different set stops the loop. Same as 4.3.6. |
| Q6 | `FindActC` | note | `ActivationIterations` runs up to `it + 1` iterations (`i <= it`), as in 4.3.6. |
| Q7 | `ListActivations` | intentional | On a 50/50 label tie the first label in scan order wins (`label fraction` 0.5). |
| Q8 | Muscle-BIDS, `FindActivations` | fixed 5.4 | The defaults were `{3., .65}` (Muscle-BIDS) and `{3.0, 0.75}` (`FindActivations`), now the validated paper values `{3.0, 0.6}`. On MOTION, fr 0.7 instead of 0.6 gives 1.4-2x more events. |

## 6. Recipes

- **Load and test** in a fresh kernel (CodeStyle §12), test scripts in the session scratchpad. Private:
  `QMRITools`FasciculationTools`Private`FindActC`.
- **Simulation** (the method used in §3): per condition 2000 voxel time series of 45 or 100 repetitions, baseline 100,
  Rician noise at SNR 5/10/20; fasciculation voxels with single-sample drops of 30-90% at rates 1-50%; pure-noise
  voxels; vessel voxels `100 (1 + a Sin[φ])`, `a` 0.2-0.5 per voxel, random phase per repetition. Score sensitivity,
  precision and the false-positive % in noise and vessel voxels, for the threshold pairs 0.5/2.25, 0.6/3.0, 0.6/3.5,
  0.75/3.75 (the paper's three plus Newcastle).
  - Guard such scripts: check that the simulated signal and every detection are packed numeric arrays before scoring,
    wrap calls in `MemoryConstrained`, and print tables instead of exporting rows. A list-length mismatch in the vessel
    model once left a symbolic expression that was threaded over all samples and wrote 144 GB to disk.
- **Compiled code**: after editing `FindActC`, check `CompilePrint` for `MainEvaluate` and `CompiledFunctionCall`
  (none now; `=!=` on lengths caused `MainEvaluate`, use `!=`).
- **Diagnosis on real data**: per voxel in the muscle seg, the temporal CV (SD/mean) against the MAD-based CV
  (1.4826 MAD/median). A ratio well above 1 means drops inflate the SD; a high MAD CV means the data itself is noisy
  (SNR = 1/CV). Shoulder data of the Newcastle study had CV 0.22-0.33 (SNR 3-4.5), below the validated range.

## 7. Function index

| Function | ~Line | Role |
| --- | --- | --- |
| `NormalizeFasciculationData` | 137 | normalise DWI data, optionally via the tensor fit |
| `FindActivations` / `FindActivationsI` | 172 | normalise, run `FindActC`, optional `mn`/`thr` output |
| `FindActC` | 219 | compiled detection per time series (§2) |
| `MeanThresh` | 275 | compiled `mn`/`thr` output |
| `SelectActivations` / `SelectActivationI` | 309 | size selection within masks |
| `AnalyzeActivations` / `AnalyzeActivationsI` | 398 | summary per mask: Amount, Amount Obs, chances, sizes |
| `ListActivations` | 484 | one row per event (§4) |
| `RemoveActivations` | 570 | remove events by key (§4) |
| `EvaluateActivation` | 599 | interactive viewer of the detection (`Module`, legacy) |
