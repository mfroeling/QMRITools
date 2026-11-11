# QMRITools data conventions: agent guide

These are the array layouts, axis orders and units that the ~950 public functions assume. Getting one wrong rarely
errors; it gives **silently wrong numbers**. Check this page before passing data between functions, writing a new
fit, or reading or writing NIfTI.

- Index of all agent docs: [AGENTS.md](../QMRITools/AGENTS.md). Related: [ToolboxStructure.md](ToolboxStructure.md)
  (where functions live), [MuscleBidsTools.md](MuscleBidsTools.md) (the pipeline that chains these conventions).
- Sources: the code and usage strings as of 2026-09-22. Where a usage string disagrees with the code, both are
  listed. Trust the code and ask the user before "fixing" either.

---

## 1. Array layout

| Data | Layout | Notes |
| --- | --- | --- |
| 3D volume | `{z, y, x}` = `{slices, rows, columns}` | the first index is the slice |
| 4D volume | **`{z, t, y, x}`**: slices first, then volumes/echoes/diffusion directions/channels | not `{t, z, y, x}`. Many helpers `Transpose` to `{t, z, y, x}` internally and back |
| 2D image | `{y, x}` | |
| Multi-channel NN data | `{z, c, y, x}` | same slot as time. `GetTrainData` picks a random `c` (SegmentationTools.md Q3) |
| Tensor (in memory) | `{6, z, y, x}` with components **`{xx, yy, zz, xy, xz, yz}`** | `TensMat` converts a component vector to the 3×3 matrix |
| Tensor (on disk / nii) | `{z, 6, y, x}` | the BIDS pipeline does `Transpose[tens]` before `ExportNii` and after `ImportNii` |
| Split segmentation | `{z, n, y, x}` masks (0/1, SparseArray) + `labels` list | output of `SplitSegmentations`, input of `MergeSegmentations` |

Useful helpers: `Transpose[data]` swaps the first two levels (`{z,t,…}` ↔ `{t,z,…}`), and `RotateDimensionsLeft` and
`RotateDimensionsRight` move the first axis to the end and back (used before per-voxel fits). `data[[All, 1]]` takes
the first volume of 4D data. This is how `SegmentData` and `ApplySegmentationNetwork` treat 4D input.

## 2. Voxel size

- **`vox` = `{z, y, x}` in mm**, in the same order as the array. `ImportNii` returns
  `Reverse[pixdim[[2 ;; 4]]]` and `ExportNii` writes it back reversed.
- Example: `{6., 1.5, 1.5}` means 6 mm slices with 1.5 mm in-plane. `SegmentData` rescales to exactly this for 3D
  networks.
- Resampling functions take `{voxFrom, voxTo}`: `RescaleData[data, {vox1, vox2}]` and
  `RescaleSegmentation[seg, {vox1, vox2}]`. The output grid changes size, so pad or crop (`PadToDimensions`) when two
  arrays must match.
- The `PlotData` and `PlotData3D` usage strings say vox is "(z,x,y)". That is a naming slip: the value is the same
  `{z, y, x}` vector used everywhere.

## 3. NIfTI I/O

| Function | Returns / does | Defaults to know |
| --- | --- | --- |
| `ImportNii[file]` | `{data, vox}` in the layout above; orientation flips are applied from the header (`NiiFlip`) | **`NiiScaling -> False`**: raw stored values, no slope/intercept. Pass `NiiScaling -> True` for quantitative maps that were stored scaled. `NiiMethod` gives header, TR, rotation, etc. |
| `ExportNii[data, vox, file]` | writes nii / nii.gz | `CompressNii -> True`. The BIDS code writes uncompressed on non-Windows systems and runs `CompressNiiFiles` afterwards |
| `ImportNiiDiff[file]` | `{data, grad, bval, vox}` (bvec/bval with the same base name) | **`FlipBvec -> True`**: gradients are converted from the bvec file with `{1, -1, 1} RotateLeft[g]`. The BIDS pipeline instead uses `FlipBvec -> False` and then `FlipGradientOrientation[grad, flip, perm]` from the config |
| `ExportBval`, `ExportBvec`, `ExportBvalvec` | write the FSL-style text files | |
| `ImportNiiDix`, `ImportNiiT2`, `ImportNiiT1` | vendor-specific corrections for scanner-exported maps | |

Orientation: don't flip axes by hand after `ImportNii`. Use the header-aware import, and fix gradient or tensor
orientation with `FlipGradientOrientation` / `FlipTensorOrientation`, or the `TensorFlips` / `TensorPermutations`
options.

## 4. Diffusion

- **Gradients** `grad`: a list of unit vectors `{{x, y, z}, ..}`, one per volume. The unweighted volume is
  `{0, 0, 0}`.
- **b-values** `bval`: s/mm², one per volume, 0 for unweighted.
- `Bmatrix[bval, grad]` gives the 7-element form `{-bxx, -byy, -bzz, -bxy, -bxz, -byz, 1}`. `Bmatrix[{bval, grad}]`
  gives the 6-element form. `BmatrixConv` converts between them.
- **`TensorCalc`** takes bval in s/mm² and returns a tensor in **mm²/s** (layout `{6, z, y, x}`).
- **`ParameterCalc[tens]`** returns `{l1, l2, l3, MD, FA}` with the eigenvalues and MD in **10⁻³ mm²/s** (so
  muscle MD ≈ 1.5). FA is unitless, 0–1.
- **IVIM** (`IVIMCalc`, `IVIMCorrectData`): f is a fraction 0–1, and D and D* (pD) are in **mm²/s**. The BIDS code
  multiplies ADC by 1000 before export, so the stored `adci` is in 10⁻³ mm²/s.
- `SortDiffusionData`, `ConcatenateDiffusionData` and `SelectBvalueData` keep data, grad and bval in step. Use them
  instead of indexing volumes by hand.
- Gradient-nonlinearity correction (`MakeGradientDerivatives[vox, "WA1"|"WA2"]`) needs the scanner `Offset` that the
  BIDS conversion writes into the json.

## 5. Relaxometry and Dixon units

- **`T2Fit`, `T1Fit` and `EPGT2Fit` return times in the units of the input times.** Pass echo times in ms to get
  T2 in ms. The BIDS pipeline passes `1000 echos` (json `EchoTime` is in seconds), so the stored T2 maps are in ms.
  EPG pulse angles are in degrees.
- **Dixon** (`DixonReconstruct`, `DixonPhase`, `SimulateDixonSignal`):
  - fractions (`watfr`, `fatfr`) are 0–1;
  - B0 is in Hz;
  - the input phase data are in radians (−π…π), the fitted phase maps (`phi0`, bipolar) are in cycles;
  - field strength is in tesla (`DixonFieldStrength`).
- **Echo times for Dixon are in seconds** (the usage strings were corrected to say so in 2026-11).
  - The model multiplies `2π · echo` by fat frequencies in Hz (field × γ × ppm).
  - The initial T2\* is clipped to `{0, 0.25}`, i.e. 250 ms in seconds.
  - The BIDS pipeline passes the json `EchoTime` (seconds) directly and later multiplies `t2star` by 1000 for
    display.

  So T2\* and R2\* come out in s and 1/s, confirmed by the user and a synthetic round-trip. Dixon phases
  (`phi0`, bipolar) are in cycles. `SimulateDixonSignal[echo, fr, ...]` takes the *water* fraction `fr`.
- With complex amplitudes (`DixonConstrainPhase -> False`, the default and the BIDS setting) `phi0` is not fitted, because
  the amplitude phases already carry it. The `phi0` output is the phase of W+F at TE = 0. With real amplitudes
  (`DixonConstrainPhase -> True`) `DixonReconstruct` first fits complex amplitudes and uses that fit as the start.
- `DixonPhase` measures B0 from a pair of in-phase echoes, so B0 is only known modulo 1/dt (500 Hz for the motion
  protocol, echoes 2 and 4). Per connected object it picks the alias whose water-fat model residual is lowest; a wrong
  alias is a water-fat swap. Local swaps inside one object (e.g. the end slices of a calf stack) are not resolved.
- `DixonToPercent` uses magnitude discrimination (Liu, 10.1002/mrm.21301) in both modes: the minor component is
  `1 - major/|W+F|`, so noise can make it slightly negative. `clip` only sets the range (`{0,1}` or `{-0.1,1.1}`).
  For ROI means keep the negatives (`DixonClipFraction -> False`); clipping at 0 biases low fractions upward.
- Philips raw Dixon scaling used in conversion: magnitude `1000 x / 2047`, real/imag `1000 (x − 2047) / 2047`, phase
  `π (x − 2047) / 2047`.
- In the BIDS analysis output (xlsx), `fatfr` is shown ×100 (percent) and Dixon `t2star` ×1000 (ms). The stored maps
  keep the raw units above.

## 6. Masks and segmentations

- **Mask**: numeric 0/1 array with the data's 3D layout. `MaskData[data, mask]` accepts 2D/3D masks for 2D/3D/4D
  data.
- `Mask[data, min]` / `Mask[data, {min, max}]` threshold in **data units**. The pipeline usually normalises first
  (`Mask[NormalizeData[data], 5 | 10 | 15, …]`), so thresholds like 10 refer to the normalised scale (§7).
- **Segmentation**: 3D integer label array, **0 = background**, labels follow an ITK-SNAP label file
  (`ImportITKLabels`, `MuscleLabelToName`, `MuscleNameToLabel`).
- In the muscle label sets, muscles use low numbers and bones sit above a cut-off: `BoneLabel` (default 100) in the
  tractography config, `n` in `analysis.Segmentation.Labels` (for example 90). Bones are taken as the next 30
  labels.
- Networks: training targets are `ClassEncoder[seg, nClass]` with `nClass = max label + 1`. Network output goes
  through `ClassDecoder` back to labels 0..nClass−1. `ReplaceLabels` then maps training labels to output labels by
  muscle name (SegmentationTools.md §9).
- Segmentations are resampled with `RescaleSegmentation` (label-safe), never `RescaleData`.

## 7. Intensity normalisation

- `NormalizeData[data]` (default method `"Set"`) scales to the mean signal. For 4D data it uses the first volume.
  `NormalizeMeanData` does the same on the mean over the 4th dimension.
- `NormalizeMethod -> "Uniform"` maps the histogram to a uniform 0–1 distribution, with 0 treated as background.
  **All network inputs use this** (`ApplySegmentationNetwork`, `GetTrainData`, `MakeTestData`).
- Augmentation assumes data scaled so that roughly Q99 = 1. That is why salt/pepper noise uses 1./0. (CodeStyle.md
  §11).

## 8. Time, files and misc

- TR and TE in json sidecars (dcm2niix) are in **seconds**. Convert explicitly where a function expects ms (§5).
- Tract files: `.trk` with `ExportTracts` / `ImportTracts` (`{tracts, vox, dim, seeds}`). The tract maps
  (`TractDensityMap` and others) take `vox` and `dim` of the target grid. The coordinate frame of the tract points
  has not been checked; read `TractographyTools.wl` before relying on it.
- Check files, logs and BIDS names: see MuscleBidsTools.md §4 and §7.

## 9. Quick sanity checks (typical healthy-muscle values)

- `Dimensions[data]` for a 4D diffusion set should read `{nSlices, nVolumes, nRows, nCols}`, and `Length[bval]` should
  equal `nVolumes`.
- `Length[vox] == 3`, and `Dimensions[data][[{1, -2, -1}]] vox` gives the field of view in mm, z first.
- Muscle DTI after `ParameterCalc`: MD ≈ 1.3–1.7 and FA ≈ 0.15–0.4. Values around 0.0015 mean the 10⁻³ scaling was
  skipped. Values around 1500 mean it was applied twice.
- Muscle T2 from EPG ≈ 25–40 (ms). Values around 0.03 mean the echoes were passed in seconds.
