# SegmentationTools: agent guide

This is a map of [QMRITools/Kernel/SegmentationTools.wl](../QMRITools/Kernel/SegmentationTools.wl) (~3200 lines). The
file holds CNN muscle segmentation: inference, the training pipeline (including parallel batch producers and masked
self-supervised pretraining), training-data preparation, augmentation, evaluation metrics and visualisation. Read this
before changing that file.

- Code state as of 2026-09-22. Line numbers (`L123`) drift over time; grep for the function name if one no longer
  matches.
- House style: [CodeStyle.md](CodeStyle.md). The main consumer of `SegmentData` is the BIDS pipeline:
  [MuscleBidsTools.md](MuscleBidsTools.md) §6.5.
- Network *construction* (`MakeUnet`, `AddLossLayer`, `ChangeNetDimensions`, `NetDimensions`, `ClassEncoder`/`ClassDecoder`/
  `ClassConfidence`, `MakeClassifyImage`, loss layers) lives in **NeuralNetworkTools.wl**. Segmentation label helpers
  (`SplitSegmentations`, `MergeSegmentations`, `ReplaceSegmentations`, `RemoveMaskOverlaps`, `RescaleSegmentation`,
  `GetSegmentationLabels`, `SmoothSegmentation`) live in **MaskingTools.wl**. Cropping and side splitting
  (`AutoCropData`, `ReverseCrop`, `FindCrop`, `ApplyCrop`, `CutData`, `FindMiddle`, `MonitorFunction`) live in
  **GeneralTools.wl**.

---

## 1. The five subsystems

| # | Subsystem | Entry points | Section |
| --- | --- | --- | --- |
| A | **Inference**: anatomy-aware whole-body segmentation | `SegmentData`, `SplitDataForSegmentation`, `ApplySegmentationNetwork`, `ClassifyData` | §3 |
| B | **Patching**: data ↔ patches, merging | `DataToPatches`, `PatchesToData`, `FindPatchDim` | §4 |
| C | **Training**: network training loop and batch producers | `TrainSegmentationNetwork`, `GetTrainData`, `AugmentTrainingData`, `FreezeEncoderLayers` | §5–§7 |
| D | **Data prep and QC**: nii pairs → `.wxf` training sets | `PrepareTrainingData`, `CheckSegmentation`, `Make*Image/Grid`, `ShowTrainLog` | §8 |
| E | **Metrics, labels and external tools** | `DiceSimilarity`, `JaccardSimilarity`, `SurfaceDistance`, `MakeDistanceMap`, `ImportITKLabels`, `MuscleLabelToName`, `RunMuscleMap`, `SegmentDataGUI` | §9 |

Typical user flows:

```wolfram
(* inference *)
seg = SegmentData[{data, vox}, "LegsHip", TargetDevice -> "GPU", SegmentationDimension -> "3D", Monitor -> True];

(* training *)
PrepareTrainingData[{labFol, datFol}, trainFol, InputLabels -> ..., OutputLabels -> ..., TrainVoxelSize -> {6, 1.5, 1.5}];
TrainSegmentationNetwork[{trainFol, outFol}, "Start", PatchSize -> {32, 96, 96}, BatchSize -> 2, RoundLength -> 256];
TrainSegmentationNetwork[{trainFol, outFol}, outFol];                    (* continue from the last saved _itt_ net *)
ShowTrainLog[outFol]
CopyTrainedNetwork[outFol, "UpperLeg", "3D"]                               (* install as a paclet asset *)

(* masked self-supervised pretraining, then transfer learning *)
TrainSegmentationNetwork[{trainFol, preFol}, "Start", MaskedPretraining -> True];
TrainSegmentationNetwork[{trainFol, outFol}, preFol <> "/..._final.wlnet", FreezeEncoderDepth -> 3];
```

---

## 2. Central tables: locations, groups, networks, labels (L398–461)

Everything anatomy-related is driven by three private tables. **Add or change anatomy here, not in code branches.**

```wolfram
$BodyPositionClasses = {"LowerLegs","Knee","UpperLegs","Hip","Torso","Shoulder","HeadNeck"};  (* classifier output order = index *)

$SegmentationLocations = <| loc -> <|"Net2D"->asset, "Net3D"->asset, "TrainLabels"->asset,
                                    "PositionClasses"->{classes that count as this loc}, "Offset"->{startShift, endShift}|> |>
$SegmentationGroups    = <| what -> <|"Locations"->{loc..}, "Split"->"Find"|"Auto",
                                     "Classify"->"Position"|"Side"|"None", "OutputLabels"->asset|> |>
```

| Location | Net2D / Net3D asset (file) | TrainLabels asset | PositionClasses | Offset |
| --- | --- | --- | --- | --- |
| LowerLegs | `SegLegMuscle2D/3D` (N6_LowerLeg_*) | `LegLowerTrainLabels` | LowerLegs, Knee | {0,0} |
| UpperLegs | `SegThighMuscle2D/3D` (N5_UpperLeg_*) | `LegUpperTrainLabels` | Knee, UpperLegs, Hip | {0,0} |
| Hip | `SegHipMuscle2D/3D` (N4_Hip_*) | `HipTrainLabels` | Hip, Torso | {-25,5} |
| Torso | *NotImplemented* | `TorsoTrainLabels` | Hip, Torso, Shoulder | {-5,0} |
| Shoulder | `SegShoulderMuscle2D/3D` (N2_Shoulder_*) | `ShoulderTrainLabels` | Torso, Shoulder, HeadNeck | {-5,0} |
| HeadNeck | *NotImplemented* | `HeadNeckTrainLabels` | Shoulder, HeadNeck | {-5,0} |
| Arm | 3D only `SegArmMuscle3D` (N7_Arm_3D) | `ArmTrainLabels` | – (no position classifier) | – |

`what` groups (the second argument of `SegmentData`, and the `Segment.Location` value in BIDS configs):

| Output label set | Groups (Locations / Split / Classify) |
| --- | --- |
| `MuscleLabels` (full body, S1) | `Body` (LowerLegs, UpperLegs, Hip, Shoulder / Auto / Position), `LegsBody`, `LegsHipBody`, `HipBody`, `UpperLegsBody`, `LowerLegsBody`, `ShoulderBody` |
| `MuscleLegLabels` (legacy legs, S2) | `Legs`, `LegsHip` (Find / Position), `UpperLegs`, `LowerLegs` (Find / Side), `Hip` (Auto / Side) |
| `MuscleShoulderLabels` (S3) | `Shoulder` (Auto / Side) |
| `MuscleArmLabels` (S4) | `Arm` (Auto / None) |

Assets are declared in [PacletInfo.wl](../QMRITools/PacletInfo.wl) (`"Body"` = `Body_Pos_Side.wlnet`, the
position/side classifier). They are resolved by `GetAssetLocation`. The label files are ITK-SNAP label text files.

`GetNeuralNet[name]` (L495) loads a net from a file path or an asset name and **caches it per session** through the
memoised helper `GetNeuralNetI`. `GetNeuralNet["Clear"]` resets the cache. NetGraph and NetChain inputs pass through
unchanged. A failed load is cached as `$Failed` too, so clear the cache after fixing a missing asset.

---

## 3. Inference pipeline (`SegmentData`, L888–988)

```text
SegmentData[{data, vox}, what, opts]            overloads: data | {data,vox} | data,vox ; what default "Body"
  SetMXenvironment["StartSegment"]               MXNet env vars (Reset at the end)
  4D -> data[[All,1]]
  rescale vox -> {6,1.5,1.5} (3D) or {vox1,1.5,1.5} (2D)   ONLY if vox =!= {1,1,1}
  Mask[..,10] -> MaskData -> AutoCropData
  SplitDataForSegmentation[data, what]  -> {{patches, pts, dim}, loc}         (§3.1)
  check Net2D/Net3D exists for every loc, else Return[$Failed]
  per part: ApplySegmentationNetwork[part, netAsset, NetworkOutput->...]      (§3.2)
            ReplaceLabels[seg, {loc, side}, what]                             train labels -> output labels
  PatchesToData[segs, pts, dim, allLabels]      per-label largest component, overlaps resolved (§4)
  ReverseCrop -> RescaleSegmentation back to original dims
  NetworkOutput -> "Both" returns {seg, conf}; anything else returns seg only
```

Key facts:

- **Without `vox` no rescaling happens.** `SegmentData[data]` assumes `data` is already at network resolution. The
  BIDS pipeline always passes `{data, vox}`.
- An invalid `SegmentationDimension` silently falls back to `"2D"`.
- A group containing a location without a net for the chosen dimension (`Torso`, `HeadNeck`, or `Arm` in 2D) returns
  `$Failed`.

### 3.1 `SplitDataForSegmentation` (L1034–1093)

1. The group's `Classify` setting decides how the body is classified:
   - `"Position"` → `ClassifyData[data,"Body"]` → `{side, {{loc, {startSlice, endSlice}}..}}`
   - `"Side"` → side only, and every location covers all slices
   - `"None"` → `"Both"`, all slices
2. If the side is `"Both"`, the data is cut into right/left with `CutData`. `Split "Find"` uses `FindMiddle` (legs).
   `"Auto"` cuts in the middle of the image. The overlap is `SplitOverlap` (0.05) × width.
3. For each side × location, the slices are selected and then `CropPart` autocrops the part. `pts` holds the part's
   ranges in full-volume coordinates.
4. Output: `{{parts, pts, dim}, {{loc, side}..}}`. The `data, seg` overloads split a segmentation identically, for
   making training data.

`ClassifyData` (L587) runs the `"Body"` classifier on 2D slice images (`MakeClassifyImage`). It returns per-slice
`"Side"` and `"Position"` values. The side is the most frequent value. **The classifier always runs on CPU.**
`SplitDataForSegmentation` hard-codes `TargetDevice->"CPU"` for it.

`FindBodyPos` (L619) fits a monotone integer staircase to the noisy per-slice position labels with
`LinearOptimization` (an L1 fit with jumps of 0 or 1, at most one jump per 4 consecutive slices, restricted to the
largest gap-free run of observed labels). Each location's slice range is the union of its `PositionClasses`, shifted
by `Offset`. Legacy class names `"Lower"`/`"Upper"` are mapped to `LowerLegs`/`UpperLegs`.

### 3.2 `ApplySegmentationNetwork` (L1121, L1152)

- Single dataset: accepts an array, a nii path, or a 4D array (takes `[[All,1]]`). It autocrops, pads by
  `DataPadding`, then calls `FindPatchDim` to pick the **largest patch that fits `MaxMemorySize`** (default 32 GB on
  CPU, 8 GB on GPU). `ChangeNetDimensions` resizes the net, `DataToPatches` cuts the data, each patch is normalised
  (`NormalizeData[.., "Uniform"]`), and the patches go through the net (`WorkingPrecision` is `"Mixed"` on a Windows
  GPU, otherwise `"Real32"`).
- `NetworkOutput`:
  - `"Segmentation"` → `ClassDecoder`, merged per label
  - `"Confidence"` → `ClassConfidence`, averaged
  - `"Both"`
  - `"Volume"` → a raw averaged volume (for pretrained reconstruction nets)
- The third argument `node` returns the activation of a named layer for the **first patch only** (for debugging and
  feature inspection).
- Folder mode `ApplySegmentationNetwork[{datFol,outFol} | {{datFol,outFol},{inTag,outTag}} | {{..},startIndex}, net]`
  processes `*<inTag>.nii.gz` (default tags `data` → `label_NN`) and also writes a `.png` grid per file.

`FindPatchDim` (L1241): patch sizes are multiples of the net's downsampling factor
(`Input / MinEncodingOut`). A 2D net uses the full slice. For a 3D net, when over budget, it shrinks z first (down to
about ¼ of the in-plane size), then the larger in-plane dimension (down to 75 % of the ratio), then the other. The
reported memory is divided by 4 (an empirical factor).

---

## 4. Patching (L699–831)

- `GetPatchRangeI`: the number of patches per dimension is `Ceiling[(dim-2pad)/(patch-2pad)] + PatchNumber`, spread
  evenly. If `dim <= patch`, there is one range `{1,dim}` and the patch is **zero-padded on the right** by `GetPatch`.
- `DataToPatches[dat, patch, nPatch|"All", PatchNumber->, PatchPadding->]` returns `{patches, ranges}`. Passing
  `ranges` instead extracts at known positions. That is how the segmentation is cut identically to the data.
- `PatchesToData[patches, ranges, dim]` **averages** overlaps (SparseArray sums divided by unitized counts) and clips
  the right-side padding.
- `PatchesToData[patches, ranges, dim, labels]` works **per label**: it merges each label's binary masks, keeps only
  the largest connected component (`SmoothMask[MaskComponents->1]`), resolves overlaps (`RemoveMaskOverlaps`), and
  calls `MergeSegmentations`. **This is where "one cluster per muscle" is enforced.**

---

## 5. Training: `TrainSegmentationNetwork` (L1353–1671)

### 5.1 Sequence

1. `SetMXenvironment["StartTrain"]`. Options are split into `MakeUnet` options (`netOpts`, via `FilterRules`) and
   training options. All option values are written to `<name>_settings.txt`.
2. Training files: `*.wxf` in `inFol` (a string or a list of folders). Each file is `{data, seg, vox}` as written by
   `PrepareTrainingData`.
3. **Test set**: made by `MakeTestData` from the **first file** (2D: about 9 slices; 3D: a slab of 2×patch depth) and
   saved as `<name>_testSet.nii`. It is reused if its dimensions match the patch.
4. **Network properties from the first file** (intentional: all training files must hold every label, as guaranteed by
   `PrepareTrainingData`): `nClass = Max[label]+1` (or 1 when pretraining).
   `nChan = 1` always. `MultiChannel` has no effect yet; 4D data trains on a random channel per sample (Q3).
5. `netCont` can be:
   - `"Start"` → `MakeUnet[nChan, nClass, patch, netOpts]`
   - a `NetGraph`
   - a `.wlnet` file
   - a previous output folder: the last `*_itt_NNNN.wlnet` is loaded and `ittTrain` is parsed from its name
   Fewer than 5 remaining rounds → `::itt`. `MaxTrainingRounds` is an **absolute** round count.
6. Loss list:
   - `All` = `{Dice, MSD, Tversky, CE, Jaccard, Focal}`
   - pretraining = `{MSDM, MAEM}` (masked losses)
   - otherwise validated against `{Dice, MSD, MAE, MSDM, MAEM, Tversky, CE, Jaccard, Focal, TopK}`
7. `ChangeNetDimensions` (patch, channels, classes), then `NetInitialize[Kaiming]`. Only *uninitialised* arrays are
   initialised, so continued weights are kept.
8. `FreezeEncoderDepth` (only when `netCont =!= "Start"`) → `FreezeEncoderLayers[net, n, chanIn === nChan]` → the
   `"start"` + `enc_1..enc_n` learning-rate multipliers are set to 0. The `"start"` layer is left trainable if the
   channel count changed. `AddLossLayer` nests the net under `"net"`, so the multiplier keys are prefixed with
   `{"net", …}`.
9. `monitorFunction` runs every `MonitorInterval` rounds. It increments `ittTrain`, runs the test set on CPU, exports
   `_itt_NNNN.nii/.png/.wlnet`, and deletes the `.wlnet` from 2 intervals earlier.
10. `batchFunction[n]` = `GetTrainData[data, n, {patch, nClass}, PatchesPerSet, AugmentData, PadData, MaskedPretraining]`.
11. `OneCycleSchedule[roundLength/batch, rounds, ittTrain]` gives a learning-rate multiplier: a cosine warm-up from
    0.2→1 over 15 %, a plateau to 50 %, a cosine decay to 0.1 at 95 %, then 0.1. By default a continued run resumes
    mid-cycle. `RestartLearningCycle->True` fits a new cycle to the remaining rounds.
12. Validation set: `Min[50, 0.2 RoundLength]` samples, saved as `<name>_validation.wxf`. It is reused only when
    continuing (`ittTrain>0`) and the dimensions match.
13. `NetTrain[net, {batchSource, "RoundLength"->…}, All, ADAM(L2, schedule, β 0.9/0.998, ε 1e-4, clip 1),
    WorkingPrecision "Mixed", TrainingProgressReporting -> <ISO-date>.json]`.
14. Export `<name>_trained.wxf` (the full NetTrain result), `_final.wlnet`, and `_final.onnx`.

Output folder content (`<name>` = the output folder's base name): `<name>_settings.txt`, `_testSet.nii`,
`_validation.wxf`, `_itt_NNNN.{nii,png,wlnet}` (only the last 2 wlnets are kept), `<yyyymmddThhmmss>.json` logs (read
by `ShowTrainLog`), `_trained.wxf`, `_final.wlnet`, `_final.onnx`.

### 5.2 Batch sources: `UseParallelKernels`

| Value | Path | How batches reach NetTrain |
| --- | --- | --- |
| `False` | in-process | NetTrain calls `batchFunction` directly (serial augmentation) |
| `True` (**default**) or `{True, n}` | **raw WSTP links** (`LinkLaunch`) | `LaunchTrainingKernels` starts `$ProcessorCount` (or n) kernels, each loads `QMRIToolsDev``/`QMRITools`` + `SetMXenvironment`. `LoadTrainData` sends each link its file share (`PartitionProducerFiles`: files repeated so each link gets ≥ 3, randomly partitioned) and defines `batchFunctionL[n]` remotely with patch/opts spliced in through `With`. NetTrain's source is `GetFromBatchQueueL`: round-robin `LinkReadyQ` polling that re-requests a batch right after reading one |
| `"Parallel"` or `{"Parallel", n}` | legacy Parallel framework | `CloseKernels[]`, `LaunchKernels[n+1]`, `SetSharedVariable` queues `queue1..n` plus `produced/used/ready/…`. `GetKernels` submits n producer loops and 1 trainer kernel that runs `trainFunc` and polls `ready[[i]]` |

Robustness rules shared by both parallel paths:

- A producer returning `Null` 10 times closes that link (`::nullbatch`). `$Failed` marks it dead immediately.
- 20 s without any batch → `Abort[]`.
- The link path wraps training in `CheckAbort` and always `LinkClose`s.
- There is no relauncher for dead producers (a deliberate decision).
- **The `index`, `deadL`, `nullCount`, `used`, `produced`, `activeProducers` counters are `TrainSegmentationNetwork`
  `Block` locals** that `GetFromBatchQueueL` (and the `producerStatus` display) read and mutate through **dynamic
  scoping**. The same holds for `monitorFunction`, which mutates `ittTrain`, `im`, `testSeg` and others. This works
  *because* the code uses `Block`. Do not convert these to `Module` or move the helpers out of the call chain without
  passing the state explicitly.
- Anything sent with `LinkWrite[link, Unevaluated[...]]` must have its master values spliced with `With` (CodeStyle.md
  §5). On the remote kernel, `data` and `batchFunctionL` are intentionally global.
- Machine gotcha: set the front-end option "Launch parallel kernels" to *not* "At startup". Otherwise every linked
  kernel spawns its own pool.

---

## 6. `GetTrainData` (L2227–2304) and patch sampling

```text
GetTrainData[dataSets, nBatch, {patch, nClass}, PatchesPerSet, AugmentData, PadData, MaskedPretraining]
  repeat until nBatch samples:
    pick a random dataset: a .wxf path (Import) | {dat.nii, seg.nii} | in-memory {dat, seg, vox}
    4D data -> ONE random channel (multichannel is not really supported, see Q3)
    AugmentTrainingData (on the whole volume)
    PatchTrainingData -> PatchesPerSet random patches (PatchNumber->2 overlap grid; 2D: full-slice patch, random slices)
    drop all-zero patches
  random subsample to nBatch; optional AddPadding (PadData: zero a random 0..p top/bottom slab in 30% of samples)
  NormalizeData "Uniform" per sample
  segmentation mode: target = ClassEncoder[seg, nClass] (or seg+1 if nClass False), Byte NumericArray
  pretraining mode:  mask = MakeBlockMask; input = data(1-mask) + mask*(noise or 0, 50/50); target = clean*mask
  returns {NumericArray[{input}] -> NumericArray[target] ..}   (input gets a leading channel axis of 1)
```

`AugmentData`: `True`/`False`, `"2D"` (in-plane only: no out-of-plane rotation and no z-scaling), or `"3D"`.

---

## 7. Augmentation (`AugmentTrainingData`, L2004–2189)

The switches are `AugmentationDefaults` = `Flip, Rotate, Scale, Noise, Blur, Bias` (all True). Pass `True`/`False`
or a partial association. Each step fires with probability 0.5 (`Coin[]`):

| Step | What it does |
| --- | --- |
| Flip | `ReverseC`: reverses the **last axis only** (left-right) |
| Rotate | z ±30°. y/x ±15° (not in 2D mode) |
| Scale | per axis 0.6–1.6 (z not in 2D mode). Transforms use `DataTransformation` (data: order 1, seg: order 0) followed by a crop |
| Blur | Gaussian blur σ 0.5–2, or an unsharp-mask sharpen of the foreground (50/50) |
| Bias | a smooth multiplicative field |
| Noise | `AddSaltAndRice`: Rician noise at SNR 5–50, plus salt/pepper in 50 % of cases (`SaltAndRiceC` writes 1./0.) |

- In `"MaskedPretraining"` mode it returns `{corrupted, clean}`, where *clean* is the data **after the geometric steps
  but before blur/bias/noise**.
- Data-only calls (`seg === 1`) return only the data.
- `MakeBlockMask[dim, n]` / `makeBoxC`: randomly rotated boxes plus scattered noise points until the target count
  (`2·RandomReal[{.1,.5}]·voxels`, before deduplication) is reached.

**Intentional designs (do not "fix")**: salt/pepper at 1./0. (data is pre-scaled so Q99 = 1); last-axis-only flip; wide 0.6–1.6 scale range; the `2x` factor
in `MakeBlockMask`; target = `clean·mask` without dilation (dilation was tried and reverted); mask fill applied to the
background too. Details are in CodeStyle.md §11.

---

## 8. Training data preparation and QC

`PrepareTrainingData[{labFol, datFol}, outFol]` (L2376):

1. Finds `*<LabelTag>.nii.gz` files and matches each to a `*<DataTag>.nii.gz` file by replacing the tag in the name.
2. Checks that voxel sizes and dimensions agree.
3. `PrepTrainData`: rescales to `TrainVoxelSize`, crops to the dilated body mask, and remaps labels
   `InputLabels → OutputLabels` (Automatic = keep).
4. `SelectTrainData`: keeps the slice range where ≥ `SegmentationsPerSlice` labels are present, ±8 slices.
5. `CheckSegmentation` flags labels with more than one component or with holes. The colour legend is red = both,
   purple = n>1, blue = hole.
6. Optional `SmoothSegmentation` cleanup (`CleanUpSegmentations`).
7. Writes `name[_OutputTag]_data.nii`, `_label.nii`, `.png`, and `.wxf` (`{Real32 data, Integer16 seg, vox}`), plus
   `summary.png`. `TestRun->True` only analyses.

Visualisation (L2541–2660): `MakeChannelImage` (middle slice, 1–99 % window), `MakeClassImage` (RomaO colours with
even/odd interleaving, middle slice), `MakeChannelClassImage` (40 % overlay), and `MakeChannelClassGrid` /
`MakeChannelGrid` (an n×n or {n,m} grid of slices through the shared `GridLayout`). These grids produce the training
monitor PNGs and the BIDS `_grid` images.

`ShowTrainLog[folder, minEntries]` (L2881): a Manipulate over all `*.json` NetTrain logs in the folder, with
`LoadLog` stitching the runs together. Gridlines mark run boundaries. It offers loss/metric toggles, log scale,
Gaussian smoothing, browse/reload, and export to `TrainLogPlot.png`.

---

## 9. Metrics, labels, external tools

- `DiceSimilarity` / `JaccardSimilarity` (compiled): per class, with **+1 smoothing** in numerator and denominator,
  so two empty segmentations score 1.
- `SurfaceDistance` gives symmetric surface distances between `GetEdge` perimeters (scaled by vox). `Method` is
  `Mean | Median | RMS | HD | HD95` (default) `| Std`, or a list of them. It returns `"noSeg"` if either
  segmentation is empty.
- `MakeDistanceMap[mask, vox, DistanceRange->Automatic|All|0|n]` returns a signed distance to the mask edge:
  positive inside the mask, negative outside.
- `ImportITKLabels[file|asset, "Labels"|"Names"|"List"]` parses ITK-SNAP label files. Names become
  `Capitalized_Words_Side`. `MuscleLabelToName` and `MuscleNameToLabel` build on it (default asset `MuscleLegLabels`).
- `ReplaceLabels` (L995) maps network labels to output labels **by name**: train label → name via the location's
  `TrainLabels` → the number of `name` or `name_<side>` in the group's `OutputLabels`. **Label names in the train and
  output label files must correspond.**
- `CopyTrainedNetwork[fileOrFolder, loc, dim]` copies the newest `.wlnet` into
  `QMRITools/NeuralNetworks/N#_<loc>_<dim>.wlnet`. `loc` uses the file-name vocabulary (`UpperLeg`, `LowerLeg`,
  `Arm`, `PosSide`, …), which is **not** the same as the `$SegmentationLocations` keys. The `dim` argument is required
  for segmentation nets.
- `RunMuscleMap[file | {data, vox}]` shells out to an external MuscleMap conda environment (auto-detected in
  `~/.conda|miniconda3|anaconda3/envs/MuscleMap`; the script is found through `pip show scripts`). It works in
  `$TemporaryDirectory/QMRIToolsMM` and writes a log. Labels are converted to the QMRITools `MuscleLabels` through
  name matching (`ImportMMLabels` alias table) unless `MuscleMapLabels->True`. It returns `{seg, json}`.
- `SegmentDataGUI[]` is a simple dialog (Legs/UpperLegs/LowerLegs/Shoulder, CPU only). See Q2.

---

## 10. Known quirks and likely bugs

✔ marks items confirmed with wolframscript on 2026-09-22. Ask the user before fixing any of them.

| # | Where | Issue |
| --- | --- | --- |
| Q1 ✔ | `MakeDistanceMap` usage | **Fixed 2026-09-22 (usage only).** The code gives positive inside and negative outside (9³ cube test: centre +2, outside −2). The usage now says so. |
| Q2 ✔ | `SegmentDataGUI` L3046 | **Fixed 2026-09-22.** It copied the non-existent asset `"MusclesLegLabels"`; it now copies `$SegmentationGroups[what, "OutputLabels"]` (leg or shoulder labels). **The GUI is work in progress** (a planned project): 4 regions only, CPU only, no vox passed to `SegmentData`. Don't polish it piecemeal; ask the user first. |
| Q3 ✔ | `TrainSegmentationNetwork`, `GetTrainData` L2263 | **Resolved 2026-09-22 (option 1):** the broken detection was replaced by `nChan = 1`, `MultiChannel` is kept but its usage says it has no effect yet, and real multichannel support was deferred. Background: `MultiChannel->True` was a **silent no-op**. `depth = ArrayDepth@testDataRaw` is taken on the ragged `{dat, seg, vox}`, so it is always 1 and `nChan` is always 1. If it did fire, `Length@First@testDataRaw` would give the slice count, not the channel count (the 4D layout is `{z, c, y, x}`). The rest of the chain is single-channel anyway: `GetTrainData` picks one random channel per sample, `MakeTestData` and `ApplySegmentationNetwork`/`SegmentData` take channel 1, and the input is wrapped as `{#}`. The real behaviour for 4D data is **random-channel training (channel augmentation) with first-channel inference**. |
| Q4 ✔ | `GetTrainData` overloads L2227–2231 | **Fixed 2026-09-22.** A direct call with a 2D patch and no nClass used to be misread as patch = 96, nClass = 96. The core overload now requires `{patch:{__Integer}, nClass_}`. Keep the patch integer-typed when editing. |
| Q5 | `TrainSegmentationNetwork` L1433 | **Intentional, not a bug.** `nClass` (and the test set) come from the **first file**. Training files are expected to be consistent: every file contains all labels, which `PrepareTrainingData` ensures and checks. Do not add a scan over all files. |
| Q6 | `SegmentData` / `SplitDataForSegmentation` usage | **Fixed 2026-09-22.** Both now list all 14 `what` groups. The defaults are `"Body"` (SegmentData) and `"Legs"` (SplitDataForSegmentation), as in the code. The non-existent `{what, netFile}` form was removed. Keep both lists in sync with `$SegmentationGroups`. |
| Q7 | `FindBodyPos` L655 comment | **Fixed 2026-09-22.** The comment said "within 6 slices"; it now says "at most one jump within 4 slices", matching `Total[dVars[[i ;; i + 3]]] <= 1`. |
| Q8 | leaked private globals | **Fixed 2026-09-22:** `GetNetwork`, `makeTest`, `datC`, `plot` are now Block locals. `ti` is gone: `OneCycleSchedule` now returns `OneCycleValue[#1 + it, n]&` (checked identical on every batch). **Intentional globals (keep):** `segmentWindow` (lets the GUI close its previous window), `GetNeuralNetI` (session cache), `data`/`batchFunctionL` on producer link kernels, `data` and `queue1..n` on Parallel kernels (must persist between calls). |
| Q9 | `SegmentData` L945 | `Return[$Failed]` for a missing net happens before `SetMXenvironment["Reset"]`, so the MXNet environment variables stay on the "StartSegment" settings. |
| Q10 | `SplitDataForSegmentation` L1054 | `TargetDevice` is ignored for the classifier (hard-coded CPU). This may be intentional (a small net). |
| Q11 | `MakeTrainData` / `NormDat` L2672 | Private, unused anywhere. `NormDat` also uses Block initialisers (old style). |
| Q12 | `CopyTrainedNetwork` without `dim` | Produces `N5_UpperLeg_.wlnet` (trailing underscore), which matches no asset. |
| Q13 | `GetNeuralNet` | A failed load is memoised as `$Failed` until `GetNeuralNet["Clear"]`. |

---

## 11. Recipes

**Add a body location or a new network:**

1. Prepare the data (`PrepareTrainingData` with `OutputLabels` in a new `<Loc>TrainLabels` label file) and train.
2. `CopyTrainedNetwork[outFol, "<FileLoc>", "3D"]`. If the location is new, extend the `Switch` in
   `CopyTrainedNetwork` with the next `N#_` prefix.
3. Register the net and label assets in `PacletInfo.wl`.
4. Add an entry to `$SegmentationLocations` (`Net2D`/`Net3D`, `TrainLabels`, `PositionClasses`, `Offset`) and add the
   location to the relevant `$SegmentationGroups`, or create a group. No code branches are needed: `SegmentData`,
   `SplitDataForSegmentation`, `FindBodyPos` and `ReplaceLabels` are all table-driven.
5. Make sure the train label names exist, with or without `_Left`/`_Right`, in the group's `OutputLabels` file.
6. If it needs a new *position class*, the `"Body"` classifier must be retrained and `$BodyPositionClasses` updated
   (order = classifier index).
7. Update the `SegmentData` usage string and MuscleBidsTools.md (`Segment.Location`).

**Add an augmentation:** add a key to `AugmentationDefaults` and a `Coin[]`-gated block in `AugmentTrainingDataI`
that keeps data and seg consistent (geometric steps must transform both; intensity steps touch data only and come
after `datC` is captured for pretraining).

**Add a loss:** implement it in NeuralNetworkTools (`AddLossLayer`) and add the name to the validation list in
`TrainSegmentationNetwork` (L1472) and to the `::loss` message and usage.

**Debug training throughput:** use `UseParallelKernels -> {True, 2}` for a small-scale run, watch the
producer-status grid (Produced/Used per link), and check `ReadLinkL` timeouts (120 s at launch).

---

## 12. Function index

| Function | Line | Role |
| --- | --- | --- |
| `$BodyPositionClasses`, `$SegmentationLocations`, `$SegmentationGroups` | L398–461 | anatomy tables |
| `CopyTrainedNetwork` | L468 | install a trained net |
| `GetNeuralNet` / `GetNeuralNetI` / `NeuralNetFunc` | L495–510 | cached net loading |
| `ImportITKLabels`, `MuscleLabelToName`, `MuscleNameToLabel` | L523–568 | label files |
| `ClassifyData`, `FindBodyPos` | L587, L619 | body position/side |
| `PatchesToData`, `PatchesToDataI` | L704, L743 | patches → volume |
| `DataToPatches`, `GetPatch`, `GetPatchRanges`, `GetPatchRangeI` | L779–831 | volume → patches |
| `SetMXenvironment` | L842 | MXNet env presets `StartSegment`/`StartTrain`/`Reset` |
| `SegmentData` | L888 | inference entry |
| `ReplaceLabels` | L995 | train → output labels |
| `SplitDataForSegmentation`, `CropPart` | L1034, L1100 | anatomical splitting |
| `ApplySegmentationNetwork` | L1121, L1152 | run a net on data or a folder |
| `FindPatchDim` | L1241 | memory-bounded patch size |
| `TrainSegmentationNetwork` | L1353 | training entry |
| `MakeTestData` | L1678 | monitor test slab |
| `FreezeEncoderLayers` | L1723 | transfer-learning LR multipliers |
| `OneCycleSchedule`, `OneCycleValue` | L1740 | LR schedule closure + phase function |
| `LaunchTrainingKernels`, `ReadLinkL` | L1759, L1820 | producer kernel start / link read |
| `LoadTrainData`, `PartitionProducerFiles` | L1843, L1906 | data distribution + validation |
| `GetFromBatchQueueL` | L1919 | link batch source |
| `GetKernels` | L1959 | legacy Parallel producers + trainer |
| `AugmentTrainingData(I)`, `Coin`, `CoinN`, `ReverseC`, `AddSaltAndRice`, `SaltAndRiceC` | L2004–2150 | augmentation |
| `MakeBlockMask`, `makeBoxC` | L2157 | pretraining masks |
| `AugmentImageData` | L2196 | 2D image augmentation (classifier training) |
| `GetTrainData`, `AddPadding`, `PatchTrainingData` | L2227–2350 | batch generation |
| `PrepareTrainingData`, `SelectTrainData`, `PrepTrainData`, `CheckSegmentation` | L2376–2533 | data prep and QC |
| `GridLayout`, `MakeChannelClassGrid`, `MakeChannelGrid`, `MakeChannelClassImage`, `MakeClassImage`, `MakeChannelImage` | L2547–2660 | images |
| `NormDat`, `MakeTrainData` | L2672 | unused |
| `DiceSimilarity(C)`, `JaccardSimilarity(C)`, `SurfaceDistance`, `SufDistFunc`, `GetEdge` | L2695–2819 | metrics |
| `MakeDistanceMap`, `DistFun` | L2833, L2868 | signed distance |
| `ShowTrainLog`, `LoadLog` | L2881, L2982 | training log viewer |
| `SegmentDataGUI` | L3007 | dialog |
| `RunMuscleMap`, `FindMuscleMap`, `FindMuscleMapEnv`, `ImportMMLabels` | L3079–3210 | external MuscleMap |
