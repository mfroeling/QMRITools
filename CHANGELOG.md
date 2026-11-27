# QMRITools major versions

This page summarises what each major version of QMRITools brought, with a short entry for each minor release since
5.0. The full release notes and the paclet files are on the
[releases page](https://github.com/mfroeling/QMRITools/releases).

| Version | Released | Requires | Public functions and options |
| --- | --- | --- | --- |
| 5.0 | September 2026 | Wolfram Language 15.0+ | ~1000 |
| 4.0 | August 2024 | Wolfram Language 14.0+ | ~875 |
| 3.0 | December 2022 | Wolfram Language 13.0+ | ~750 |
| 2.0 | January 2019 | Mathematica 11.0+ | ~480 |
| 1.x | 2016–2018 (as DTITools) | Mathematica 11.3 | ~440 |

------------------------------------------------------------------------

## QMRITools 5.x updates

- **5.5** (unreleased): Muscle-BIDS `Segment.Method` takes a list, so several segmentations of one map are made next
  to each other, each with its own name and check file. New method `"MuscleMap"` runs `RunMuscleMap` and saves
  `seg_musclemap_…` with the labels in QMRITools numbering. New `Segment` key `MuscleMapEnvironment`. `"Automatic"`
  may now be written as text in the config, and an unknown method is reported instead of skipped. The
  `"Registration"` method registers the whole target once (it registered three single slices before) and has its own
  check file `seg_reg_…`. **Existing Registration segmentations are redone once.** Convert skips a scan whose
  magnitude and phase were put in one file (`_phMag`, classic DICOM with `MergeEchos` true) and says how to fix it;
  before, that file was converted as phase only. The Dixon method `Dixon-A` now also reads Siemens data (scanner
  water, fat, fat fraction, T2\* and goodness of fit maps, and the magnitude echoes) from old and new dcm2niix, and is
  the advised method for all Dixon data; the other Dixon methods stay for older data. Repeated Dixon source scans are
  no longer stacked into one series: the last series per echo time is used. **Dixon-A source data with repeated scans
  (more echoes than acquired) changes on a rerun.**
- **5.4** (unreleased): `FindActivations` uses the published fasciculation detection again (Heskamp et al. 2025,
  QMRITools 4.3.6). Since 4.5.3 the standard deviation also included undetected drops, so fewer fasciculations were
  found, most of all in muscles with many. **Fasciculation counts from 4.5.3-5.3.0 change on a rerun.** New
  `ListActivations` lists every fasciculation with its slice, volume, size, position, label and signal drop, so
  events can be checked, for example for motion, and `RemoveActivations` removes the rejected ones from the map.
  `AnalyzeActivations` now also counts a slice with a single active voxel (only maps without `SelectActivations`).
  The default `ActivationThreshold` of `FindActivations` and of the Muscle-BIDS fasciculation detection is now the
  published `{3.0, 0.6}` (was `{3.0, 0.75}` and `{3., .65}`). Muscle-BIDS merge, segmentation, tractography and
  analysis compare the dates in the check files. They redo a subject when an input was finished after their own
  result, and skip it (with a log line and an error check) when the input itself is out of date, for example after a
  failed merge. Before, they could run on old merged data. Dixon: the bipolar phase is now set per leg from the
  echoes themselves, and the graph cut keeps the start field in water. This fixes the remaining calf water-fat
  swaps (patches or a whole leg) and a +0.5 pp left-right calf fat fraction bias of 5.3.0. **Dixon fat fractions
  from 5.3.0 change on a rerun, mostly in the calf.** `SegmentData` checks that the networks load before it starts,
  and a failed network load is no longer cached. New `Process` key `DixonUnwrapSplit` for `Dixon-B` on a single limb.
- **5.3** (November 2026): Dixon water-fat swaps fixed (alias selection, graph cut start field) and a precession
  check. EPG T2 models fat as 90% lipid and 10% water, pools B1 over a 3D box and corrects for short TR. Every
  Muscle-BIDS output json records the settings it was made with, and the analysis can weight tract-based values
  by tract density. Settings that were hard-coded are now config keys. **The EPG fat fraction is now the lipid
  fraction, about 10% lower than before. `Tractography.BoneLabel` defaults to `"Body"`; configs for older leg models
  need `"Legs"`.**
- **5.2** (November 2026): one debug switch for all toolboxes (`$Debug`, `DebugPrint`), messages instead of prints,
  safer conversion of multi-echo spin echo DICOM in Muscle-BIDS, and a faster config and check GUI.
- **5.1** (October–November 2026): `MuscleBidsGUI` and check tools that follow the config file, Rician bias
  correction for denoising and DTI, decompression of compressed DICOM, and segmentation comparison tools.

------------------------------------------------------------------------

## QMRITools 5.0

About two years of development since 4.0. It adds three new toolboxes, makes deep learning segmentation a complete
train-and-deploy workflow, and turns Muscle-BIDS into an end-to-end pipeline from DICOM to per-muscle results.

### New toolboxes

- **NeuralNetworkTools**: general building blocks for neural networks, split out of SegmentationTools. It includes
  loss layers (cross-entropy, top-K, overlap losses), configurable activation layers, freezing encoder layers for
  transfer learning, copying trained weights between networks, and inspecting network nodes and features.
- **ShapeTools**: statistical shape models and meshes of segmented structures. You can build a shape model
  (`MakeShapeModel`), fit it to new data (`FitShapeModel`, `ApplyShapeModel`), and evaluate and plot it
  (`EvaluateModel`, `PlotShapeVariation`). It also builds and splits region meshes (`MakeRegionMesh`, `SplitRegionMesh`), makes muscle
  templates (`MakeMuscleTemplate`, `TemplateToVolume`), and plots meshes (`PlotMesh`, `MeshGridPlot`).
- **AmaresTools**: AMARES spectral fitting for MR spectroscopy, with basis construction, starting values, analytic
  Jacobians, and a parameter matrix you can edit to constrain the fit.

### Muscle segmentation

- New network architectures: UNET++ and strided-convolution downscaling.
- Training on several contrasts at once, and a better data preparation step (`MakeTrainData`).
- Masked pretraining (`MaskedPretraining`) with encoder freezing (`FreezeEncoderLayers`, `FreezeEncoderDepth`), for
  transfer learning to new contrasts or anatomies with little labelled data.
- Faster training: data loading and augmentation run in parallel background kernels (`UseParallelKernels`), and
  learning-rate restart cycles are supported (`RestartLearningCycle`).
- More anatomy, including shoulder segmentation. Anatomy groups and label handling were reworked, with configurable
  output labels, label checks (`CheckConfigLabels`), and merging of several segmentations (`JoinSegmentations`).
- Per-class confidence maps (`ClassConfidence`) and segmentation QC plots
  (`SegmentationCrossSection`, `SegmentationsPerSlice`).
- Distance-map-based mask tools (`MaskToDistanceMap`, `MaskFromDistanceMap`), `MaskVolume`, and more robust masking
  of data with very bright pixels.

### Muscle-BIDS

- `MuscleBidsAnalysis`: per-muscle extraction of quantitative values from processed data and segmentations.
- Tractography inside the BIDS pipeline (`BidsTractographyMethod`).
- T2 mapping in `MuscleBidsProcess`, plus better handling of multi-stack Dixon and DTI.
- Much better DICOM → BIDS conversion, including multiple sessions, repeated scans and protocol naming
  (`GetProtocolNames`, `ViewProtocolNames`).
- Configs can be split over several files and combined (`MergeConfig`).
- Command-line and batch use (`ParseCommandLine`, `ProcessSubjects`).
- Integration with MuscleMap (`RunMuscleMap`).

### Diffusion, IVIM and gradients

- Random permeable barrier model (RPBM) fitting with a dictionary approach (`CreateRPBMDictionary`,
  `FitRPBMDictionary`, `GetRPBMValues`).
- Correction for gradient non-linearity (`GradientCoilTensor`, `MakeGradientMaps`, `MakeGradientDerivatives`).
- Export of b-values and b-vectors (`ExportBvalvec`).

### Tractography

- Import and export of MRtrix `.tck` tract files (`ImportTCK`, `ExportTCK`).
- Tract measures: length, curvature and curvature maps (`TractLength`, `TractCurvature`, `TractCurvatureMap`).
- Endpoint density maps, filtering tracts by length, and fitting tract segments.
- Choice of integration method for each tract step (`TractMethod`).

### Fasciculation analysis

- Activation size can be set in voxels or as an area in mm² (`SizeMethod`), and with a voxel size the analysis also
  reports ROI volume, chance per dm³ and activation areas.

### Denoising and reconstruction

- `P2SDenoise` (self-supervised neural network denoising, formerly `NNDeNoise`), harmonic denoising of tensor data (`HarmonicDenoiseTensor`), and
  denoising of spectroscopic imaging data (`DenoiseCSI`).
- Noise prewhitening, SENSE smoothing and weighting, k-space averaging and zero-padding options in reconstruction.

### Spectroscopy and simulation

- `SpectraSimulator`, phase shifting of FIDs (`PhaseShiftFid`), and sinc profile plots.
- Extended phase graph (EPG) simulation now accounts for refocusing pulse phase (`EPGRefocussingPhase`).

### Plotting and general

- `LoessPlot` for smoothed trend plots with prediction intervals, `MeanStdRange`, `MinMaxRange`, and legend images.
- Proper support for the notebook dark mode.
- Log-file handling (`SetLogFile`, `SaveLogFile`), demo notebook helpers (`OpenDemonstrationNotebook`,
  `SetDemoDirectory`), and temporary file management (`ClearQMRIToolsTemp`).
- Performance updates.

### Breaking changes

- **Requires Wolfram Language 15.0 or newer.**
- **Renamed and removed functions.** Several function and option names were corrected to consistent spelling and a few functions were replaced. Code that uses the old names must be updated.
- **Muscle-BIDS config files.** Segmentation label handling and config merging changed. Check existing study configs
  against the examples in `QMRITools/BIDS Example/`.

------------------------------------------------------------------------

## QMRITools 4.0

This version introduced deep learning muscle segmentation.

- **SegmentationTools** (new): a CNN segmentation framework built on UNET-style networks, with
  `TrainSegmentationNetwork`, patch-based training and inference, data augmentation, `SegmentData` for automatic
  muscle segmentation, and scripts for running segmentation in batch.
- **ScientificColorData** (new): perceptually uniform scientific colour maps (Crameri) for all plotting functions.
- **Legacy** (new): old and rarely used functions moved here to keep the main toolboxes lean.
- Dixon reconstruction was rewritten to be more flexible: three-point phase correction, DCT phase unwrapping, a more
  stable bipolar reconstruction, and multi-stack Dixon.
- T2 mapping added to Muscle-BIDS processing, and scripts for multi-stack Dixon and DTI.
- Faster tractography and new tractography features.
- An interactive gradient design tool that runs in the free Wolfram Player.
- Requires Wolfram Language 14.0 or newer.

------------------------------------------------------------------------

## QMRITools 3.0

This version moved the toolbox to the modern paclet structure and introduced Muscle-BIDS.

- **MuscleBidsTools** (new): the Muscle-BIDS data structure and pipeline to convert, process and merge muscle MRI
  data from DICOM to analysis-ready data.
- **LoggingTools** (new): logging of processing steps for reproducible batch pipelines.
- **FasciculationTools** (new): detection and analysis of muscle fasciculations.
- **TractographyTools** (added in 2.x): fibre tractography directly in the Wolfram Language.
- **SpectroTools** and **ReconstructionTools** (added in 2.x): MR spectroscopy fitting and plotting, basis spectra,
  CSI tools, Hankel SVD, and reconstruction of raw list data with coil combination and SNR maps.
- **TaggingTools** (added in 2.x): analysis of tagged MRI.
- Neural-network-based denoising, multi-contrast registration with Elastix, and better Dixon handling of bipolar
  acquisitions.
- Refactored to the new paclet layout (`Kernel/*.wl`) with rebuilt documentation and guide pages.
- Requires Wolfram Language 13.0 or newer.

------------------------------------------------------------------------

## QMRITools 2.0

With this version DTITools was renamed to **QMRITools**, to reflect a scope wider than diffusion.

- Renamed from DTITools to QMRITools, with complete documentation and a full documentation rebuild.
- **TensorTools**, **JcouplingTools** and **CoilTools** added.
- Demonstration notebook and example data.
- Slice visualisation tools, rescaling of segmentations, and a complex EPG signal model that accounts for off-resonance.
- Published in the Journal of Open Source Software
  ([Froeling, JOSS 2019](https://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)).

------------------------------------------------------------------------

## DTITools 1.x

The original toolbox, focused on diffusion tensor imaging of muscle, with toolboxes for cardiac, denoising, Dixon,
Elastix registration, gradients, IVIM, masking, NIfTI, plotting, processing, relaxometry and simulation.
