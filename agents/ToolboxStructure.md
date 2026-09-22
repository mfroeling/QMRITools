# QMRITools toolbox structure and loader: agent guide

This guide shows how the paclet is organised, how it loads, how the ~28 subpackages depend on each other, and how the
documentation and build work. Read it before touching `Kernel/QMRITools.wl`, `PacletInfo.wl`, the package list, or
the docs, and whenever something "works on the second load only".

- Index of all agent docs: [AGENTS.md](../QMRITools/AGENTS.md). Style: [CodeStyle.md](CodeStyle.md) (§2 has the file
  anatomy of a single subpackage). Code maps: [MuscleBidsTools.md](MuscleBidsTools.md),
  [SegmentationTools.md](SegmentationTools.md).
- Numbers below were **measured with wolframscript on 2026-09-22** (version 4.11.0) unless marked otherwise.

---

## 1. Repository and paclet layout

```text
QMRITools/                         repo root
  QMRITools/                       the paclet (PacletDirectoryLoad points here)
    PacletInfo.wl                  name, version (4.11.0), WolframVersion "15.0+", extensions + assets
    Kernel/                        QMRITools.wl (loader) + 29 subpackage .wl files
    Documentation/English/         Guides/ (28 .nb), ReferencePages/Symbols/ (1030 .nb), no Tutorials
    Resources/                     demo notebooks, All-Functions.nb/.pdf, DemoData.zip, SCMv8txt.zip (colour maps), GradientGUI
    NeuralNetworks/                segmentation/classifier nets + label .txt files   (*.wlnet git-ignored)
    Applications/<SystemID>/       elastix, transformix, dcm2niix, pigz binaries      (git-ignored)
    BIDS Example/                  example config.json + notebook for the Muscle-BIDS pipeline
    AGENTS.md, CLAUDE.md           agent entry point
  agents/                          agent guides (this folder)
  docs/                            qmritools.com Jekyll site (index.md + images), not the paclet docs
  scripts/                         Install_QMRITools.wls, Segment_Nii.wls (CLI segmentation), segmentGUI.nb
  buildscripts.nb                  build/release notebook (§6)
  build/, build_tools/             build output and release notes (git-ignored)
  BIDS-config.docx                 Muscle-BIDS user manual (outdated, see MuscleBidsTools.md §9)
```

A fresh git clone does **not** contain the network weights or the external binaries (both git-ignored). They come
with the released `.paclet`, or have to be copied in.

### PacletInfo.wl assets (resolved by `GetAssetLocation[name]`, GeneralTools)

| Group | Asset names → files |
| --- | --- |
| Resources | `Logo`, `Functions` (All-Functions.nb), `Demo`, `Demo_Unet`, `DemoData`, `ColorData` (SCMv8txt.zip), `GradientTool` |
| NeuralNetworks | `Body` (position/side classifier), `Seg{Shoulder,Hip,Thigh,Leg}Muscle{2D,3D}`, `SegArmMuscle3D`; train label files `*TrainLabels` (N1–N7); output label files `MuscleLabels`, `MuscleLegLabels`, `MuscleShoulderLabels`, `MuscleArmLabels`, `MusclesAllLabels`, `MusclesLegAllLabels` (see SegmentationTools.md §2) |
| Windows-x86-64 | `Elastix`, `Transformix`, `ElastixLib`, `pigz`, `DcmToNii` (current dcm2niix), `DcmToNii-Own`, `DcmToNii-25/24/23/21/20/19/17` (older dcm2niix versions, picked with `DcmToNii[UseVersion->]` / BIDS `conversion.Version`) |
| MacOSX-x86-64 and MacOSX-ARM64 | `Elastix`, `Transformix`, `ElastixLib` (×2), `DcmToNii`. Both point at the same universal binaries in `Applications/MacOSX-x86-64` |
| Linux | **none declared.** Registration and DICOM conversion have no bundled binaries on Linux |

---

## 2. The loader, `Kernel/QMRITools.wl`, step by step

`Get["QMRITools`"]` (or `<< QMRIToolsDev``, see §4) runs this file:

1. **`$HistoryLength = 0`** for the whole session (saves memory with large arrays). This is a user-visible side
   effect: `Out[n]` and `%` history stop working after loading.
2. `BeginPackage["QMRITools`"]`. If `$VersionNumber < 13.3`, a warning dialog appears. Note that PacletInfo requires
   15.0+.
3. Flags (set them **before** `Get` to take effect):
   - `$Verbose`: echo every loading step and the full function list
   - `$Legacy`: also load `Legacy.wl`, 68 public symbols
   - `$Loaded`, `$LoadedColor`: internal state that marks a reload
4. `$SubPackages` is the ordered package list. The grouping comments in it ("core", "specific data types", "lots of
   dependencies") are **not** a real dependency order (§3). `$Contexts` = `"QMRITools`" <> #`.
   `$InstalledVersion` comes from `PacletFind`.
5. **Reload branch** (only if `$Loaded` is True, i.e. a second `Get` in the same kernel):
   - `Quiet[Get /@ $Contexts]` to get all names;
   - for each context: `Unprotect` + `ClearAll` all public symbols, and `Remove` any `Global`` symbol with the same
     name (created when a function name was typed before loading).
6. **Load loop**: `Get[#] & /@ $Contexts`. After `ScientificColorData`, the colour maps are unzipped from the
   `ColorData` asset and registered once per kernel (`$LoadedColor`).
7. `$Loaded = True`. `$ContextsFunctions = {context, names}`. **Every public symbol gets `Protected` +
   `ReadProtected`.**
8. Then some globals are unprotected again so users can set them: `$lastElastixTemp`, `$debugElastix`, `$debugUnet`,
   `$debugBids`, `$debugDenoise`, `$AmaresB1Regularization`, `$AmaresT1Values`, `$Log`, `$LogFile`. **A new
   user-settable `$` variable must be added to this list**, or assigning it fails with `Protected` errors.

Every subpackage starts with

```wolfram
BeginPackage["QMRITools`X`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`X`"}]]];
```

so **every package sees every other package's public symbols**. Public symbols are exactly those that get a
`::usage` in the front section. That rule held: 0 public symbols lack a usage string.

Measured: a fresh load takes about 2 s; 953 public symbols (1021 with Legacy); 29 `QMRITools` contexts on
`$ContextPath`; no public short name is defined in two packages.

---

## 3. How loading really works: the depth-first cascade and the two passes

This is non-obvious and load-bearing.

- `BeginPackage[ctx, needs]` calls `Needs` on every listed context **immediately**. So the first `Get` in the loop
  (GeneralTools) starts a **depth-first cascade**. GeneralTools begins and needs ScientificColorData, which begins and
  needs LoggingTools, and so on down `$SubPackages`. The **body** of the last package (AmaresTools) is read first,
  then the bodies unwind back up to GeneralTools. Measured `$Packages` order confirms this.
- When a deep package's body is read, the packages above it in the chain have **not yet run their usage section**, so
  their public symbols don't exist yet. A reference such as `AddToLog` inside MuscleBidsTools is then created as a
  **shadow symbol** `QMRITools`MuscleBidsTools`Private`AddToLog`.
  - Measured: **223 shadow symbols** in 17 private contexts (MuscleBidsTools 79, SegmentationTools 44, SpectroTools 19,
    ShapeTools 18, …).
- The load loop then **`Get`s every package a second time**. Now all public symbols exist, every definition is
  re-read with the correct symbols, and the pass-1 definitions are overwritten (same left-hand sides).
  - Measured: **0 live definitions reference a shadow symbol.** The shadows are inert, empty leftovers.

Consequences and rules:

- **The second pass is required, not redundant.** The loader comment "on a first load … there is no need to load
  everything twice" refers only to the reload branch. In practice every package except GeneralTools is read twice on
  every fresh load. Never "optimise" the loop to `Needs`, and don't reorder it assuming it is a dependency order.
- **Top-level code in a package (outside definitions) runs in both passes.** In pass 1, cross-package symbols are
  still unresolved. Keep top-level code to System symbols and the package's own symbols (for example
  `compress = ($OperatingSystem === "Windows")`). Put anything that needs another package inside a function, or in
  the loader after the loop (like the colour-data step).
- Seeing a symbol like `QMRITools`X`Private`SomePublicName` while debugging is normal. It is only a problem if a live
  definition uses it. The check script in §8 finds those.
- **Dev reload caveat.** The reload branch clears only **public** symbols. Private definitions survive a second `Get`:
  - A renamed or re-patterned private helper keeps its old `DownValues` next to the new ones.
  - Memoised private caches persist (for example `GetNeuralNetI`, the segmentation network cache).
  - When private helpers are removed or their argument patterns change, **restart the kernel** (or `ClearAll` that
    private context) instead of relying on the reload.

---

## 4. Development loading

- `D:\Werk\workspace\QMRIToolsDev` is a tiny paclet whose kernel file does
  `PacletDirectoryLoad["D:\\werk\\workspace\\QMRITools\\QMRITools"]; PacletDataRebuild[]; Get["QMRITools`"]`.
  `<< QMRIToolsDev`` therefore loads the working tree instead of the installed paclet.
- In scripts, use the same two lines without the rebuild. Test scripts go outside the repo.
- Code that starts sub-kernels must load the same version there. `LaunchTrainingKernels` (SegmentationTools) checks
  whether the paclet location contains `"workspace"` and loads `QMRIToolsDev`` or `QMRITools`` accordingly.
- Because of the `Protected` attribute, redefining a public function in a notebook fails unless you `Unprotect` it
  first. Reloading via `<< QMRIToolsDev`` is the normal way to pick up edits.

---

## 5. Package map and measured dependencies

Columns: non-blank lines; public **f**unctions / **o**ptions (split by whether the usage says "is an option", which
is approximate); the guide page's own description; which other packages it **uses** (references their public
symbols); and **used by** (how many packages reference it).

| Package | Lines | f / o | Description (main guide) | Uses | Used by |
| --- | --- | --- | --- | --- | --- |
| GeneralTools | 1508 | 72/14 | shared helpers (assets, crop, cut, rescale, monitor, conversions) | none | **25** |
| ScientificColorData | 118 | 2 | registers Crameri scientific colour maps at load | none | 0 (via `ColorData` names) |
| LoggingTools | 239 | 11/1 | `$Log`, `AddToLog`, check files | General | 1 |
| MaskingTools | 520 | 23/8 | masks and segmentations | Elastix, General | **18** |
| NiftiTools | 1098 | 22/13 | nii import/export, `DcmToNii` | General, Masking, Plotting, Processing | 4 |
| ElastixTools | 1363 | 12/29 | registration via elastix | General, Masking, Nifti, Tensor | 6 |
| PlottingTools | 2016 | 23/22 | visualisation | General, Masking, Tensor | **10** |
| MuscleBidsTools | 2913 | 26/7 | Muscle-BIDS pipeline | 17 packages | 0 (top-level application) |
| NeuralNetworkTools | 1063 | 22/12 | UNet building, losses, class encoding | General, Masking, Plotting | 1 |
| DixonTools | 1105 | 13/24 | Dixon reconstruction | General, Gradient, Masking, Processing, Relaxometry | 2 |
| IVIMTools | 268 | 5/8 | IVIM fitting | General, Gradient | 1 |
| DenoiseTools | 1098 | 12/19 | PCA / other denoising | Elastix, General, Masking, Processing, Tensor, Tractography | 2 |
| CardiacTools | 1885 | 25/26 | cardiac MRI | General, Masking, Plotting, Processing | 1 |
| RelaxometryTools | 1033 | 12/17 | T2 / T1 / T1rho, EPG | General, Masking, Processing | 3 |
| GradientTools | 1676 | 29/17 | diffusion gradients, b-matrices, nonlinearity | Elastix, General, Plotting, Tensor | 6 |
| TensorTools | 1196 | 30/13 | DTI fitting and parameters | General, Gradient, Masking, Processing | 8 |
| JcouplingTools | 824 | 15/10 | J-coupled spectra simulation | General | 2 |
| SpectroTools | 1245 | 30/12 | MRS fitting and processing | Cardiac, General, Jcoupling, Plotting, Reconstruction | 2 |
| ReconstructionTools | 810 | 30/14 | basic recon, coil combination | Denoise, General | 3 |
| TractographyTools | 1126 | 25/16 | fibre tractography | General, Masking, Plotting, Tensor | 3 |
| ProcessingTools | 1261 | 33/37 | general processing (`JoinSets`, `SplitSets`, …) | Elastix, General, Masking, Reconstruction, Tensor, Tractography | **10** |
| FasciculationTools | 388 | 5/7 | fasciculation detection | General, Gradient, Masking, Plotting | 1 |
| SimulationTools | 1033 | 21/11 | DWI/DTI/Dixon/EPG simulation | General, Gradient, Masking, Processing, Relaxometry, Spectro, Tensor | 1 |
| CoilTools | 251 | 7/3 | coil/SNR analysis | General, Masking, Nifti, Plotting, Processing | 0 |
| TaggingTools | 428 | 2/2 | tagging MRI | General, Masking | 0 |
| SegmentationTools | 2466 | 33/30 | CNN segmentation | General, Masking, NeuralNetwork, Nifti, Processing | 1 (MuscleBids) |
| ShapeTools | 630 | 16/12 | shape models | Elastix, General, Masking, Plotting | 0 |
| AmaresTools | 546 | 9/4 | AMARES MRS fitting | Dixon, Jcoupling, Reconstruction, Spectro | 0 |
| Legacy (optional) | 3616 | 68 total | old functions, incl. a commented-out old `MakeUnet` | – | – |

Reading the graph:

- **Foundation** (widely used): GeneralTools → MaskingTools → Plotting / Processing / Tensor / Elastix / Gradient.
- **Leaf applications** (nothing depends on them): MuscleBidsTools, CoilTools, TaggingTools, ShapeTools, AmaresTools.
  MuscleBidsTools sits on top of almost everything.
- **There are dependency cycles** (from the "Uses" column):
  - direct: MaskingTools ↔ ElastixTools, TensorTools ↔ ProcessingTools, GradientTools ↔ TensorTools;
  - longer: NiftiTools → ProcessingTools → ElastixTools → NiftiTools, and
    TractographyTools → TensorTools → ProcessingTools → TractographyTools.

  This is why the all-see-all
  `BeginPackage` plus two passes (§3) is needed. **No package order could load in a single pass.** Don't try to
  "fix" the order.
- Where a helper lives is not always where you'd guess. For example `JoinSets` and `SplitSets` are in
  ProcessingTools, `DcmToNii` is in NiftiTools, `GetAssetLocation`, `ConvertExtension` and `MonitorFunction` are in
  GeneralTools, `CheckFile` and `MakeCheckFile` are in LoggingTools. Grep for `^Name\[` in `Kernel/` before assuming.

---

## 6. Documentation and build

**Paclet documentation** (`QMRITools/Documentation/English`):

- **Main guide** `Guides/QMRITools.nb`: an abstract, then one `GuideText` cell per package (a link plus a one-line
  description, the source of the descriptions in §5).
- **Package guides** `Guides/<Package>.nb` (28, hand-maintained): a title, an abstract, and flat
  `InlineGuideFunctionListing` groups of function links. Only GeneralTools uses `GuideFunctionsSubsection`
  subsections. **Options are not listed on guides**, only functions.
- **Reference pages** `ReferencePages/Symbols/<Symbol>.nb` (1030): one per public symbol and option, generated from
  the usage strings. The notebooks carry `DocuToolsSettings` metadata, then get edited and built.
- There are no tutorials. The `docs/` folder at the repo root is the website, not paclet docs.

**Build** (`buildscripts.nb`, run by hand in the front end; uses `PacletTools`` and `DocumentationBuild``):

- `ApplicationBuild[location, opts]`:
  - copies the `.wl` files (plus the extension files with `RebuildExtentions`) and `PacletInfo.wl` to
    `build/QMRITools`;
  - `RebuildDoc` → `PacletDocumentationBuild`;
  - `MakePaclet` → `CreatePacletArchive`;
  - `MakePacletSite`;
  - `PushRelese` → GitHub release through the `gh` CLI, with `ReleaseTitle` and `ReleaseNote`;
  - `InstallApplication`.
- Helpers: `BuildHTMLDoc` (writes `build/QMRITools-html`), `CheckDocumentation`, `FuncTable`, `ShowStatus`.
- Release notes live in `build/note_<version>.md`.

**Documentation gaps measured 2026-09-22** (fixing them is doc work for the user, not code):

- The main guide links to `LecacyTools` (a misspelling, and there is no such guide). It **doesn't link** the
  `AmaresTools` or `NeuralNetworkTools` guides.
- Stale reference pages for symbols that no longer exist: `AugmentMask` (renamed `MaskedPretraining`),
  `MaxPatchSize`, `SegmentationResolution`, `ZeropadData`, `tempDir`, and a page named
  ``QMRITools`NeuralNetworkTools`$debugUnet``.
- Public symbols without a reference page: `MakeTrainData` (new), plus Legacy `PlotRespiract`, `ReadBrukerDiff`,
  `ROIMask`, `ShiftPar`, `SpectraFitResult`.
- Functions missing from their package guide:
  - AmaresTools guide lists **none** of its 9 functions.
  - SegmentationTools: `CopyTrainedNetwork`, `FreezeEncoderLayers`, `MakeChannelGrid`, `MakeTrainData`,
    `RunMuscleMap`, `SplitDataForSegmentation`.
  - NeuralNetworkTools: `ActivationLayer`, `AnalyzeNetworkFeatures`, `CELossFunction`, `CELossLayer`,
    `ClassConfidence`.
  - PlottingTools: `GenerateRotationFrames`, `LegendImage`, `LoessPlot`, `ShowLink`.
  - GeneralTools: `FindMiddle`, `FitGradientMap`, `LightDarkV`, `MonitorFunction`, `ParseCommandLine`.
  - 1–4 each in Nifti, Elastix, Denoise (`PCADeNoise`), Gradient, Spectro, Reconstruction, Tractography, Simulation,
    Masking, Logging, Cardiac, Relaxometry, Fasciculation, Coil, Tagging, Shape, MuscleBids.

  The script's option detection is approximate: a few "missing" names such as `DeNoiseKernel`, `MuscleMapPath`,
  `MeshesPerRow`, `UseMask` and `SphereColor` are really options.

---

## 7. Checklists

**Add a public function to an existing package:**

1. Add the `::usage` string in the front "Usage Notes" section, and a `::tag` message in "Error Messages" if needed.
2. Add `Options[f]` and `SyntaxInformation[f]`, then the definition in its own `(* ::Subsubsection::Closed:: *)` cell.
3. Keep top-level code free of cross-package symbols (§3).
4. The docs are the user's build step: a reference page and a link on the package guide page. Mention it; don't
   hand-edit notebooks unasked.
5. Reload with `<< QMRIToolsDev``. After changing private helper patterns, restart the kernel.

**Add a user-settable global** (`$debugX`, a tuning constant): define it in the package, then **add it to the
Unprotect list at the end of `QMRITools.wl`**.

**Add a new subpackage:**

1. Create `Kernel/NewTools.wl` with the standard `BeginPackage[... Complement[QMRITools`$Contexts, ...]]` header.
2. Add it to `$SubPackages` in `QMRITools.wl`, after the packages it mostly uses.
3. Add `Guides/NewTools.nb` and a `GuideText` link on the main guide.
4. Update §5 of this guide and AGENTS.md if it gets its own agent guide.

**Add an asset:** put the file under `Resources/`, `NeuralNetworks/` or `Applications/<SystemID>/`, register it in
`PacletInfo.wl` (per SystemID for binaries), and fetch it with `GetAssetLocation["Name"]`. Remember the Linux gap.

---

## 8. Reproducing the measurements

All checks are short wolframscript scripts that load the dev paclet
(`PacletDirectoryLoad[...]; Get["QMRITools`"]`) and then:

- **Load order**: `Select[$Packages, StringStartsQ[#, "QMRITools"] &]`. The newest is first, so the reversed list is
  the begin order.
- **Shadow symbols**: private names whose short name equals any public short name,
  `Select[Names[ctx <> "Private`*"], MemberQ[publicShortNames, Last@StringSplit[#, "`"]] &]`.
- **Live shadow references**:
  - clear `Protected` and `ReadProtected` on all names;
  - for every name, collect `DownValues`/`UpValues`/`SubValues`/`OwnValues`/`Options` with
    `Function[s, Hold @@ {{DownValues[s], …}}, HoldAll]`, keeping `Hold @@` so the values actually evaluate;
  - check `FreeQ` against each shadow symbol.
- **Dependencies**: `Cases[defs, sym_Symbol :> Context[Unevaluated[sym]], Infinity, Heads -> True]` over each
  package's public and private names, restricted to `$Contexts`.
- **Reference pages vs symbols**: compare `FileBaseName /@ FileNames["*.nb", ".../ReferencePages/Symbols"]` with the
  public short names (load with `QMRITools`$Legacy = True` first to include Legacy).
- **Guide coverage**: `StringCases[Import[guide, "Text"], "paclet:QMRITools/ref/" ~~ x : WordCharacter .. :> x]`
  compared with each package's public non-option names. `NotebookImport` of cell text needs a front end and fails
  in wolframscript.
