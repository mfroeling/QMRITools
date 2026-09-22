# MuscleBidsTools: agent guide to the Muscle-BIDS pipeline

This is a map of [QMRITools/Kernel/MuscleBidsTools.wl](../QMRITools/Kernel/MuscleBidsTools.wl) (~3700 lines), the
config-driven processing pipeline that goes from DICOM to per-muscle Excel sheets. Read it before you change
anything in that file.

- Code state: as of 2026-09-22, including the `Dixon-A` work.
- Line numbers (`L123`) drift over time. Grep for the function name if a number no longer matches.
- Index of all agent docs: [AGENTS.md](../QMRITools/AGENTS.md). Code style for this repo: [CodeStyle.md](CodeStyle.md).
  The segmentation step (§6.5) calls `SegmentData`: [SegmentationTools.md](SegmentationTools.md).
- `BIDS-config.docx` in the repo root is the user manual (v1.1, 2025-07). It is **outdated**. Section 9 lists how it
  differs from the code. When they disagree, the code wins.

---

## 1. How the user runs it

The user runs everything from a notebook in the study folder (`dir`, the folder that holds `config.json`):

```wolfram
<< QMRIToolsDev`                 (* loads the dev paclet from D:\werk\workspace\QMRITools\QMRITools *)
SetDirectory[NotebookDirectory[]];

ViewConfig[dir]                  (* tabbed view of config.json + defaults *)
SelectSubjects[dir]              (* checkbox UI over dicom folders -> copy list to clipboard *)
subs = All;

BidsDcmToNii[dir, ProcessSubjects -> subs]                          (* 1 dicom -> raw nii (dcm2niix) *)
ViewProtocolNames[dir, ProcessSubjects -> subs]                     (* check Labels vs ProtocolName *)
MuscleBidsConvert[dir, ProcessSubjects -> subs]                     (* 2 raw nii -> BIDS named nii *)
MuscleBidsProcess[dir, ProcessSubjects -> subs, VersionCheck -> False]      (* 3 fitting *)
MuscleBidsMerge[dir, ProcessSubjects -> subs, VersionCheck -> False]        (* 4 join stacks + register *)
MuscleBidsSegment[dir, ProcessSubjects -> subs, VersionCheck -> False]      (* 5 CNN segmentation *)
MuscleBidsTractography[dir, ProcessSubjects -> subs, VersionCheck -> False] (* 6 tracts *)
MuscleBidsAnalysis[dir, ProcessSubjects -> subs, BidsOutputImages -> "None"](* 7 xlsx + images *)
```

Every step is **idempotent**: it skips work that a `*_check.json` file marks as done (section 7). The user reruns
the whole chain repeatedly as new subjects arrive.

Example configs (real studies, outside the repo):

| Config | What it exercises |
| --- | --- |
| `D:\Werk\Research\7_motion_long\config.json` | Classic: 6 `Stacks` per contrast, `Dixon` (raw complex echoes), DTI with tractography plus harmonic denoise, EPGT2. **No duplicate Type+Suffix**, so there are no keys in the Targets. Flat `analysis`. |
| `D:\Werk\Research\MOTOR\config.json` | Many body regions (Head/Calf/Body/Shoulder). `Dixon-A` everywhere (both recon and raw source). `tse` type. `Chunks` and `Volumes` classes. Duplicates everywhere, so **every Target starts with a dataset key**. No `analysis` section. |
| `D:\Werk\Research\TWITCH\config.json` | `Dixon-P` with `Types` vector (4D recon nii). Segment `Dimensions: "2D"` and `VoxSize`. Duplicates, so keys. |
| `D:\Werk\Research\ext - Bochum\BIDS\config.json` | `Dixon-B` with `EchoTime`, mese with `EchoTime` (4D nii), `conversion.Version: "17"`, **nested `analysis`** (one block per key `OS`/`US`), and keys that contain `_` (see quirk Q1). |

---

## 2. Folder layout and data flow

```text
<dir>/                                   study root; functions SetDirectory here, all folders relative
  config.json
  01_sourcedata/  (folders.dicomData)    <subj> or <subj>_<ses> dicom folders (+ optional per-subject config.json)
  02_rawdata/     (folders.rawData)
    DcmToNii_<DateName>.log
    sub-<s>/ses-<e>/
      raw/                               dcm2niix output, files deleted after Convert (DeleteAfterConversion)
      dix/ dwi/ quant/ anat/ miss/       MuscleBidsConvert output (folder = BidsType[Type])
      sub-<s>_ses-<e>_BIDSConvert.log
      sub-<s>_ses-<e>_config.json        copied per-subject config (if present)
  03_derivatives/ (folders.derivedData)  sub/ses/<typefolder>/   MuscleBidsProcess output (native space)
  04_merged/      (folders.mergeData)    sub/ses/<typefolder>/ + seg/   Merge, Segment, Tractography
  05_analysis/    (folders.analysis)     sub/ses/*.xlsx|wxf|jpg + All_<DateName>.xlsx|wxf
```

| # | Public function | Reads (`folders.*`) | Writes | Worker | `Method->` string | Log file |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | `BidsDcmToNii` | dicomData | rawData/sub/ses/raw | `BidsDcmToNiiI` | `"BidsDcmToNii"` | `rawData/DcmToNii_<date>.log` (one per run) |
| 2 | `MuscleBidsConvert` | rawData (`InFolder` = `raw`) | rawData (in place) | `MuscleBidsConvertI` | `"MuscleBidsConvert"` | `rawData/sub/ses/<nam>_BIDSConvert.log` |
| 3 | `MuscleBidsProcess` | rawData | derivedData | `MuscleBidsProcessI` | `"MuscleBidsProcess"` | `derivedData/sub/ses/<nam>_BIDSProcess.log` |
| 4 | `MuscleBidsMerge` | derivedData | mergeData | `MuscleBidsMergeI` | `"MuscleBidsMerge"` | `mergeData/.../<nam>_BIDSMerge.log` |
| 5 | `MuscleBidsSegment` | mergeData | mergeData | `MuscleBidsSegmentI` | `"MuscleBidsSegment"` | `..._BIDSSegment.log` |
| 6 | `MuscleBidsTractography` | mergeData | mergeData | `MuscleBidsTractographyI` | `"MuscleBidsTractography"` | `..._BIDSTractography.log` |
| 7 | `MuscleBidsAnalysis` | mergeData | analysis | `MuscleBidsAnalysisI` | `"MuscleBidsAnalysis"` | `analysis/sub/ses/<nam>_BIDSAnalysis.log` |

`<nam>` = `GenerateBidsName[<|sub, ses|>]`, for example `sub-001_ses-1`.

Stage dependencies are carried only by files on disk. There is no in-memory state between stages. Stage N finds its
input by **regenerating the file name** that stage N-1 wrote (section 4). Changing a naming rule in one stage therefore
breaks every later stage.

---

## 3. Architecture: the call chain

Every public step uses the same three-tier overload chain, then goes through one shared loop:

```text
MuscleBidsX[dir]                          -> MuscleBidsX[dir, GetConfig[dir]]
MuscleBidsX[dir, config_Association]      -> SetDirectory[dir];
                                             MuscleBidsX[ConfigLookup[config,"folders",IN],
                                                         ConfigLookup[config,"folders",OUT],
                                                         config["datasets"]   (* "analysis" for step 7, "conversion" for step 1 *)
                                                         , opts];
                                             SetDirectory[old]
MuscleBidsX[inFol, outFol, datDis]        -> BidsFolderLoop[inFol, outFol, datDis, Method -> "MuscleBidsX", opts]
                                             (MuscleBidsAnalysis also joins all per-subject .wxf into All_<date>.xlsx afterwards)

BidsFolderLoop (L986)                     one loop, all 7 steps
  for fol in subject/session folders:            (* SelectBids[inFol,"ses"]; dicom dirs for BidsDcmToNii *)
      ass  = SubNameToBids[fol, met]              (* <|sub, ses|> *)
      out  = GenerateBidsFolderName[outFol, ass]  (* outFol/sub-x/ses-y *)
      {custom, datDis} = CheckConfig[fol, out]    (* per-subject config patch, copied to out *)
      datDis = CheckDataDescription[MergeConfig[datDisIn, datDis], met]  (* validate + add Key/HasDuplicate/InFolder/OutFolder/Class/Suffix *)
      SetLogFile/ImportLog/ShowLog
      Switch[met,
        "BidsDcmToNii"       -> BidsDcmToNiiI[fol, out, datDisIn]
        "MuscleBidsAnalysis" -> MuscleBidsAnalysisI[{fol, outFol}, #, versCheck, imOut] & /@ datDis
        _ -> for type in datDis:                 (* each dataset description *)
               if type has $Failed key -> log + skip
               for folIn in SelectBids[fol, type["InFolder"]]:   (* e.g. .../sub-1/ses-1/dix *)
                   Echo[...]; Switch[met,
                     Convert  -> MuscleBidsConvertI[folIn, type, delete]
                     Process  -> MuscleBidsProcessI[{folIn, outFol}, type, versCheck]
                     Merge    -> MuscleBidsMergeI[{folIn, outFol}, {type, datDis}, versCheck]
                     Segment  -> MuscleBidsSegmentI[{folIn, outFol}, {type, datDis}, versCheck]
                     Tracto   -> MuscleBidsTractographyI[{folIn, outFol}, {type, datDis}, versCheck, tractMet]]
      SetLogFile[]
```

Important consequences:

- **Workers get the full `datDis` list (`allType`)** in Merge, Segment and Tractography, so they can resolve
  Targets that point at *other* datasets.
- The loop order is subject → dataset (config order) → input folder. In Merge the target dataset is simply
  rebuilt from derivatives when its merged file does not exist yet, so dataset order is not critical.
- `outFol` passed to workers is the **stage root** (for example `04_merged`), not the subject folder. Workers rebuild
  sub/ses paths themselves through `PartitionBidsFolderName[folIn]` → `{root, parts}`.
- `BidsFolderLoop` options: `Method`, `ProcessSubjects`, `VersionCheck`, `DeleteAfterConversion`,
  `BidsTractographyMethod`, `BidsOutputImages`. **`BidsIncludeSession` is not among them** (see Q6).

### Subject selection (`ProcessSubjects`)

`subs` is matched by `SubNameToBids` (L304). Both sides are turned into `<|"sub"->.., "ses"->.., "suf"->{}|>` and
compared:

| Input | Result |
| --- | --- |
| `"001_2"` (dicom folder style) | `<|sub->001, ses->2|>` |
| `"001"` (no session) | `<|sub->001, ses->001|>`, **forced session `001`** |
| `"sub-001_ses-2"` | `<|sub->001, ses->2|>` |
| folder `02_rawdata\sub-001\ses-2` | `<|sub->001, ses->2|>` |

So `"001"` will **not** match `ses-2`. Pass what `SelectSubjects` copies (the dicom folder names). Subject and session
names must not contain `-`, `_`, `.` or spaces.

---

## 4. The naming system (the core of the whole file)

### Entities and types

```wolfram
bidsName  = {"sub","ses","vol","stk","chunk","rep","acq","part","type","suf"};   (* L230 *)
bidsClass = {"Volume","Volumes","Stacks","Repetitions","Chunks","Acquisitions","Mixed"};
bidsTypes = <|"T1w"|"T1w-FS"|"T2w"|"T2w-FS" -> "anat", "megre"|"tse" -> "dix",
              "mese"|"T1"|"T2"|"wT2" -> "quant", "dwi" -> "dwi", "seg" -> "seg"|>;  (* unknown -> "miss" *)
```

- `GenerateBidsName` (L380) writes entities in the fixed order `sub, ses, vol, stk, rep, chunk, acq, part`, then `type`,
  then the `suf` list, all joined with `_`.
- `GenerateBidsFileName[root, parts]` = `root/sub-x/ses-y/<BidsType[type]>/<GenerateBidsName>`, with **no extension**.
  Callers append `.nii`, `.json`, `_<output>.nii`, or use `ConvertExtension`.
- `PartitionBidsName` (L340) is the inverse. `k-v` tokens become entities. The first remaining token becomes `type`
  **only if** it is a key of `bidsTypes`. The rest go to `suf`. A token is dropped from `suf` if it contains any
  entity key as a substring (so a suffix named `subX` would vanish).
- `PartitionBidsFolderName[fol]` = `{text before "sub-", PartitionBidsName[path parts containing "-"]}`. For a bare
  folder it returns the key `"parts"->{}` instead of `"suf"`. That is harmless.

Examples, checked with wolframscript on 2026-09-22:

```text
PartitionBidsName["sub-001_ses-1_stk-DIXOS_megre_dix_wat"]
  -> <|sub->001, ses->1, stk->DIXOS, type->megre, suf->{dix, wat}|>
GenerateBidsFileName["04_merged", <|sub->001, ses->1, chunk->BodyDIX, type->seg, suf->{auto,megre,dix,outph}|>]
  -> 04_merged\sub-001\ses-1\seg\sub-001_ses-1_chunk-BodyDIX_seg_auto_megre_dix_outph
GetClassName["Chunks", "DIX-SRC-Body-1"] -> chunk -> DIXSRCBody1
```

### Which entity carries what, per stage

| Stage | Entity used | Value |
| --- | --- | --- |
| Convert → Process outputs (per acquisition) | `GetClassName[Class, label]`: Volume/Volumes→`vol`, Stacks→`stk`, Chunks→`chunk`, Repetitions→`rep`, Acquisitions→`acq` | `StringStrip[label]` (`-_. ` removed) |
| Merge / Segment / Tractography / Analysis (merged) | `stk` if Class=="Stacks", else `chunk` | `StringStrip[datasetKey]`, **only if `HasDuplicate`** |

- Raw file: `02_rawdata/sub-1/ses-1/dix/sub-1_ses-1_stk-DIXON1_megre_dix_real.nii`
- Derived: `03_derivatives/.../dix/sub-1_ses-1_stk-DIXON1_megre_dix_wat.nii` (plus `..._megre_dix_check.json`)
- Merged: `04_merged/.../dix/sub-1_ses-1[_stk-DIXOS]_megre_dix_wat.nii`

### `Key` and `HasDuplicate` (added by `CheckDataDescription`, L874)

- `Key` = `StringStrip[datasetName]`, so `"DIX_OS"` becomes `"DIXOS"`.
- `HasDuplicate` is **one global flag**. It is True for **every** dataset as soon as **any** two datasets share
  `{Type, Suffix}`. Then every merged file name carries the key, and every Target/Segmentation list in the config must
  start with a dataset key (MOTOR, TWITCH, Bochum). Without duplicates there are no keys and Targets have none
  (motion_long).

### `BuildBidsNameFromConfig` (L607–668): turning config lists into paths

It has three forms:

1. `BuildBidsNameFromConfig[{root, parts}, datType]`: the merged output base name of this dataset
   (`..._[key]_<Type>_<Suffix>`).
2. `BuildBidsNameFromConfig[{root, parts}, datType, con_String]`: the list of **per-stack input** `.nii` names for
   contrast `con` (`<Type>_<Suffix>_<con>`, one per Label; for `Volumes` only the first Label).
3. `BuildBidsNameFromConfig[{outRoot, (inRoot,) parts}, {datType, all}, targetList]`: resolves a config target list
   such as `["DIX_OS","megre","dix","wat"]` or `["megre","dix","outph"]`:
   - If `First[target]` is a known `Key`, then key = it, type = next element, suf = rest. Otherwise key = this dataset's
     own key and the whole list is type+suf.
   - Without `inRoot` it returns one merged name `.nii`. With `inRoot` it returns the per-stack derived names of the
     target dataset, used for `tarStack` in Merge.
   - Non-keyed targets pick the target dataset via `#["InFolder"] === First[suf]`. That only works when the dataset's
     `Suffix` equals its BIDS folder (in practice `dix`). See Q2.

---

## 5. Config reference (what the code actually reads)

`ConfigLookup[assoc, section, key]` = `assoc[section][key]`, falling back to `defaultConfig[section, key]` (L480).
Keys not listed below are **ignored silently**. `MergeConfig` is a recursive deep merge: the patch wins, nested
associations merge.

### Top level

| Key | Used by | Notes |
| --- | --- | --- |
| `folders.{dicomData,rawData,derivedData,mergeData,analysis}` | all | defaults `01_sourcedata` … `05_analysis` |
| `conversion.Version` | BidsDcmToNii → `DcmToNii[UseVersion->]` | `1` (default) → asset `DcmToNii`; anything else → asset `"DcmToNii-"<>ToString[v]` (for example `"Own"`, `"17"`) |
| `datasets` | steps 2–6 | one entry per dataset; the entry name becomes `Key` |
| `analysis` | step 7 | **flat** (contains `Analysis`, so it is wrapped as `{"Default"->...}` with no key) or **nested** (one block per name, where the name becomes the key used in file names) |

### Dataset description (validated by `CheckDataDescription[list, met]`, L904)

| Key | Required for | Default / derived |
| --- | --- | --- |
| `Label` | Convert (mandatory) | string or list of ProtocolNames (without `WIP `). Also used to glob files in Process/Merge (`*<StringStrip[label]>*.json`) |
| `Type` | all | must be a `bidsTypes` key, else message and folder `miss` |
| `Class` | – | default `"Volume"`. Volume needs a string Label; list classes need a list |
| `Suffix` | – | default `""` |
| `InFolder` | – | Convert: `raw`; others: `BidsType[Type]` |
| `OutFolder` | – | always `BidsType[Type]` (set, but hardly used) |
| `Merging` | Merge (mandatory) | a dataset without it is logged as "Wrong data description" and skipped in Merge |
| `Process`, `Segment`, `Tractography` | optional | a missing section means the step logs a skip for that dataset |

### `Process` (MuscleBidsProcess, and a few keys for Convert)

| Key | Default | Read in | Meaning |
| --- | --- | --- | --- |
| `Method` | – | Convert (megre/tse), Process | Dixon: `Dixon`, `Dixon-B`, `Dixon-P`, `Dixon-S`, `Dixon-A`. dwi: `DTI`. mese: `EPGT2`, `Exp` |
| `Types` | `{wat,fat,inph,outph}` (Process only) | Convert Dixon-P/S, Process branch 1 | vector = order of volumes in a 4D nii (Dixon-P); matrix `[[suffix, fileTag]..]` = separate files (Dixon-P: tag is the last `_` token of the file name; Dixon-S: `SeriesDescription == label_tag`) |
| `EchoTime` | – | Convert | Dixon-B: `[TE1, dTE]` ms. mese: echo spacing ms (switches mese conversion to "single 4D file" mode) |
| `Masking` | 5 | Process dwi/mese | `Mask` threshold |
| `FlipPermute` | `{{1,1,1},{x,y,z}}` | DTI prep | `FlipGradientOrientation` |
| `RegistrationDimension` | `"3D"` | DTI prep | `"2D"` → `RegisterCardiacData` (and disables split) |
| `RegistrationMethod` | `"bspline"` | DTI prep | `MethodReg` |
| `SplitRegistration` | False | DTI prep | `RegisterDiffusionDataSplit` (left/right separate) |
| `GradientCorrection` | False | DTI | coil name string (for example `"WA1"`), needs `Offset` in json (written by Convert) |
| `IVIMCorrection` | True | DTI | IVIM fit then `IVIMCorrectData` before the tensor fit |
| `FasciculationDetection` | False | DTI | `FindActivations` on the normalised `reg` data |
| `Settings` | – | EPGT2 | `[[exName, exAngle, {..}], [refName, refAngle, {..}]]` → `GetPulseProfile`, or plain numbers `[ex, ref]` |
| `Shift` | – | EPGT2 | `true` → fat profile shift computed from the pulse settings |

### `Merging`

| Key | Default | Meaning |
| --- | --- | --- |
| `Target` | – | target list (section 4). Its space becomes the common space |
| `Moving` | – | one suffix of *this* dataset that is registered to the target |
| `Process` | – | suffixes carried to merged space. Multi-dim ones (`tens`, `data`, `filt`, `fasc`, `reg`) are handled specially |
| `Overlap` | 0 | slices; int or `[overTarget, overNative]`. `0` → no motion correction, no padding |
| `Padding` | 0 | pad the overlap |
| `Motion` | True | `JoinSets[MotionCorrectSets->]`, and whether to register when target is same dataset |
| `Reverse` | False | stack order reversed (Body Chunks in MOTOR/TWITCH) |
| `SplitRegistration` | **True** | `RegisterDataSplit` for stacks ≠ first. **Configs write `"Split"`, which is never read** (Q3) |

### `Segment`

| Key | Default | Meaning |
| --- | --- | --- |
| `Target` | – | target list, or list of target lists → one segmentation each |
| `Location` | `"Body"` | `SegmentData` location (`Legs`, `LegsHip`, `LowerLegs`, `Shoulder`, …) |
| `Device` | `"CPU"` | `TargetDevice` |
| `Dimensions` | `"3D"` | `"2D"` or `"3D"` → `SegmentationDimension` |
| `Method` | `Automatic` (symbol) | omit it for CNN. `"Registration"` registers an existing segmentation (`Target`, `Moving`, `Segmentation`, `VoxSize`). See Q9 |
| `VoxSize` | Automatic | **only used by the Registration method** (TWITCH sets it for Automatic, where it has no effect) |

### `Tractography`

| Key | Default | Meaning |
| --- | --- | --- |
| `Target` | – | tensor, for example `["dwi","dti","tens"]`. The output base replaces the last element with `trk` |
| `Stopping` | – | `[[targetList, [min,max]], ...]`. Rescaled to the tensor voxel size if needed |
| `Segmentation` | – | seg target list (for example `["seg","auto","megre","dix","outph"]`). Needed for segmenting tracts and for HarmonicDenoise |
| `FlipPermute` | `{{1,1,1},{x,y,z}}` | `TensorFlips` / `TensorPermutations` |
| `TractLength` `{15,500}`, `TractAngle` 25, `TractSeed` 0.66 (<1 → `Scaled`), `TractStep` `"Automatic"` (non-number → auto), `SegmentLength` `{15,500}`, `BoneLabel` 100 (muscles = 1..BoneLabel, bones = BoneLabel+1..+30), `HarmonicDenoise` False |

### `analysis` block (per analysis entry)

| Key | Meaning |
| --- | --- |
| `Segmentation.Type` | seg name list (with key when nested) |
| `Segmentation.Labels` | `[nMax, what]`. `what` = `"Legs"`, `"Arm"`, or a path to a custom label file |
| `Analysis.Types` | `[[ (key,) type, suf.., [params..]] ..]`. The last element is threaded (`Thread`) |
| `Analysis.TractBased` | subset analysed only inside the tract density mask (`dwi_dti_trk_dens`, unitized) |
| `Analysis.MaskErosion` | default True. **Must sit directly in `Analysis`; the documented `Analysis.Options.MaskErosion` is ignored** (Q4) |
| `Analysis.UseFilter` | `[range, typeList]` → extra `_F` columns with the seg masked by `Mask[filterMap, range]` |
| `Analysis.Class` | `"Stacks"` → key entity `stk`, otherwise `chunk` (Q5) |
| `Images.Reference` | map used for slice positions and seg background |
| `Images.QuantImages` | `[[typeList, colorFunc?, range?, [min,max,label]?] ..]` |
| `Images.TractImages` | `..._trk_plot` (a `.wxf` scene written by Tractography) |

### Per-subject config overrides

Put `config.json` (or `<sub-x_ses-y>_config.json`) in a subject folder. Its content is a **patch on the `datasets`
section** (not a full config). `CheckConfig` finds it in the stage input folder, merges it, and **copies it to the
stage output folder**, so it propagates dicom → raw → derivatives → merged. The log then shows
`***** Using custom config *****`. `CheckConfigLabels[dir, subs]` lists Labels missing per subject, including the
custom patches. Analysis does not use custom configs.

---

## 6. The stages in detail

### 6.1 BidsDcmToNii (L1104, L1126)

- Loops over **all directories** in dicomData, not `sub-` folders. `ResetLog[]` runs once per run.
- Output `rawData/sub-x/ses-y/raw`. **Skips if that folder is not empty.** To redo a subject, empty `raw`.
- `DcmToNii[{in, out}, MonitorCalc->False, UseVersion->conversion.Version]`. On Windows the naming is
  `%s_%t_%p` (series_time_protocol). `GetProtocolNames` hides files whose name contains the `_hhmmss.fff` time
  pattern.

### 6.2 MuscleBidsConvert → `MuscleBidsConvertI` (L1178)

For each Label (`nameIn`), it globs `*<label with spaces→_>*.json` in the `raw` folder, imports all json files, and
selects files with `GetJSONPosition[json, {{key, value}..}, sortKey]` (case-insensitive, `WIP` stripped). If several
files match, `CheckPos` takes the **last** one and logs a warning. Output goes to the BIDS type folder with
`MergeJSON[{info, infoExtra}]`, where `infoExtra` holds the conversion software/version plus
`Volume/Stack/Repetition/OverLap...` and `ForthDimension`/`DataClass`. Used raw files are removed with
`DelFiles[..., del]` (nii, nii.gz, json, bval, bvec).

| Type | Method / condition | What it does | Output suffixes |
| --- | --- | --- | --- |
| megre, tse | `Dixon-P` + `Types` vector | one 4D nii, volume i → `Types[[i]]` | `Types` |
| | `Dixon-P` + `Types` matrix | separate files matched by last `_` token | `Types[[All,1]]` |
| | `Dixon-S` | Siemens, `SeriesDescription == label_tag` | `Types[[All,1]]` |
| | `Dixon-B` | Philips single 4D (mag/real/imag/phase interleaved). Needs `EchoTime`. Special case when `ImageType` ends in `FFE` (B0-only) | `""`, `ph`, `real`, `imag` |
| | `Dixon-A` | **auto**: if any json `ImageType` ends in WATER/FAT/IN_PHASE/OUT_OF_PHASE → recon path; otherwise raw MIXED/PHASE/REAL/IMAGINARY per echo, sorted by `EchoNumber`, scaled `1000/2047`, phase `π(x-2047)/2047` | `wat fat inph outph`, or `"" ph real imag` |
| | `Dixon` | per-echo files per ImageType (Mixed/Phase/Real/Imaginary). Philips R12+: `ORIGINAL` (GE-like reorder) vs `DERIVED` (MOTOR/TWITCH reorder with RotateRight). Older: Transpose | `""`, `ph`, `real`, `imag` |
| | other Method | **no default branch, silently nothing** | – |
| dwi | Class `Volumes` | several labels joined (`ConcatenateDiffusionData`) into one file named after the first label | Suffix |
| | other | Siemens → match `SeriesDescription`, else `ProtocolName` | Suffix (+ `.bval/.bvec`, `Offset` in json) |
| mese | `Process.EchoTime` present | single 4D nii (`ImportNiiT2` unless json has EchoTime+EchoTrainLength or Siemens). TE = `EchoTime·(1..n)` | Suffix |
| | otherwise | single-echo files sorted by `EchoNumber`, cut to `AcquisitionNumber`/`EchoTrainLength`. R12 reorder | Suffix |
| anat / other | – | "Unknown datatype for conversion" | – |

On non-Windows systems files are written uncompressed and `CompressNiiFiles` runs at the end (`compress` flag, L208).
This applies to every stage.

### 6.3 MuscleBidsProcess → `MuscleBidsProcessI` (L1817)

Input files are `*<StringStrip[label]>*.json` in the `InFolder`. **Sets** are the distinct `PartitionBidsName` results.
For megre/tse they are reduced to the first suffix, so all `_dix_*` files of one stack form one set.
`outFile = GenerateBidsFileName[derivedRoot, set]`, and outputs are written as `outFile_<name>.nii` + `.json`.

**Export-by-name idiom:** `ExportNii[ToExpression[con<>#], vox, outFile<>"_"<>#<>".nii"] & /@ outTypes`, where
`con = Context[con]` is the Block's private context. The **local variable name is the output suffix**. To add an
output: declare the variable in the `Block` list, give it exactly the suffix name, and append that name to `outTypes`.

| Type | Method | Check file | Outputs |
| --- | --- | --- | --- |
| megre/tse | branch 1: `Dixon-S`, `Dixon-P`, `Dixon-A` (recon files present) | `outFile` "done" | the found `Types`, plus `r2star` (if t2star), `watfr`/`fatfr` (from fatfr file or `DixonToPercent`) |
| | branch 2: `Dixon`, `Dixon-B`, `Dixon-A` (real/imag or mag/ph files present) | same | `real imag mag ph b0i t2stari b0 t2star r2star inph outph wat fat watfr fatfr itt res snr sig`, plus for Dixon/Dixon-A: `phbp phi phbpt phii phbpi` (and `dbond` if >6 echoes, CallDB fat model) |
| dwi | *prep* (always, even without Method) | `outFile_prep` "done" | `raw den reg filt sig snr0 snr` (`reg`/`filt` with bval/bvec) |
| | `DTI` (reruns when prep ran this call) | `outFile` "done" | `data tens res out s0 l1 l2 l3 md fa rd` + IVIM `mean adci fri s0i` + `field` (grad corr) + `fasc fascm norm` |
| mese | `EPGT2` | `outFile` "done" | `data t2w t2f b1 wat fat fatfr res` |
| | `Exp` | same | `data t2 s0` |

The Dixon branches are **gated `If`s, not a `Switch`**. `Dixon-A` runs both. Branch 2 overwrites branch 1's values
when raw data exists, and each branch's "files not found" guard makes it a no-op when its data is absent. This design
is intentional (see CodeStyle.md, "avoid duplication"). Export and the check file run once after both branches.

DTI prep pipeline: `SortDiffusionData` → `FlipGradientOrientation` → `Mask` → `PCADeNoise` → `SNRCalc` →
register (2D `RegisterCardiacData` / 3D `RegisterDiffusionData(Split)`) → `P2SDenoise` (`filt`).
DTI pipeline: `filt` → mask → optional gradient-nonlinearity correction (`MakeGradientDerivatives`) → optional IVIM
correction → `TensorCalc[iWLLS, RobustFit]` → `ParameterCalc` → optional fasciculation analysis.

### 6.4 MuscleBidsMerge → `MuscleBidsMergeI` (L2523)

1. Names: `outFile` (merged base), `tarFile` (merged target), `tarStack` (per-stack derived target files), `movStack`,
   `processStacks` (per Process type, per stack). A missing file → log + `Return[]`. Then a "done" `CheckFile` → skip.
2. Settings: `Overlap` gives `{overT, overM}` (target-space and native overlap). `overT==0` → `pad=0`, `motion=False`.
   Registration method by `InFolder`: `dix` → rigid, `quant` → rigid+affine, else (dwi) → rigid+affine+bspline.
3. Target: if a merged `tarFile` exists and the target is a different dataset → import and `SplitSets` back into
   stacks. Otherwise import `tarStack` and `JoinSets` (motion correct, normalise), then `SplitSets`.
4. Moving: import all process types per stack. Multi-dim types (`tens data filt fasc reg`) are flattened to 3D volumes
   (`lengthMD` remembers the 4th dimension).
5. Registration, per stack, **unless `!motion && sameType`**. `sameType` means tarStack == movStack, so the target is
   this dataset's own contrast. Register the `Moving` contrast to the target stack. Transform the scalar types into
   target space and the `native` types (`tens`, `fasc`) into the target rescaled to native voxel size. The first stack
   (last one when `Reverse`) uses non-split functions and is not registered when sameType.
6. `JoinSets` per type (`overM` for native types). `nonQuant` types (`inph outph wat fat s0 mean`) are `Ramp`ed and
   normalised. The multi-dim types are unflattened again.
7. Export `outFile_<type>.nii` + json (source json + `Merge *` settings) and write the check file.

### 6.5 MuscleBidsSegment → `MuscleBidsSegmentI` (L2865)

- Check file: `..._[key]_seg_auto_<Type>`, **one per dataset, not per target**.
- `Method` Automatic: for each target, `SegmentData[{data, vox}, Location, SegmentationDimension, TargetDevice]` →
  `seg/…_[key]_seg_auto_<target...>.nii`. A missing input sets the check status to `"error"`, so the dataset reruns
  next time.
- `Method "Registration"`: registers an existing `Segmentation` using `Target`/`Moving` into target space →
  `seg_reg_…`. This branch looks experimental (Q9).
- Manual segmentations follow the pattern `seg_man_…` in the `seg` folder. Analysis and Tractography use whatever the
  config names.

### 6.6 MuscleBidsTractography → `MuscleBidsTractographyI` (L3014)

The check file holds a **state machine** on `outFile = ..._<target minus last>_trk`:

```text
(none) --FiberTractography--> "track" --SegmentTracts+maps--> "seg" --> "done"
```

- `BidsTractographyMethod` controls the steps: `"Full"` runs both, `"Tractography"` only tracking, `"Segmentation"`
  only segmenting (needs an existing `.trk`).
- Tracking writes `_trk.trk`. With `HarmonicDenoise` (needs the segmentation) it also writes
  `_trk_con/amp/tensH.nii.gz` and `_trk_har.trk`, and the segmentation step then uses the `_har` tracts.
- Segmenting writes `_trk_seg.trk`, the maps `_trk_{dens,leng,ang,seed,curv}.nii.gz`, and the plot scene
  `_trk_plot.wxf`, which Analysis `TractImages` imports.

### 6.7 MuscleBidsAnalysis → `MuscleBidsAnalysisI` (L3294)

This step runs per analysis entry: `"Default"` (flat) or each nested key.

- **xlsx part** (check `…_xls`): segmentation → `SelectSegmentations[1..n]` → volume and cross-section → label names
  from `GetAssetLocation["MuscleLegLabels"]` and similar → per `Analysis.Types` file: `GetMaskData` with the median+IQR
  (mean+SD for `trk_seed`/`trk_dens`). Scaling: `fatfr`×100 (dix, t2), `dix t2star`×1000. TractBased types use
  `seg ∩ Unitize[dens]`. `MaskErosion` → `DilateMask[seg,-1]`. Output `analysis/sub/ses/<sub_ses[_chunk-KEY]>.xlsx`
  + `.wxf`.
- **image part** (check `…_img`): `BidsOutputImages` = `All | Quantitative | Segmentation | Tractography | None`.
  Slice positions come from `Images.Reference`. Output: 2D quant jpgs with optional `LegendImage`, seg 2D / 3D
  (`_vol`) / `_grid` jpgs, tract 3D jpg.
- After all subjects: `All_<DateName>.xlsx/.wxf` joins every per-subject `.wxf` whose name does not start with `All`.

---

## 7. Idempotency: check files and VersionCheck

- `MakeCheckFile[base, rules]` writes `base_check.json` with `Check`, the rules, `ProcessingSoftware`, `Version`
  (`$InstalledVersion`) and `Date`. It lives in LoggingTools.wl.
- `CheckFile[base, status, verCheck]` is True when the file exists, `Check === status`, and (if `verCheck`) the
  stored version is not older than the installed one.
- **To force a redo, delete the relevant `*_check.json`** (or run with `VersionCheck->True` after a version bump). The
  exceptions are DcmToNii (empty the `raw` folder) and Convert (it has no check file: it reruns while raw files exist,
  and `DeleteAfterConversion` removes them).

| Stage | Check file base | Status |
| --- | --- | --- |
| Process | `derived/.../<set>` and `<set>_prep` (dwi) | `done` |
| Merge | merged `outFile` | `done` |
| Segment | `seg/..._seg_auto_<Type>` | `done` / `error` |
| Tractography | `..._trk` | `track` → `seg` → `done` |
| Analysis | `analysis/.../<name>_xls`, `_img` | `done` |

---

## 8. Debugging

- `QMRITools`MuscleBidsTools`$debugBids = True` makes every `debugBids[...]` print (paths, dimensions, positions).
- The log is written continuously: `AddToLog` exports the whole `$Log` on every call. `ShowLog[]` opens a live
  window. Levels 0–5 are indentation, and `True` adds a timestamp.
- `Echo[DateString[], folIn <> " - " <> Key]` prints once per input folder per dataset.
- Helpers: `ViewConfig[dir]`, `ViewProtocolNames[dir]` (duplicate protocol names shown in red), `CheckConfigLabels[dir, All]`,
  `GetProtocolNames[dir, sub]`.
- Most "nothing happened" cases are a **name mismatch**: a Target list that does not resolve to an existing file.
  Turn on `$debugBids` and compare the printed expected path with the files on disk.

---

## 9. Differences from BIDS-config.docx (v1.1)

| Manual says | Code does |
| --- | --- |
| Classes are `Chunks` / `Acquisitions` | `Volume` (default), `Volumes` (dwi concat), `Stacks`, `Chunks`, `Repetitions`, `Acquisitions`, `Mixed` |
| After merging the key becomes `acq-<KEY>` | the key is `stk-<KEY>` (Stacks) or `chunk-<KEY>` (everything else), and only when `HasDuplicate` |
| `megre` → dix | `tse` → dix as well |
| Dixon method `Dixon` only | `Dixon`, `Dixon-B`, `Dixon-P`, `Dixon-S`, `Dixon-A`. mese: `EPGT2`, `Exp` |
| `Analysis.Options.{MaskErosion,TractWeighting}` | `Options` is not read. `MaskErosion` must sit directly under `Analysis`. Tract weighting is hard-coded False |
| Merging key `SplitRegistration` (defaults) | correct. The example configs use `Split`, which is ignored |
| `ViewConfigFile` | the function is `ViewConfig` |
| not documented | Segment `Dimensions`/`Method`/`VoxSize`; Tractography `TractLength`/`TractAngle`/`TractSeed`/`TractStep`/`SegmentLength`/`BoneLabel`/`HarmonicDenoise`; Analysis `UseFilter`; nested analysis blocks; per-subject configs; `conversion.Version` strings |

---

## 10. Known quirks and likely bugs

These were verified by reading the code. The ones marked ✔ were also confirmed with wolframscript. Ask the user before
"fixing" any of them, because some may be relied upon.

| # | Where | Issue |
| --- | --- | --- |
| Q1 ✔ | `BuildBidsNameFromConfig` L641 vs `CheckDataDescription` L895 | `Key` is `StringStrip`ped (`DIX_OS`→`DIXOS`), but `isKey` compares the **unstripped** `First[target]`. Dataset names containing `_ - . space` are therefore never recognised as keys in Merge/Segment/Tractography targets (the Bochum config). Analysis strips correctly. Both lines come from the same refactor commit (509747d6). |
| Q2 | `BuildBidsNameFromConfig` L661 | Non-keyed targets resolve the target dataset through `InFolder === First[suf]`, which only works when Suffix == BIDS folder (`dix`). |
| Q3 | Merge L2606 | Reads `Merging.SplitRegistration` (default True). All example configs write `"Split"`, so split registration is always on. |
| Q4 | Analysis L3308 | `MaskErosion` is read from `Analysis`, not `Analysis.Options`, so the configs' `Options` block does nothing. `tractWeighting` is hard-coded False (L3309), and the config key is spelled `TractWeigthing`. |
| Q5 | Analysis L3319 | Key entity is `stk` only if `Analysis.Class == "Stacks"`. Merge uses the **dataset** Class. With Stacks + duplicates (Bochum) and no `Analysis.Class`, the analysis looks for `chunk-KEY` while the files are named `stk-KEY`. |
| Q6 | `BidsDcmToNii` options | `BidsIncludeSession` is declared but never forwarded (`BidsFolderLoop` does not know it, and `SubNameToBids` always uses the default True). |
| Q7 | Merge L2566 | The "moving files" existence check tests `tarStack` again instead of `movStack`. |
| Q8 ✔ | `CheckDataDescription` L940–944 | Class `Chunks` is missing from the label/class `Switch`, so the result is unevaluated, `If[!cls,…]` does not fire, and it passes by accident. `ListQ[ass["Label"] && Length[ass]>1]` has a misplaced bracket and works only because `And[list, True]` reduces to `list`. |
| Q9 | Segment "Registration" L2960 | `Table[..., {i,1,3}]` over `tar[[{i}]]` looks like debug/experimental code. The same function logs `datType["key"]` (lowercase) at L2891, which shows `Missing`. If `Segment.Method` is given as the string `"Automatic"`, no branch matches but the check file is still written as `done`. |
| Q10 | Analysis L3373 | `segF` is not in the `Block` locals, so it leaks as a global. |
| Q11 | Convert dwi | A missing bval/bvec only logs a warning, but the export then runs with an undefined `data`/`vox`. |
| Q12 | Convert megre/tse | An unknown `Process.Method` has no default `Switch` branch, so conversion silently does nothing. Anat types are not converted at all ("Unknown datatype"). |
| Q13 | Process/Convert globbing | `*<StringStrip[label]>*` also matches longer labels (`DIXON1` matches `DIXON10`). |
| Q14 | `BidsFolderLoop` for BidsDcmToNii | `CheckDataDescription` runs on the `conversion` block, which is not a dataset description. The result is unused. |

---

## 11. Recipes

**Add a processing method** (for example a new mese fit):

1. Add a `Switch` case under the type in `MuscleBidsProcessI`. Copy the guard structure (CheckFile → json exists →
   nii exists → import → compute → export → MakeCheckFile → compress).
2. Name each result variable exactly as its output suffix, add it to the Block locals, and set `outTypes`.
3. If conversion differs, add the case in `MuscleBidsConvertI` under the same `Type`.
4. If the new outputs must be merged, the user lists them in `Merging.Process`. If they are multi-dim or must stay
   native, extend `multDim` / `native` / `nonQuant` in `MuscleBidsMergeI`.
5. If the method shares behaviour with an existing one, prefer gated `If`s or a broader alternation over copying a
   body (the Dixon-A precedent).

**Add a config key:** give it a default in `defaultConfig` (L480) and read it with `ConfigLookup[datType, section, key]`.
Update this file's section 5.

**Add a new `Type`:** add it to `bidsTypes`, add a `mandatory` list for it in `CheckDataDescription` if needed, then add
cases in Convert and Process.

**Change naming:** stop and trace every stage that regenerates the name (Convert output → Process sets → Merge
`BuildBidsNameFromConfig` → Segment/Tractography targets → Analysis `fileName`). Existing study folders on disk depend
on the current names.

---

## 12. Function index

| Function | Line | Role |
| --- | --- | --- |
| `debugBids`, `$debugBids`, `dataToLog`, `compress` | L195–208 | debug/log helpers, compression flag |
| `WipStrip`, `StringStrip` | L240–243 | remove `WIP`, remove `-_. ` |
| `BidsType`, `BidsValue`, `BidsString`, `GetClassName` | L250–289 | entity helpers |
| `SubNameToBids` | L296 | sub/ses from a name or folder (subject matching) |
| `PartitionBidsName`, `PartitionBidsFolderName` | L340, L366 | name → parts |
| `GenerateBidsName/FolderName/FileName` | L380–419 | parts → name/path |
| `SelectBids`, `SelectBidsFolders/Subjects/Sessions` | L430–469 | find `sub-*/ses-*/<folder>` dirs |
| `defaultConfig`, `ConfigLookup`, `CheckConfig`, `GetConfig`, `MergeConfig` | L480–600 | config |
| `BuildBidsNameFromConfig` | L607–668 | config target lists → file names |
| `ImportJSON`, `GetJSONPosition`, `MergeJSON` | L679–709 | json selection/merge |
| `ViewConfig`, `MakeTable`, `ProcessImage` | L720–757 | config viewer UI |
| `ViewProtocolNames`, `GetProtocolNames`, `SelectSubjects`, `CheckConfigLabels` | L764–865 | inspection UI |
| `CheckDataDescription` | L874, L904 | validate and complete dataset descriptions |
| `BidsFolderLoop` | L986 | the shared driver |
| `BidsDcmToNii` / `I` | L1100 / L1126 | step 1 |
| `MuscleBidsConvert` / `I`, `CheckPos`, `DelFiles` | L1152 / L1178, L1768, L1778 | step 2 |
| `MuscleBidsProcess` / `I` | L1791 / L1817 | step 3 |
| `MuscleBidsMerge` / `I` | L2497 / L2523 | step 4 |
| `MuscleBidsSegment` / `I` | L2839 / L2865 | step 5 |
| `MuscleBidsTractography` / `I` | L2988 / L3014 | step 6 |
| `MuscleBidsAnalysis` / `I` | L3256 / L3294 | step 7 |

External dependencies: `AddToLog`, `SetLogFile`, `ImportLog`, `ShowLog`, `ResetLog`, `MakeCheckFile`, `CheckFile`
(LoggingTools.wl); `DcmToNii`, `ImportNii*`, `ExportNii`, `CompressNiiFiles` (NiftiTools.wl); `ConvertExtension`,
`NiiFileExistQ`, `DateName` (GeneralTools.wl); `JoinSets`/`SplitSets` (ProcessingTools.wl); `SegmentData`
(SegmentationTools.wl); `FiberTractography` and the tract maps (TractographyTools.wl); `DixonReconstruct`,
`EPGT2Fit`, `TensorCalc`, `IVIMCalc` and the registration functions from their respective files.
