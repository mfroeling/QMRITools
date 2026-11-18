# Muscle-BIDS vs BIDS and ORMIR-MIDS: deviations and possible fixes

Where the files that [MuscleBidsTools.wl](../QMRITools/Kernel/MuscleBidsTools.wl) writes differ from
[BIDS](https://bids-specification.readthedocs.io/) and from ORMIR-MIDS, and what could be changed later. Nothing here is
decided: every fix is the author's call. The code map is [MuscleBidsTools.md](MuscleBidsTools.md), the naming rules are
in its section 4.

- **ORMIR-MIDS**: Santini et al., *JBMR Plus* 2026;10(3):ziag013, doi 10.1093/jbmrpl/ziag013 (Froeling co-author).
  Compared against its Table 1 and text on 2026-11-18. It extends BIDS and MIDS for musculoskeletal imaging.
- **MIDS** (Saborit-Torres et al., arXiv:2010.00434) adds a modality level above the BIDS datatype folders
  (`mr/anat`). ORMIR-MIDS replaces that with `mr-anat` and `mr-quant` at the datatype level, so MIDS is covered by the
  ORMIR-MIDS column.
- **BIDS** points are from the BIDS entity rules as known on 2026-11-18, not yet checked line by line against the
  current specification. Check them before a publication cites them.
- Scope: the converted data in `02_rawdata` and, where noted, the derivatives. ORMIR-MIDS does not cover the later
  pipeline folders (derived maps per stack, merged data, analysis).

## 1. Where Muscle-BIDS already agrees

| Aspect | Muscle-BIDS | BIDS | ORMIR-MIDS |
| --- | --- | --- | --- |
| Formats | NIfTI + one JSON sidecar per image | yes | yes |
| `sub-<label>/ses-<label>/<datatype>/` tree | always, default session `001` | yes | yes (session level left to the user) |
| Complex parts in separate files | magnitude, phase, real, imaginary each one file | yes | yes |
| `EchoTime` as an array in the sidecar | yes | – | yes |
| `4D (x, y, z, echo)` for megre/mese | yes | no (see D3) | yes |
| `FourthDimension` key | yes since 2026-11-18 (was `ForthDimension`; Process copies either spelling) | – | yes |
| Diffusion with `.bval`/`.bvec` in `dwi/` | yes | yes | no (see D8) |
| `seg` folder for segmentations | yes | – | yes |

## 2. Deviations

| # | What | Muscle-BIDS now | BIDS | ORMIR-MIDS |
| --- | --- | --- | --- | --- |
| D1 | Complex-part naming | trailing words `_ph`, `_real`, `_imag`; magnitude has none | entity `part-mag\|phase\|real\|imag` before the suffix (`_phase` suffix deprecated) | `part-` entity, as BIDS |
| D2 | Suffix | several words after the entities: type + dataset `Suffix` + map (`megre_dix_real`), type in lower case | exactly one suffix at the end (`MEGRE`, `MESE`, `dwi`, `T2map`) | one suffix (`megre`, `mese`, `t2`) |
| D3 | Echoes | 4D file, echoes in the 4th dimension | one file per echo, entity `echo-<index>` | 4D, as Muscle-BIDS |
| D4 | Class entities | `stk-`, `vol-`, `rep-`, `chunk-`, `acq-` with the label text (`chunk-DIXSRCBody1`) | `chunk-<index>`, `run-<index>`, `acq-<label>`; `stk`, `vol`, `rep` do not exist | not defined |
| D5 | Datatype folders | `dix` (megre, tse), `quant` (mese, T1, T2, wT2), `dwi`, `anat`, `seg`, `miss` | `anat`, `dwi`, `fmap`, ... (`MEGRE`, `MESE` in `anat`) | `mr-anat` (megre, mese, DESS), `mr-quant` (t1, t2, wt2, planned diff and vel), `seg` |
| D6 | `EchoTime` unit | seconds (dcm2niix) | seconds | **ms** |
| D7 | Sequence fields in the sidecar | only what dcm2niix writes; not added by the conversion | sequence-specific | megre: `WaterFatShift` (pixels), `MagneticFieldStrength`; mese: `RefocusingFlipAngle` |
| D8 | Diffusion layout | BIDS `dwi/` with `.bval`/`.bvec` | yes | planned: `mr-quant`, suffix `diff`, directions in the JSON (`FourthDimension: EncodingDirection`) |
| D9 | Extra metadata files | none | – | optional full DICOM header JSON and patient JSON |
| D10 | Segmentation labels | 3D label map, names from a separate label file, not in the JSON | – | 3D or 4D, JSON `Labels` (array of names) |
| D11 | Dataset files | no `dataset_description.json`, no `participants.tsv` | both at the root (required / recommended) | inherited from BIDS |
| D12 | Derivatives | own folders `03_derivatives`, `04_merged`, `05_analysis`; map names as trailing words (`_fatfr`, `_t2star`) | `derivatives/<pipeline>/` with its own `dataset_description.json` (`GeneratedBy`), entities `desc-`, `space-`, qMRI suffixes (`T2starmap`, `MWFmap`, ...) | maps in `mr-quant` |
| D13 | Folder names of the study root | `01_sourcedata` ... `05_analysis` (configurable) | `sourcedata/`, raw at the root, `derivatives/` | as BIDS |

## 3. Possible fixes

Effort: **S** additive, no existing file changes; **M** a code change in a few places; **L** renames files or json
keys that later steps and existing studies read, so it needs a migration path (see section 4).

| Fix | For | Effort | What to change | Risk |
| --- | --- | --- | --- | --- |
| F1 Add `dataset_description.json` (raw and each derived folder, with `GeneratedBy` = QMRITools version and step) | D11, D12 | S | write once per run in `BidsFolderLoop` or the public step | none |
| F2 Add `participants.tsv` from the subject/session folders | D11 | S | one export after `BidsDcmToNii` or `MuscleBidsConvert` | none |
| F3 Add `WaterFatShift`, `MagneticFieldStrength`, `RefocusingFlipAngle` to the converted sidecar when the DICOM json has them (and log when not) | D7 | S | `infoExtra` in `MuscleBidsConvertI` | none |
| F4 Add the label names (`Labels` array) to the segmentation json | D10 | S | `MuscleBidsSegmentI` json, names from `MuscleLabelToName` | none |
| F5 Add `EchoTimeUnits: "s"`, or a second key `EchoTimeMs`, instead of changing the unit | D6 | S | conversion `infoExtra` | none; switching the unit itself would be L (Process reads `EchoTime` in seconds, see DataConventions.md) |
| F6 `.bidsignore` for the Muscle-BIDS-only files (logs, check files, per-subject configs) | D11 | S | written with F1 | none |
| F7 Use `part-mag\|phase\|real\|imag` for complex parts | D1 | L | `GenerateBidsName` entity order, `PartitionBidsName`, the export suffixes in Convert, the imports in Process, `BidsCheckFiles` | file names of every Dixon study |
| F8 Map the classes to BIDS entities with an index: `stk` → `chunk-<n>`, `rep` → `run-<n>`, `vol` → none; keep the label in the json | D4 | L | `GetClassName`, `bidsClass`, naming in Merge; index from the order of `Label` | all file names and check files; Analysis `Class` key |
| F9 One BIDS suffix: type as suffix (`MEGRE`, `MESE`), dataset `Suffix` to `desc-` (derivatives) or `acq-` (raw) | D2, D12 | L | `GenerateBidsName`, `PartitionBidsName`, every `Target` list in configs | all configs and file names |
| F10 Datatype folders as ORMIR-MIDS (`mr-anat`, `mr-quant`) | D5 | M–L | `bidsTypes` only, but every study changes folder | folder names of every study; option instead of default? |
| F11 Echo layout | D3 | – | keep 4D (ORMIR-MIDS); an export to per-echo BIDS files could be an extra tool | none if kept |
| F12 Diffusion in ORMIR-MIDS form | D8 | – | wait until ORMIR-MIDS finalises it; BIDS `dwi` is the safer choice now | – |
| F13 Derivatives as `derivatives/qmritools/` with `desc-`/`space-` and qMRI suffixes | D12, D13 | L | folder defaults and map names; the `folders` config already allows renaming | analysis column names follow map names |
| F14 Optional DICOM header JSON (ORMIR "extra" file) | D9 | M | dcm2niix cannot write it; needs a DICOM header dump at `BidsDcmToNii` | storage size, patient data in files |

Suggested order: F1–F6 first (additive, no study breaks, and they close the most visible gaps for a reviewer). Then
decide on the naming group F7–F9 together, because each renames the same files and one migration is cheaper than three.

## 4. If names change: migration

Any L fix needs, in the same release:

- the old and new names readable for one major version (as was done for `ForthDimension`), or
- a one-off rename tool per study (`MuscleBidsMigrate[dir]`?) that renames files and check files and rewrites the json
  keys, plus a version gate so `VersionCheck` does not redo all processing.
- the manual (`QMRITools/BIDS Example/BIDS-config.md`, *How the output files are named*), the example configs (their
  `Target` lists) and `BidsCheckFiles` updated together.
